#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <map>

#define USE_UNORDERED_MAP
#ifdef USE_UNORDERED_MAP
#include <unordered_map>
#else
#include <sparsehash/dense_hash_map>
#endif

#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

// The low pointer bits are masked out in the hash function.
#define LOW_BITS (20)
#define MAX_COLLISIONS (6)
#define MAX_REBALANCE_COLLISIONS (3)

// NUM_POOLS = 1 + (max block size / 8) = 1 + 256 / 8 = 33.
#define NUM_POOLS (33)

#ifdef ALLOCATION_CHECK

typedef std::map<void*, void*> allocation_map;
static allocation_map allocations;

static bool iterating = false;

static void add_to_allocation_map(void* allocation, size_t size);
static void remove_from_allocation_map(void* allocation);
static void* lookup_in_allocation_map(void* tentative_pointer);
#endif

typedef std::pair<uintptr_t, uintptr_t> keytype;

struct PointerHash {
    std::size_t operator()(const keytype& key) const {
        uintptr_t pointer = key.first >> LOW_BITS;
        unsigned low = pointer & 0xFFFFFFFF;
        unsigned hi = (pointer >> 32) & 0xFFFFFFFF;
        unsigned hash = low ^ hi;
        low = hash & 0xFFFF;
        low = (low >> 12) | (low << 4);
        hi = (hash >> 16) & 0xFFFF;
        hash = low ^ hi ^ (pointer & (~0xFFFF));
        return hash;
    }
};

struct PointerEquality {
    bool operator()(const keytype& a, const keytype& b) const {
        // Returns true if a is inside b or b is inside a:
        return (a.first >= b.first && a.second < b.second)
           || (b.first >= a.first && b.second < a.second);
    }
};

struct Allocation {
    uintptr_t data;
    unsigned size_log2;

    Allocation(): data(0), size_log2(0) { }
    Allocation(uintptr_t data, unsigned size): data(data), size_log2(size) { }
};

#ifdef USE_UNORDERED_MAP
std::unordered_map<keytype, Allocation, PointerHash, PointerEquality> sparse_map;
#else
google::dense_hash_map<keytype, Allocation, PointerHash, PointerEquality> sparse_map;
#endif

static void allocerr(const char* message) {
    fprintf(stderr, "ERROR: %s\n", message);
    abort();
}

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

using namespace std;

// Tracks heap bounds for fast pointer rejection during lookup.
void* heap_start = reinterpret_cast<void*>(UINTPTR_MAX);
void* heap_end = reinterpret_cast<void*>(0);

struct FreeNode {
    FreeNode* next;
};

struct LargeSuperblock {
    u32 block_size;
    u32 next_untouched;
    u64 free[]; // Free bitset.
};

struct SuperblockFreeNode {
    SuperblockFreeNode* next;
    unsigned block;
    LargeSuperblock* superblock;
};

unsigned sparse_bits = 16;
unsigned num_sparse_buckets = 1 << sparse_bits;
unsigned hash_mask = num_sparse_buckets - 1;

static BlockPool* pools[NUM_POOLS];

static unsigned probe_func(unsigned hash, int i, unsigned hash_mask) {
    return (hash + i + 1) & hash_mask;
}

SuperblockFreeNode* free_set_sb[64];
FreeNode* free_set[64];
LargeSuperblock* large_superblocks[18]; // Superblocks with untouched blocks.

uintptr_t deleted_bucket = UINTPTR_MAX;

inline int first_free(u64 bitset) {
    if (bitset == 0) {
        return -1;
    }
    // Use builtin count trailing zeroes if available.
#if __has_builtin(__builtin_ctzll)
    return __builtin_ctzll(bitset);
#else
    // This is an inefficient implementation of finding the number of
    // trailing zeroes.  We rely on __builtin_ctzll above for performance.
    // This version can be replaced to improve the case where __builtin_ctzll
    // is not available.
    int log2 = 0;
    while (!(bitset & u64{1})) {
        log2 += 1;
        bitset >>= 1;
    }
    return log2;
#endif
}

static unsigned next_log2_32(int size) {
    if (size == 0) {
        return 0;
    }
    int log2 = 0;
    int temp = size;
    if (temp & 0xFFFF0000) {
        log2 += 16;
        temp >>= 16;
    }
    if (temp & 0xFF00) {
        log2 += 8;
        temp >>= 8;
    }
    if (temp & 0xF0) {
        log2 += 4;
        temp >>= 4;
    }
    if (temp & 0xC) {
        log2 += 2;
        temp >>= 2;
    }
    if (temp > 0) {
        log2 += temp - 1;
    }
    if ((size & (1 << log2)) && (size & ((1 << log2) - 1))) {
        log2 += 1;
    }
    return log2;
}

// Megaarena is 2 gig = 31 bits.
#define MEGASIZE (1 << 30)

#define SB_BITS (18)
#define SB_SIZE (1 << SB_BITS)

// NOTE: We use a fixed superblock header size, regardless of block size. This
// leaves some unused bitset entries for larger block sizes.
#define SB_HEADER_SIZE (1040)

// 1 Mbyte for medium block arenas.
#define LARGE_ARENA_SIZE_LOG2 (19)
#define LARGE_ARENA_SIZE (1 << LARGE_ARENA_SIZE_LOG2)
#define LARGE_HEADER_SIZE (1040)

void* megaarena = nullptr;
uintptr_t sbstart;
uintptr_t sbend;
uintptr_t sbnext;

BlockPool::Superblock* BlockPool::SuperblockFromPointer(void* p) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(p);
    if (pointer >= sbstart && pointer < sbnext) {
        return reinterpret_cast<Superblock*>(pointer & (~static_cast<uintptr_t>(SB_SIZE - 1)));
    } else {
        return nullptr;
    }
}

BlockPool::BlockPool(size_t block_size, size_t superblock_bytes)
    : m_block_size(block_size),
    m_next_untouched(0),
    m_next_superblock(0),
    m_last_freed(nullptr) {
        m_superblock_size = (superblock_bytes - 1040) / block_size;
        m_bitset_entries = (m_superblock_size + 63) / 64;
}

void BlockPool::Initialize() {
    megaarena = sbrk(MEGASIZE);
    uintptr_t start = (uintptr_t) megaarena;
    uintptr_t end = start + MEGASIZE;
    uintptr_t pad = SB_SIZE - (start & (SB_SIZE - 1));
    uintptr_t num_sb = (end - sbstart) >> SB_BITS;
    sbstart = start + pad;
    sbend = sbstart + num_sb * SB_SIZE;
    sbnext = sbstart;

    for (int i = 0; i < NUM_POOLS; ++i) {
        pools[i] = nullptr;
    }

    for (int i = 0; i < 64; ++i) {
        free_set_sb[i] = nullptr;
        free_set[i] = nullptr;
    }

    for (int i = 0; i < 18; ++i) {
        large_superblocks[i] = nullptr;
    }
}

void add_sparse_block(uintptr_t pointer, unsigned size_log2) {
    uintptr_t data = pointer;
    uintptr_t block_end = pointer + (uintptr_t{1} << size_log2);
    uintptr_t key = pointer >> LOW_BITS;
    uintptr_t key_end = (pointer + (1L << size_log2) + ((1L << LOW_BITS) - 1)) >> LOW_BITS;
    do {
        sparse_map.insert({ keytype(pointer, block_end), { data, size_log2 } });
        key += 1;
        pointer = key << LOW_BITS;
        data &= ~uintptr_t{1}; // Mask away first flag.
    } while (key < key_end);
}

void add_free_block(uintptr_t data, unsigned size);

bool remove_sparse_block(uintptr_t pointer) {
    keytype keyobj(pointer, pointer);
    auto iter = sparse_map.find(keyobj);
    if (iter != sparse_map.end()) {
        uintptr_t data = iter->second.data;
        unsigned size_log2 = iter->second.size_log2;

        if (data & 2) {
            data &= ~uintptr_t{3}; // Mask away flag bits.
            // This is a big superblock.
            LargeSuperblock* superblock = reinterpret_cast<LargeSuperblock*>(data);
            uintptr_t first_block = data + SB_HEADER_SIZE;
            unsigned index = (pointer - first_block) >> superblock->block_size;
            unsigned bitset = index / 64;
            superblock->free[bitset] |= u64{1} << (index & 63);
            SuperblockFreeNode* node = reinterpret_cast<SuperblockFreeNode*>(pointer);
            node->block = index;
            node->superblock = superblock;
            node->next = free_set_sb[superblock->block_size];
            free_set_sb[superblock->block_size] = node;
            return true;
        } else {
            data &= ~uintptr_t{3}; // Mask away flag bits.
            // This is a separate allocation.
            // Remove all internal mapped pointers.
            uintptr_t key = pointer >> LOW_BITS;
            uintptr_t key_end = (pointer + (1L << size_log2) + ((1L << LOW_BITS) - 1)) >> LOW_BITS;
            do {
                sparse_map.erase(keytype(pointer, pointer));
                key += 1;
                pointer = key << LOW_BITS;
            } while (key < key_end);
            add_free_block(data, size_log2);
        }
        return true;
    }
    return false;
}

void add_free_block(uintptr_t pointer, unsigned size_log2) {
    FreeNode* new_node = reinterpret_cast<FreeNode*>(pointer);
    new_node->next = free_set[size_log2];
    free_set[size_log2] = new_node;
}

void* remove_free_block(unsigned size_log2) {
    if (size_log2 <= 17) {
        SuperblockFreeNode* node = free_set_sb[size_log2];
        if (node) {
            unsigned index = node->block;
            unsigned bitset = index / 64;
            node->superblock->free[bitset] &= ~(u64{1} << (index & 63)); // Mark block as allocated.
            free_set_sb[size_log2] = node->next;
        }
        return static_cast<void*>(node);
    } else {
        FreeNode* node = free_set[size_log2];
        if (node) {
            free_set[size_log2] = node->next;
        }
        add_sparse_block(reinterpret_cast<uintptr_t>(node), size_log2);
        return static_cast<void*>(node);
    }
}

static void update_heap_bounds(void* allocation, size_t size) {
    if (allocation < heap_start) {
        heap_start = allocation;
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    if (allocation_end > heap_end) {
        heap_end = allocation_end;
    }
}

void* BlockPool::AllocBlock(size_t bytes) {
    void* result;
    unsigned block_bytes;
    if (bytes <= 256) {
        int pool_index = (bytes + 7) / 8;
        if (pool_index < 4) {
            // Ensure at least 32-byte blocks. This is required both to to fit a
            // FreeListEntry (16 bytes), and to reduce the number of bytes needed
            // for the fixed size bitset (as part of the constant SB_HEADER_SIZE).
            pool_index = 4;
        }
        block_bytes = pool_index * 8;
        BlockPool* pool = pools[pool_index];
        if (!pool) {
            pool = new BlockPool(block_bytes, SB_SIZE);
            pools[pool_index] = pool;
        }
        result = pool->AllocSmall();
        if (!result) {
            allocerr("ran out of space for small allocations");
        }
    } else {
        // Default to separate allocation if block size is larger than small block threshold.
        unsigned log2 = next_log2_32(bytes);
        block_bytes = 1 << log2;
        result = AllocLarge(log2);
        update_heap_bounds(result, block_bytes);
    }
    if (!result) {
        allocerr("returning null from alloc");
    }

#ifdef ALLOCATION_CHECK
    if (lookup_in_allocation_map(result)) {
        allocerr("reusing live allocation");
    }
    add_to_allocation_map(result, block_bytes);
#endif
#ifdef ALLOCATION_CHECK
    Lookup(result); // Check lookup table consistency.
#endif
    return result;
}

void* BlockPool::AllocLarge(unsigned log2) {
    void* result = nullptr;
    result = remove_free_block(log2);
    if (!result) {
        if (log2 <= 17) {
            unsigned num_blocks = (LARGE_ARENA_SIZE - LARGE_HEADER_SIZE) >> log2;
            LargeSuperblock* superblock;
            if (large_superblocks[log2]) {
                // Reuse existing superblock for this allocation size.
                superblock = large_superblocks[log2];
            } else {
                // Allocate a midsize superblock.
                void* arena = new double[LARGE_ARENA_SIZE / sizeof(double)];
                superblock = reinterpret_cast<LargeSuperblock*>(arena);
                superblock->block_size = log2;
                superblock->next_untouched = 0;
                unsigned bitset_entries = ((1 << (LARGE_ARENA_SIZE_LOG2 - log2)) + 63) / 64;
                for (int i = 0; i < bitset_entries; ++i) {
                    superblock->free[i] = ~0ull;
                }
                add_sparse_block(reinterpret_cast<uintptr_t>(arena) | 2, LARGE_ARENA_SIZE_LOG2);
                large_superblocks[log2] = superblock;
            }
            unsigned index = superblock->next_untouched;
            unsigned bitset = index / 64;
            superblock->free[bitset] &= ~(u64{1} << (index & 63));
            result = reinterpret_cast<char*>(superblock) + LARGE_HEADER_SIZE + (index << log2);
            superblock->next_untouched += 1;
            if (superblock->next_untouched == num_blocks) {
                large_superblocks[log2] = nullptr;
            }
        } else {
            result = new double[(1L << log2) / sizeof(double)];
            add_sparse_block(reinterpret_cast<uintptr_t>(result), log2);
        }
    }
    return result;
}

void BlockPool::FreeBlock(void* p) {
#ifdef ALLOCATION_CHECK
    if (!Lookup(p)) {
        allocerr("can not free unknown/already-freed pointer");
    }
    remove_from_allocation_map(p);
#endif
    Superblock* superblock = SuperblockFromPointer(p);
    if (superblock) {
        superblock->pool->FreeSmall(p, superblock->index);
    } else if (!remove_sparse_block(reinterpret_cast<uintptr_t>(p))) {
        allocerr("failed to free pointer - unallocated or double-free problem");
    }
}

void* BlockPool::AllocSmall() {
    if (m_last_freed) {
        u32 index = m_last_freed->block;
        u32 superblock_index = m_last_freed->superblock_index;
        m_last_freed = m_last_freed->next;
        return AllocateBlock(m_superblocks[superblock_index], index);
    } else {
        u32 index = m_next_untouched;
        u32 superblock_index = m_next_superblock;
        if (index + 1 < m_superblock_size) {
            m_next_untouched += 1;
        } else {
            m_next_untouched = 0;
            m_next_superblock += 1;
        }
        if (superblock_index < m_superblocks.size()) {
            return AllocateBlock(m_superblocks[superblock_index], index);
        } else {
            // Allocate new superblock.
            Superblock* superblock = AddSuperblock();
            if (superblock) {
                return AllocateBlock(superblock, index);
            } else {
                // No more room for allocations of this size!
                return nullptr;
            }
        }
    }
}

// Free a pointer inside a given superblock. The block MUST be in the given superblock.
void BlockPool::FreeSmall(void* pointer, unsigned superblock_index) {
    size_t block = reinterpret_cast<size_t>(pointer);
    Superblock* superblock = m_superblocks[superblock_index];
    size_t superblock_start = reinterpret_cast<size_t>(superblock) + SB_HEADER_SIZE;
    size_t index = (block - superblock_start) / m_block_size;
    size_t bitset = index / 64;
    superblock->free[bitset] |= u64{1} << (index & 63);

    // Use the block as a free list node and prepend to the free list.
    FreeListEntry* node = reinterpret_cast<FreeListEntry*>(pointer);
    node->next = m_last_freed;
    node->block = index;
    node->superblock_index = superblock_index;
    m_last_freed = node;
}

BlockPool::Superblock* BlockPool::AddSuperblock() {
    if (sbnext >= sbend) {
        return nullptr;
    }
    Superblock* superblock = reinterpret_cast<Superblock*>(sbnext);
    superblock->block_size = m_block_size;
    superblock->index = m_superblocks.size();
    superblock->pool = this;
    sbnext += SB_SIZE;
    // Here we mark all bitset entries as free, later we don't have to do
    // precise range checking while iterating allocated blocks.
    for (int i = 0; i < m_bitset_entries; ++i) {
        superblock->free[i] = ~0ull;
    }
    m_superblocks.push_back(superblock);
    return superblock;
}

// Tag a block as allocated.
void* BlockPool::AllocateBlock(Superblock* superblock, int block) {
    int bitset = block / 64;
    superblock->free[bitset] &= ~(u64{1} << (block & 63));
    return reinterpret_cast<char*>(superblock) + SB_HEADER_SIZE + block * m_block_size;
}

void BlockPool::ApplyToPoolBlocks(std::function<void(void*)> fun) {
    for (auto superblock : m_superblocks) {
        uintptr_t block = reinterpret_cast<uintptr_t>(superblock) + SB_HEADER_SIZE;
        uintptr_t block_end = reinterpret_cast<uintptr_t>(superblock) + SB_SIZE;
        for (int i = 0; i < m_bitset_entries; ++i) {
            u64 bitset = superblock->free[i];
            if (bitset != ~0ull) {
                for (int index = 0; index < 64; ++index) {
                    if (!(bitset & (u64{1} << index))) {
#ifdef ALLOCATION_CHECK
                        if (!lookup_in_allocation_map(reinterpret_cast<void*>(block))) {
                            allocerr("apply to all blocks iterating over non-alloc'd pointer");
                        }
#endif
                        fun(reinterpret_cast<void*>(block));
                    }
                    block += m_block_size;
                }
            } else {
                block += m_block_size * 64;
            }
        }
    }
}

void BlockPool::ApplyToAllBlocks(std::function<void(void*)> fun) {
#ifdef ALLOCATION_CHECK
    iterating = true;
#endif
    for (BlockPool* pool : pools) {
        if (pool) {
            pool->ApplyToPoolBlocks(fun);
        }
    }
    for (auto iter = sparse_map.begin(); iter != sparse_map.end(); ++iter) {
        uintptr_t data = iter->second.data;
        unsigned size_log2 = iter->second.size_log2;
        if (data & 1) {
            void* pointer = reinterpret_cast<void*>(data & ~uintptr_t{3}); // Mask away flags.
            if (data & 2) {
                // This is a large superblock.
                LargeSuperblock* superblock = reinterpret_cast<LargeSuperblock*>(pointer);
                uintptr_t block = reinterpret_cast<uintptr_t>(pointer) + LARGE_HEADER_SIZE;
                uintptr_t block_end = reinterpret_cast<uintptr_t>(pointer) + LARGE_ARENA_SIZE;
                unsigned block_size = 1 << superblock->block_size;
                unsigned bitset_entries = ((1 << (LARGE_ARENA_SIZE_LOG2 - superblock->block_size)) + 63) / 64;
                for (int i = 0; i < bitset_entries; ++i) {
                    if (superblock->free[i] != ~0ull) {
                        for (int index = 0; index < 64; ++index) {
                            if (!(superblock->free[i] & (u64{1} << index))) {
#ifdef ALLOCATION_CHECK
                                if (!lookup_in_allocation_map(reinterpret_cast<void*>(block))) {
                                    allocerr("apply to all blocks iterating over non-alloc'd pointer");
                                }
#endif
                                fun(reinterpret_cast<void*>(block));
                            }
                            block += block_size;
                        }
                    } else {
                        block += block_size * 64;
                    }
                }
            } else {
                fun(pointer);
            }
        }
    }
#ifdef ALLOCATION_CHECK
    iterating = false;
#endif
}

void* BlockPool::Lookup(void* candidate) {
    void* result = nullptr;
    uintptr_t candidate_uint = reinterpret_cast<uintptr_t>(candidate);
    Superblock* superblock = SuperblockFromPointer(candidate);
    if (superblock) {
        uintptr_t first_block = reinterpret_cast<uintptr_t>(superblock) + SB_HEADER_SIZE;
        if (candidate_uint < first_block) {
            // The candidate pointer points inside the superblock header.
            return nullptr;
        }
        unsigned index = (candidate_uint - first_block) / superblock->block_size;
        unsigned bitset = index / 64;
        if (!(superblock->free[bitset] & (u64{1} << (index & 63)))) {
            result = reinterpret_cast<char*>(first_block) + index * superblock->block_size;
        } else {
            // The block is not allocated.
            result = nullptr;
        }
    } else if (candidate >= heap_start && candidate < heap_end) {
        keytype keyobj(candidate_uint, candidate_uint);
        auto iter = sparse_map.find(keyobj);
        if (iter != sparse_map.end()) {
            uintptr_t data = iter->second.data;
            unsigned size_log2 = iter->second.size_log2;
            uintptr_t pointer = data & ~uintptr_t{3};
            if (data & 2) {
                uintptr_t first_block = pointer + LARGE_HEADER_SIZE;
                if (first_block <= candidate_uint && pointer + LARGE_ARENA_SIZE > candidate_uint) {
                    // TODO make size dynamic.
                    LargeSuperblock* superblock = reinterpret_cast<LargeSuperblock*>(pointer);
                    unsigned index = (candidate_uint - first_block) >> superblock->block_size;
                    unsigned bitset = index / 64;
                    if (!(superblock->free[bitset] & (u64{1} << (index & 63)))) {
                        result = reinterpret_cast<char*>(first_block) + (index << superblock->block_size);
                    }
                }
            } else if (pointer <= candidate_uint && (pointer + (1L << size_log2)) > candidate_uint) {
                result = reinterpret_cast<void*>(pointer);
            }
        }
    }
#ifdef ALLOCATION_CHECK
    if (result != lookup_in_allocation_map(candidate)) {
        allocerr("allocation map mismatch");
    }
#endif
    return result;
}

#ifdef ALLOCATION_CHECK
void add_to_allocation_map(void* allocation, size_t size) {
    if (iterating) {
        allocerr("allocating node during GC pass");
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    allocations[allocation] = allocation_end;
}

void remove_from_allocation_map(void* allocation) {
    if (iterating) {
        allocerr("freeing node during GC pass");
    }
    allocations.erase(allocation);
}

void* lookup_in_allocation_map(void* tentative_pointer) {
    // Find the largest key less than or equal to tentative_pointer.
    allocation_map::const_iterator next_allocation = allocations.upper_bound(tentative_pointer);
    if (next_allocation != allocations.begin()) {
        allocation_map::const_iterator allocation = std::prev(next_allocation);

        // Check that tentative_pointer is before the end of the allocation.
        void* allocation_end = allocation->second;
        if (tentative_pointer < allocation_end) { // Less-than-or-equal handles one-past-end pointers.
            return allocation->first;
        }
    }
    return nullptr;
}
#endif

void print_sparse_table() {
    printf("######## SPARSE TABLE\n");
    for (auto iter = sparse_map.begin(); iter != sparse_map.end(); ++iter) {
        uintptr_t data = iter->second.data;
        unsigned size_log2 = iter->second.size_log2;
        printf("%p .. %p -> %p (%d)\n", reinterpret_cast<void*>(iter->first.first),
                reinterpret_cast<void*>(iter->first.second),
                reinterpret_cast<void*>(data), size_log2);
    }
}

void BlockPool::DebugPrint() {
    for (auto pool : pools) {
        if (pool) {
            pool->DebugPrintPool();
        }
    }
    print_sparse_table();
}

void BlockPool::DebugPrintPool() {
    printf(">>>>>>>>>> POOL (blocksize=%zu, num block=%zu)\n", m_block_size, m_superblock_size);
    for (int i = 0; i < m_superblocks.size(); ++i) {
        Superblock* superblock = m_superblocks[i];
        printf(">>>>> SUPERBLOCK %d\n", i);
        for (int i = 0; i < m_bitset_entries; ++i) {
            int num_free = 0;
            for (int index = 0; index < 64; ++index) {
                if (superblock->free[i] & (u64{1} << index)) {
                    num_free += 1;
                }
            }
            switch (num_free) {
                case 0:
                    printf("#"); // Full.
                    break;
                case 64:
                    printf("."); // Empty
                    break;
                default:
                    printf("%d", (64 - num_free) / 10);
                    break;
            }
            if (i + 1 < m_bitset_entries) {
                printf(" ");
            }
        }
        printf("\n");
    }
}

void BlockPool::DebugRebalance(int new_bits) {
    print_sparse_table();
}

