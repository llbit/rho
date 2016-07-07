#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

#ifdef ALLOCATION_CHECK
typedef std::map<void*, void*> allocation_map;
static allocation_map allocations;

static bool iterating = false;

static void add_to_allocation_map(void* allocation, size_t size);
static void remove_from_allocation_map(void* allocation);
static void* lookup_in_allocation_map(void* tentative_pointer);
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

// Hash bucket for the sparse block table and free lists.
struct HashBucket {
    uintptr_t data;
    unsigned size;
};

struct FreeNode {
    FreeNode* next;
};

struct SuperblockFreeNode {
    SuperblockFreeNode* next;
    unsigned block;
    uintptr_t superblock;
};


struct BigSuperblock {
    u32 block_size;
    u32 next_untouched;
    SuperblockFreeNode* free_list;
    u64 free[]; // Free bitset.
};

unsigned sparse_bits = 16;
unsigned num_sparse_buckets = 1 << sparse_bits;
unsigned hash_mask = num_sparse_buckets - 1;

// The low pointer bits are masked out in the hash function.
#define LOW_BITS (15)
#define MAX_COLLISIONS (6)
#define MAX_REBALANCE_COLLISIONS (3)

// NUM_POOLS = 1 + (max block size / 8) = 1 + 256 / 8 = 33.
#define NUM_POOLS (33)

static BlockPool* pools[NUM_POOLS];

static unsigned hash_ptr(uintptr_t ptr, unsigned hash_mask) {
    unsigned low = (ptr << 4) & 0xFFFFFFFF;
    unsigned hi = (ptr >> 32) & 0xFFFFFFFF;
    unsigned hash = low ^ hi;
    low = hash & 0xFFFF;
    hi = (hash >> 16) & 0xFFFF;
    hash = low ^ hi ^ (ptr & (~0xFFFF));
    return hash & hash_mask;
}

static unsigned probe_func(unsigned hash, int i, unsigned hash_mask) {
    return (hash + i + 1) & hash_mask;
}

HashBucket* sparse_buckets = nullptr;
FreeNode* free_set[64];
BigSuperblock* large_superblocks[18];

uintptr_t deleted_bucket = UINTPTR_MAX;

unsigned add_collision = 0;
unsigned remove_collision = 0;
unsigned lookup_collision = 0;

// Allocate buckets in batches to get better cache locality when iterating over buckets.
static HashBucket* alloc_sparse_bucket();

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
#define MEDIUM_ARENA_SIZE_LOG2 (19)
#define MEDIUM_ARENA_SIZE (1 << MEDIUM_ARENA_SIZE_LOG2)

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
        free_set[i] = nullptr;
    }

    for (int i = 0; i < 18; ++i) {
        large_superblocks[i] = nullptr;
    }

    sparse_buckets = new HashBucket[num_sparse_buckets];
    for (int i = 0; i < num_sparse_buckets; ++i) {
        sparse_buckets[i].data = 0;
    }
}

bool rebalance_sparse_table(unsigned new_sparse_bits) {
    if (new_sparse_bits > 29) {
        allocerr("sparse table is too large");
    }
    unsigned new_table_size = 1 << new_sparse_bits;
    unsigned new_mask = new_table_size - 1;
    HashBucket* new_table = new HashBucket[new_table_size];
    for (int i = 0; i < new_table_size; ++i) {
        new_table[i].data = 0;
    }
    // Build new table.
    for (int i = 0; i < num_sparse_buckets; ++i) {
        HashBucket& bucket = sparse_buckets[i];
        if (bucket.data != 0 && bucket.data != deleted_bucket) {
            unsigned hash = hash_ptr(bucket.data >> LOW_BITS, new_mask);
            bool placed = false;
            for (int j = 0; j < MAX_REBALANCE_COLLISIONS; ++j) {
                if (new_table[hash].data) {
                    new_table[hash].data = bucket.data;
                    new_table[hash].size = bucket.size;
                    placed = true;
                    break;
                }
                hash = probe_func(hash, j, new_mask);
            }
            if (!placed) {
                delete new_table;
                return false;
            }
        }
    }
    // Replace the old table.
    delete sparse_buckets;
    sparse_buckets = new_table;
    sparse_bits = new_sparse_bits;
    num_sparse_buckets = new_table_size;
    hash_mask = new_mask;
    return true;
}

void add_sparse_block(uintptr_t pointer, size_t size) {
    uintptr_t key = pointer >> LOW_BITS;
    uintptr_t key_end = (pointer + (1L << size) + ((1L << LOW_BITS) - 1)) >> LOW_BITS;
    bool first = true;
    do {
        unsigned hash = hash_ptr(key, hash_mask);
        bool added = false;
        for (int i = 0; i < MAX_COLLISIONS; ++i) {
            HashBucket& bucket = sparse_buckets[hash];
            if (bucket.data == 0 || (bucket.data == deleted_bucket)) {
                if (first) {
                    // Tag this as the first hash entry.
                    first = false;
                    sparse_buckets[hash].data = pointer | 1;
                } else {
                    sparse_buckets[hash].data = pointer;
                }
                sparse_buckets[hash].size = size;
                added = true;
                break;
            }
            add_collision += 1;
            hash = probe_func(hash, i, hash_mask);
        }
        if (!added) {
            // Too many collisions, rebalance the hash table.
            unsigned new_sparse_bits = sparse_bits;
            do {
                new_sparse_bits += 1;
            } while (!rebalance_sparse_table(new_sparse_bits));
        }
        key += 1;
    } while (key < key_end);
}

void add_free_block(uintptr_t data, unsigned size);

bool remove_sparse_block(uintptr_t pointer) {
    uintptr_t key = pointer >> LOW_BITS;
    unsigned hash = hash_ptr(key, hash_mask);
    unsigned size = 0;
    HashBucket* bucket = nullptr;
    for (int i = 0; i < MAX_COLLISIONS; ++i) {
        bucket = &sparse_buckets[hash];
        if (bucket->data != 0) {
            uintptr_t data = bucket->data & ~u64{3};
            uintptr_t first_block = data + SB_HEADER_SIZE;
            if ((bucket->data & 2) && first_block <= pointer && (data + MEDIUM_ARENA_SIZE) > pointer) {
                BigSuperblock* superblock = reinterpret_cast<BigSuperblock*>(data);
                unsigned index = (pointer - first_block) >> superblock->block_size;
                unsigned bitset = index / 64;
                superblock->free[bitset] |= u64{1} << (index & 63);
                SuperblockFreeNode* node = reinterpret_cast<SuperblockFreeNode*>(pointer);
                node->block = index;
                node->superblock = reinterpret_cast<uintptr_t>(superblock);
                if (large_superblocks[superblock->block_size]) {
                    superblock = large_superblocks[superblock->block_size];
                } else {
                    large_superblocks[superblock->block_size] = superblock;
                }
                node->next = superblock->free_list;
                superblock->free_list = node;
                return true;
            } else if (data == pointer) {
                bucket->data = deleted_bucket;
                size = bucket->size;
                add_free_block(pointer, bucket->size);
                break;
            }
        } else {
            // Free failed.
            return false;
        }
        remove_collision += 1;
        hash = probe_func(hash, i, hash_mask);
    }
    key += 1;
    uintptr_t key_end = (pointer + (1L << size) + ((1L << LOW_BITS) - 1)) >> LOW_BITS;
    while (key < key_end) {
        unsigned hash = hash_ptr(key, hash_mask);
        for (int i = 0; i < MAX_COLLISIONS; ++i) {
            HashBucket& bucket = sparse_buckets[hash];
            if (bucket.data == pointer) {
                bucket.data = deleted_bucket;
                break;
            }
            if (!bucket.data) {
                // Free failed.
                return false;
            }
            remove_collision += 1;
            hash = probe_func(hash, i, hash_mask);
        }
        key += 1;
    }
    return true;
}

void add_free_block(uintptr_t data, unsigned size_log2) {
    FreeNode* new_node = reinterpret_cast<FreeNode*>(reinterpret_cast<void*>(data));
    new_node->next = free_set[size_log2];
    free_set[size_log2] = new_node;
}

void* remove_free_block(unsigned size_log2) {
    FreeNode* node = reinterpret_cast<FreeNode*>(free_set[size_log2]);
    if (node) {
        free_set[size_log2] = node->next;
    }
    return static_cast<void*>(node);
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
    if (log2 <= 17) {
        unsigned num_blocks = (MEDIUM_ARENA_SIZE - SB_HEADER_SIZE) >> log2;
        BigSuperblock* superblock;
        if (large_superblocks[log2]) {
            // Reuse existing superblock for this allocation size.
            superblock = large_superblocks[log2];
        } else {
            // Allocate a midsize superblock.
            void* arena = new double[MEDIUM_ARENA_SIZE];
            superblock = reinterpret_cast<BigSuperblock*>(arena);
            superblock->block_size = log2;
            superblock->free_list = nullptr;
            superblock->next_untouched = 0;
            unsigned bitset_entries = ((1 << (MEDIUM_ARENA_SIZE_LOG2 - log2)) + 63) / 64;
            for (int i = 0; i < bitset_entries; ++i) {
                superblock->free[i] = ~0ull;
            }
            add_sparse_block(reinterpret_cast<uintptr_t>(arena) | 2, MEDIUM_ARENA_SIZE_LOG2);
            large_superblocks[log2] = superblock;
        }
        if (superblock->free_list) {
            SuperblockFreeNode* free_node = superblock->free_list;
            unsigned index = free_node->block;
            unsigned bitset = index / 64;
            BigSuperblock* parent = reinterpret_cast<BigSuperblock*>(free_node->superblock);
            parent->free[bitset] &= ~(u64{1} << (index & 63));
            superblock->free_list = free_node->next;
            if (!superblock->free_list && superblock->next_untouched == num_blocks) {
                large_superblocks[log2] = nullptr;
            }
            result = free_node;
        } else {
            unsigned index = superblock->next_untouched;
            unsigned bitset = index / 64;
            superblock->free[bitset] &= ~(u64{1} << (index & 63));
            result = reinterpret_cast<char*>(superblock) + SB_HEADER_SIZE + (index << log2);
            superblock->next_untouched += 1;
            if (superblock->next_untouched == num_blocks) {
                large_superblocks[log2] = nullptr;
            }
        }
    } else {
        result = remove_free_block(log2);
        if (result) {
            add_sparse_block(reinterpret_cast<uintptr_t>(result), log2);
        } else {
            result = new double[1L << log2];
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
            return AllocateBlock(superblock, index);
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
        allocerr("out of superblock space");
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
    return reinterpret_cast<char*>(superblock)
        + SB_HEADER_SIZE + block * m_block_size;
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
    for (int i = 0; i < num_sparse_buckets; ++i) {
        HashBucket& bucket = sparse_buckets[i];
        if (bucket.data && bucket.data != deleted_bucket) {
            if (bucket.data & 1) {
                void* pointer = reinterpret_cast<void*>(bucket.data & ~uintptr_t{3});
                if (bucket.data & 2) {
                    // This is a large superblock.
                    BigSuperblock* superblock = reinterpret_cast<BigSuperblock*>(pointer);
                    uintptr_t block = reinterpret_cast<uintptr_t>(pointer) + SB_HEADER_SIZE;
                    uintptr_t block_end = reinterpret_cast<uintptr_t>(pointer) + MEDIUM_ARENA_SIZE;
                    unsigned block_size = 1 << superblock->block_size;
                    unsigned bitset_entries = ((1 << (MEDIUM_ARENA_SIZE_LOG2 - superblock->block_size)) + 63) / 64;
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
        unsigned hash = hash_ptr(candidate_uint >> LOW_BITS, hash_mask);
        for (int i = 0; i < MAX_COLLISIONS; ++i) {
            HashBucket& bucket = sparse_buckets[hash];
            if (bucket.data && bucket.data != deleted_bucket) {
                uintptr_t pointer = bucket.data & ~uintptr_t{3};
                if (bucket.data & 2) {
                    uintptr_t first_block = pointer + SB_HEADER_SIZE;
                    if (first_block <= candidate_uint && pointer + MEDIUM_ARENA_SIZE > candidate_uint) {
                        BigSuperblock* superblock = reinterpret_cast<BigSuperblock*>(pointer);
                        unsigned index = (candidate_uint - first_block) >> superblock->block_size;
                        unsigned bitset = index / 64;
                        if (!(superblock->free[bitset] & (u64{1} << (index & 63)))) {
                            result = reinterpret_cast<char*>(first_block) + (index << superblock->block_size);
                        }
                        break;
                    }
                } else if (pointer <= candidate_uint && (pointer + (1L << bucket.size)) > candidate_uint) {
                    result = reinterpret_cast<void*>(pointer);
                    break;
                }
            }
            if (!bucket.data) {
                break;
            }
            lookup_collision += 1;
            hash = probe_func(hash, i, hash_mask);
        }
        // Block not found in separate allocation list.
        //return nullptr;
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

static void print_sparse_table() {
    printf(">>>>> SPARSE TABLE\n");
    int size = 0;
    for (int i = 0; i < num_sparse_buckets; i += 16) {
        int count = 0;
        for (int j = 0; j < 16 && i + j < num_sparse_buckets; ++j) {
            HashBucket& bucket = sparse_buckets[i + j];
            if (bucket.data && bucket.data != deleted_bucket) {
                count += 1;
            }
        }
        if (count) {
            printf("%X", count);
        } else {
            printf(".");
        }
        size += count;
    }
    printf("\n");
    printf("table size: %d\n", size);
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
    rebalance_sparse_table(new_bits);
    print_sparse_table();
}

