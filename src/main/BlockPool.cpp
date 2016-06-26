#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

#define NO_LOG_ALLOCS

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
    HashBucket* next;
    void* data;
    size_t size;
};

#define BUCKET_BITS (9)
#define NUM_FREE_BUCKETS (1 << BUCKET_BITS)

unsigned sparse_bits = 16;
unsigned num_sparse_buckets = 1 << sparse_bits;
unsigned hash_mask = num_sparse_buckets - 1;

#define MAX_COLLISIONS 10

// NUM_POOLS = 1 + (max block size / 8) = 1 + 256 / 8 = 33.
#define NUM_POOLS (33)

static BlockPool* pools[NUM_POOLS];

static unsigned hash_ptr(uintptr_t ptr, unsigned hash_mask) {
    unsigned low = ptr & 0xFFFFFFFF;
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

HashBucket** sparse_buckets = nullptr;
HashBucket* free_set[NUM_FREE_BUCKETS];

HashBucket* deleted_bucket = reinterpret_cast<HashBucket*>(-1);

unsigned add_collision = 0;
unsigned remove_collision = 0;
unsigned lookup_collision = 0;

unsigned sparse_size() {
    unsigned size = 0;
    for (int i = 0; i < num_sparse_buckets; ++i) {
        HashBucket* bucket = sparse_buckets[i];
        if (bucket && bucket != deleted_bucket) {
            size += 1;
        }
    }
    return size;
}

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
    while (!(bitset & 1ull)) {
        log2 += 1;
        bitset >>= 1;
    }
    return log2;
#endif
}

#ifndef NO_LOG_ALLOCS
FILE* logfile;
#endif

HashBucket* alloc_sparse_bucket() {
    static HashBucket* buffer = nullptr;
    static int available = 0;

    if (available == 0) {
        available = 300;
        buffer = new HashBucket[available];
    }

    HashBucket* result = buffer;
    buffer += 1;
    available -= 1;
    return result;
}

// Megaarena is 2 gig = 31 bits.
#define MEGASIZE (1 << 30)

#define SB_BITS (18)
#define SB_SIZE (1 << SB_BITS)

// NOTE: We use a fixed superblock header size, regardless of block size. This
// leaves some unused bitset entries for larger block sizes.
#define SB_HEADER_SIZE (1040)

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
#ifndef NO_LOG_ALLOCS
    logfile = fopen("/usr/local/google/home/joqvist/foo.log", "w");
#endif

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

    for (int i = 0; i < NUM_FREE_BUCKETS; ++i) {
        free_set[i] = nullptr;
    }

    sparse_buckets = new HashBucket*[num_sparse_buckets];
    for (int i = 0; i < num_sparse_buckets; ++i) {
        sparse_buckets[i] = nullptr;
    }
}

bool rebalance_sparse_table(unsigned new_sparse_bits) {
    unsigned new_table_size = 1 << new_sparse_bits;
    unsigned new_mask = new_table_size - 1;
    HashBucket** new_table = new HashBucket*[new_table_size];
    for (int i = 0; i < new_table_size; ++i) {
        new_table[i] = nullptr;
    }
    // Build new table.
    for (int i = 0; i < num_sparse_buckets; ++i) {
        HashBucket* bucket = sparse_buckets[i];
        if (bucket && bucket != deleted_bucket) {
            unsigned hash = hash_ptr(reinterpret_cast<uintptr_t>(bucket->data), new_mask);
            bool placed = false;
            for (int j = 0; j < MAX_COLLISIONS; ++j) {
                if (!new_table[hash]) {
                    new_table[hash] = bucket;
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

void add_sparse_block(HashBucket* new_bucket) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(new_bucket->data);
    while (true) {
        unsigned hash = hash_ptr(pointer, hash_mask);
        for (int i = 0; i < MAX_COLLISIONS; ++i) {
            HashBucket* bucket = sparse_buckets[hash];
            if (!bucket || (bucket == deleted_bucket)) {
                sparse_buckets[hash] = new_bucket;
                return;
            }
            add_collision += 1;
            hash = probe_func(hash, i, hash_mask);
        }
        // Too many collisions, rebalance the hash table.
        unsigned new_sparse_bits = sparse_bits;
        do {
            new_sparse_bits += 1;
        } while (!rebalance_sparse_table(new_sparse_bits));
    }
}

void add_free_block(HashBucket* new_bucket);

bool remove_sparse_block(void* data) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(data);
    unsigned hash = hash_ptr(pointer, hash_mask);
    for (int i = 0; i < MAX_COLLISIONS; ++i) {
        HashBucket* bucket = sparse_buckets[hash];
        if (bucket && bucket != deleted_bucket && bucket->data == data) {
            sparse_buckets[hash] = deleted_bucket;
            add_free_block(bucket);
            return true;
        }
        if (!bucket) {
            break;
        }
        remove_collision += 1;
        hash = probe_func(hash, i, hash_mask);
    }
    return false;
}

// Get next bucket of same size.
HashBucket* get_next_size(HashBucket* bucket) {
    //return bucket->next_size;
    return *(reinterpret_cast<HashBucket**>(bucket->data));
}

// Set link to next bucket of same size.
void set_next_size(HashBucket* bucket, HashBucket* next) {
    //bucket->next_size = next;
    *(reinterpret_cast<HashBucket**>(bucket->data)) = next;
}


void add_free_block(HashBucket* new_bucket) {
    uintptr_t hash = new_bucket->size & (NUM_FREE_BUCKETS - 1);

    HashBucket* bucket = free_set[hash];
    HashBucket* pred = nullptr;
    while (bucket && new_bucket->size > bucket->size) {
        pred = bucket;
        bucket = bucket->next;
    }
    if (bucket && bucket->size == new_bucket->size) {
        set_next_size(new_bucket, bucket);
        new_bucket->next = bucket->next;
    } else {
        set_next_size(new_bucket, nullptr);
        new_bucket->next = bucket;
    }
    if (pred) {
        pred->next = new_bucket;
    } else {
        free_set[hash] = new_bucket;
    }
}

HashBucket* remove_free_block(size_t size) {
    uintptr_t hash = size & (NUM_FREE_BUCKETS - 1);
    HashBucket* bucket = free_set[hash];
    HashBucket* pred = nullptr;
    while (bucket && size > bucket->size) {
        pred = bucket;
        bucket = bucket->next;
    }
    if (bucket && size == bucket->size) {
        if (pred) {
            if (get_next_size(bucket)) {
                pred->next = get_next_size(bucket);
            } else {
                pred->next = bucket->next;
            }
        } else {
            if (get_next_size(bucket)) {
                free_set[hash] = get_next_size(bucket);
            } else {
                free_set[hash] = bucket->next;
            }
        }
        if (get_next_size(bucket)) {
            get_next_size(bucket)->next = bucket->next;
        }
        return bucket;
    }
    return nullptr;
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
    if (bytes > 256) {
        // Default to separate allocation if block size is larger than page size.
        void* result = AllocLarge(bytes);
        if (!result) {
            allocerr("returning null from alloc");
        }
#ifndef NO_LOG_ALLOCS
        fprintf(logfile, "alloc %zu -> %p\n", bytes, result);
#endif
#ifdef ALLOCATION_CHECK
        if (lookup_in_allocation_map(result)) {
            allocerr("reusing live allocation");
        }
        add_to_allocation_map(result, bytes);
#endif
        update_heap_bounds(result, bytes);
#ifdef ALLOCATION_CHECK
        Lookup(result); // Check lookup table consistency.
#endif
        return result;
    }

    int pool_index = (bytes + 7) / 8;
    if (pool_index < 4) {
        // Ensure at least 32-byte blocks. This is required both to to fit a
        // FreeListEntry (16 bytes), and to reduce the number of bytes needed
        // for the fixed size bitset (as part of the constant SB_HEADER_SIZE).
        pool_index = 4;
    }
    int block_bytes = pool_index * 8;
    BlockPool* pool = pools[pool_index];
    if (!pool) {
        pool = new BlockPool(block_bytes, SB_SIZE);
        pools[pool_index] = pool;
    }
    void* result = pool->AllocSmall();
    if (!result) {
        allocerr("returning null from alloc");
    }
#ifndef NO_LOG_ALLOCS
    fprintf(logfile, "alloc %zu -> %p\n", bytes, result);
#endif
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

static void redundant_alloc(size_t bytes) {
    HashBucket* bucket = alloc_sparse_bucket();
    bucket->data = new double[(bytes + 7) / 8];
    bucket->size = bytes;
    add_free_block(bucket);
}

void* BlockPool::AllocLarge(size_t bytes) {
    HashBucket* bucket = remove_free_block(bytes);
    if (bucket) {
        add_sparse_block(bucket);
    } else {
        bucket = alloc_sparse_bucket();
        bucket->data = new double[(bytes + 7) / 8];
        bucket->size = bytes;
        add_sparse_block(bucket);
        redundant_alloc(bytes);
        redundant_alloc(bytes);
    }
    return bucket->data;
}

void BlockPool::FreeBlock(void* p) {
#ifndef NO_LOG_ALLOCS
    fprintf(logfile, "free %p\n", p);
#endif
#ifdef ALLOCATION_CHECK
    if (!Lookup(p)) {
        allocerr("can not free unknown/already-freed pointer");
    }
    remove_from_allocation_map(p);
#endif
    Superblock* superblock = SuperblockFromPointer(p);
    if (superblock) {
        superblock->pool->FreeSmall(p, superblock->index);
    } else if (!remove_sparse_block(p)) {
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
    superblock->free[bitset] |= 1ull << (index & 63);

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
    superblock->free[bitset] &= ~(1ull << (block & 63));
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
                    if (!(bitset & (1ull << index))) {
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
        HashBucket* bucket = sparse_buckets[i];
        if (bucket && bucket != deleted_bucket) {
            fun(bucket->data);
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
        if (!(superblock->free[bitset] & (1ull << (index & 63)))) {
            result = reinterpret_cast<char*>(first_block) + index * superblock->block_size;
        } else {
            // The block is not allocated.
            result = nullptr;
        }
    } else if (candidate >= heap_start && candidate < heap_end) {
        unsigned hash = hash_ptr(candidate_uint, hash_mask);
        for (int i = 0; i < MAX_COLLISIONS; ++i) {
            HashBucket* bucket = sparse_buckets[hash];
            if (bucket && bucket != deleted_bucket && bucket->data <= candidate
                    && (reinterpret_cast<uintptr_t>(bucket->data) + bucket->size) > candidate_uint) {
                result = bucket->data;
                break;
            }
            if (!bucket) {
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


void BlockPool::DebugPrint() {
    for (auto pool : pools) {
        if (pool) {
            pool->DebugPrintPool();
        }
    }

    printf(">>>>> SPARSE TABLE\n");
    for (int i = 0; i < num_sparse_buckets; i += 16) {
        int count = 0;
        for (int j = 0; j < 16 && i + j < num_sparse_buckets; ++j) {
            HashBucket* bucket = sparse_buckets[i + j];
            if (bucket && bucket != deleted_bucket) {
                count += 1;
            }
        }
        if (count) {
            printf("%X", count);
        } else {
            printf(".");
        }
    }
    printf("\n");
}

void BlockPool::DebugPrintPool() {
    printf(">>>>>>>>>> POOL (blocksize=%zu, num block=%zu)\n", m_block_size, m_superblock_size);
    for (int i = 0; i < m_superblocks.size(); ++i) {
        Superblock* superblock = m_superblocks[i];
        printf(">>>>> SUPERBLOCK %d\n", i);
        for (int i = 0; i < m_bitset_entries; ++i) {
            int num_free = 0;
            for (int index = 0; index < 64; ++index) {
                if (superblock->free[i] & (1ull << index)) {
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

void BlockPool::DebugRebalance(int low_bits) {
}

