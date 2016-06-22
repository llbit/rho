#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>

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

// Hash bucket for the dense block table.
struct HashBucket {
    HashBucket* next;
    BlockPool* pool;
    uintptr_t start;
    uintptr_t end;
    unsigned superblock_index;
};

// Hash bucket for the sparse block table and free lists.
struct SparseHashBucket {
    SparseHashBucket* next;
    void* data;
    size_t size;
};

#define LOW_BITS (17)
#define BUCKET_BITS (9)
#define NUM_BUCKET (1 << 9)

// NUM_POOLS = 1 + (max block size / 8) = 1 + 256 / 8 = 33.
#define NUM_POOLS (33)

static BlockPool* pools[NUM_POOLS];

uintptr_t hash_ptr(uintptr_t ptr) {
    return (ptr >> LOW_BITS) & (NUM_BUCKET - 1);
}

HashBucket* buckets[NUM_BUCKET];
SparseHashBucket* sparse_buckets[NUM_BUCKET];
SparseHashBucket* free_set[NUM_BUCKET];

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

static HashBucket* bucket_from_pointer(void* p);

#ifndef NO_LOG_ALLOCS
FILE* logfile;
#endif

void BlockPool::Initialize() {
#ifndef NO_LOG_ALLOCS
    logfile = fopen("/usr/local/google/home/joqvist/foo.log", "w");
#endif

    for (int i = 0; i < NUM_POOLS; ++i) {
        pools[i] = nullptr;
    }

    for (int i = 0; i < NUM_BUCKET; ++i) {
        buckets[i] = nullptr;
        sparse_buckets[i] = nullptr;
        free_set[i] = nullptr;
    }
}

void add_sparse_block(SparseHashBucket* new_bucket) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(new_bucket->data);
    uintptr_t hash = hash_ptr(pointer);

    SparseHashBucket* bucket = sparse_buckets[hash];
    SparseHashBucket* pred = nullptr;
    while (bucket && new_bucket->data > bucket->data) {
        pred = bucket;
        bucket = bucket->next;
    }
    new_bucket->next = bucket;
    if (pred) {
        pred->next = new_bucket;
    } else {
        sparse_buckets[hash] = new_bucket;
    }
}

void add_free_block(SparseHashBucket* new_bucket);

bool remove_sparse_block(void* data) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(data);
    uintptr_t hash = hash_ptr(pointer);
    SparseHashBucket* bucket = sparse_buckets[hash];
    SparseHashBucket* pred = nullptr;
    while (bucket && data > bucket->data) {
        pred = bucket;
        bucket = bucket->next;
    }
    if (bucket && data == bucket->data) {
        if (pred) {
            pred->next = bucket->next;
        } else {
            sparse_buckets[hash] = bucket->next;
        }
        add_free_block(bucket);
        return true;
    }
    return false;
}

void add_free_block(SparseHashBucket* new_bucket) {
    uintptr_t hash = new_bucket->size & (NUM_BUCKET - 1);

    SparseHashBucket* bucket = free_set[hash];
    SparseHashBucket* pred = nullptr;
    while (bucket && new_bucket->size >= bucket->size) {
        pred = bucket;
        bucket = bucket->next;
    }
    new_bucket->next = bucket;
    if (pred) {
        pred->next = new_bucket;
    } else {
        free_set[hash] = new_bucket;
    }
}

SparseHashBucket* remove_free_block(size_t size) {
    uintptr_t hash = size & (NUM_BUCKET - 1);
    SparseHashBucket* bucket = free_set[hash];
    SparseHashBucket* pred = nullptr;
    while (bucket && size > bucket->size) {
        pred = bucket;
        bucket = bucket->next;
    }
    if (bucket && size == bucket->size) {
        if (pred) {
            pred->next = bucket->next;
        } else {
            free_set[hash] = bucket->next;
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
    if (pool_index < 2) {
        // Ensure at least 16-byte blocks (to fit a FreeListEntry).
        pool_index = 2;
    }
    int block_bytes = pool_index * 8;
    BlockPool* pool = pools[pool_index];
    if (!pool) {
        size_t superblock_size;
        superblock_size = (32 * 4096) / block_bytes;
        pool = new BlockPool(block_bytes, superblock_size);
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
    update_heap_bounds(result, block_bytes);
#ifdef ALLOCATION_CHECK
    Lookup(result); // Check lookup table consistency.
#endif
    return result;
}

void* BlockPool::AllocLarge(size_t bytes) {
    SparseHashBucket* bucket = remove_free_block(bytes);
    if (bucket) {
        add_sparse_block(bucket);
    } else {
        bucket = new SparseHashBucket();
        bucket->data = new double[(bytes + 7) / 8];
        bucket->size = bytes;
        add_sparse_block(bucket);
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
    HashBucket* bucket = bucket_from_pointer(p);
    if (bucket) {
        bucket->pool->FreeSmall(p, bucket->superblock_index);
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
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    Superblock* superblock = m_superblocks[superblock_index];
    size_t superblock_start = reinterpret_cast<size_t>(superblock) + block_offset;
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

/**
 * Inserts new superblock in hash table.
 */
void BlockPool::RegisterSuperblock(int id) {
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    Superblock* superblock = m_superblocks[id];
    uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + block_offset;
    uintptr_t superblock_end = superblock_start + m_block_size * m_superblock_size;
    uintptr_t hash = hash_ptr(superblock_start);
    uintptr_t hash_end = hash_ptr(superblock_end - 1);
    while (true) {
        HashBucket* bucket = buckets[hash];
        HashBucket* pred = nullptr;
        while (bucket && superblock_start > bucket->start) {
            pred = bucket;
            bucket = bucket->next;
        }
        HashBucket* new_bucket = new HashBucket();
        new_bucket->pool = this;
        new_bucket->start = superblock_start;
        new_bucket->end = superblock_end;
        new_bucket->superblock_index = id;
        new_bucket->next = bucket;
        if (!pred) {
            buckets[hash] = new_bucket;
        } else {
            pred->next = new_bucket;
        }
        if (hash == hash_end) {
            break;
        }
        hash = (hash + 1) % NUM_BUCKET;
    }
}

HashBucket* bucket_from_pointer(void* p) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(p);
    uintptr_t hash = hash_ptr(pointer);
    HashBucket* bucket = buckets[hash];
    HashBucket* next = bucket;
    while (next && pointer >= next->start) {
        bucket = next;
        next = next->next;
    }
    if (bucket && pointer >= bucket->start && pointer < bucket->end) {
        return bucket;
    }
    return nullptr;
}

BlockPool::Superblock* BlockPool::AddSuperblock() {
    Superblock* superblock = (Superblock*) new char[std::max(sizeof(int), alignof(u64*))
        + m_bitset_entries * 8 + m_superblock_size * m_block_size];
    for (int i = 0; i < m_bitset_entries; ++i) {
        superblock->free[i] = ~0ull;
    }
    uintptr_t superblock_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + superblock_offset;
    uintptr_t superblock_end = superblock_start + m_block_size * m_superblock_size;
    m_superblocks.push_back(superblock);
    RegisterSuperblock(m_superblocks.size() - 1);
    return superblock;
}

// Tag a block as allocated.
void* BlockPool::AllocateBlock(Superblock* superblock, int block) {
    int bitset = block / 64;
    superblock->free[bitset] &= ~(1ull << (block & 63));
    return reinterpret_cast<char*>(superblock)
            + std::max(sizeof(int), alignof(u64*))
            + (m_bitset_entries * 8) + block * m_block_size;
}

void BlockPool::ApplyToPoolBlocks(std::function<void(void*)> fun) {
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    for (auto superblock : m_superblocks) {
        uintptr_t block = reinterpret_cast<uintptr_t>(superblock) + block_offset;
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
    for (int i = 0; i < NUM_BUCKET; ++i) {
        SparseHashBucket* bucket = sparse_buckets[i];
        while (bucket) {
            fun(bucket->data);
            bucket = bucket->next;
        }
    }
#ifdef ALLOCATION_CHECK
    iterating = false;
#endif
}

void* BlockPool::Lookup(void* candidate) {
    void* result = nullptr;
    if (candidate >= heap_start && candidate < heap_end) {
        HashBucket* bucket = bucket_from_pointer(candidate);
        if (bucket) {
            BlockPool* pool = bucket->pool;
            Superblock* superblock = pool->m_superblocks[bucket->superblock_index];

            unsigned index = (reinterpret_cast<uintptr_t>(candidate) - bucket->start) / pool->m_block_size;
            unsigned bitset = index / 64;
            if (!(superblock->free[bitset] & (1ull << (index & 63)))) {
                result = reinterpret_cast<char*>(bucket->start) + index * pool->m_block_size;
            } else {
                // The block is not allocated.
                result = nullptr;
            }
        } else {
            uintptr_t hash = hash_ptr(reinterpret_cast<uintptr_t>(candidate));
            SparseHashBucket* bucket = sparse_buckets[hash];
            SparseHashBucket* next = bucket;
            while (next && candidate >= next->data) {
                bucket = next;
                next = next->next;
            }
            if (bucket && bucket->data <= candidate && (reinterpret_cast<uintptr_t>(bucket->data) + bucket->size) > reinterpret_cast<uintptr_t>(candidate)) {
                result = bucket->data;
            }
            // Block not found in separate allocation list.
            //return nullptr;
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


void BlockPool::DebugPrint() {
    for (auto pool : pools) {
        if (pool) {
            pool->DebugPrintPool();
        }
    }
    printf(">>>>>>>>>> DENSE TABLE\n");
    unsigned total_size = 0;
    for (int i = 0; i < NUM_BUCKET; i += 20) {
        printf("%03d:", i);
        for (int j = 0; j < 20 && j + i < NUM_BUCKET; ++j) {
            int bucket_size = 0;
            HashBucket* bucket = buckets[i + j];
            while (bucket) {
                bucket_size += 1;
                bucket = bucket->next;
            }
            printf(" %*d", 2, bucket_size);
            total_size += bucket_size;
        }
        printf("\n");
    }
    printf("size: %d\n", total_size);
    printf(">>>>>>>>>> SPARSE TABLE\n");
    total_size = 0;
    for (int i = 0; i < NUM_BUCKET; i += 20) {
        printf("%03d:", i);
        for (int j = 0; j < 20 && j + i < NUM_BUCKET; ++j) {
            int bucket_size = 0;
            SparseHashBucket* bucket = sparse_buckets[i + j];
            while (bucket) {
                bucket_size += 1;
                bucket = bucket->next;
            }
            printf(" %*d", 2, bucket_size);
            total_size += bucket_size;
        }
        printf("\n");
    }
    printf("size: %d\n", total_size);
    printf(">>>>>>>>>> FREE TABLE\n");
    total_size = 0;
    for (int i = 0; i < NUM_BUCKET; i += 20) {
        printf("%03d:", i);
        for (int j = 0; j < 20 && j + i < NUM_BUCKET; ++j) {
            int bucket_size = 0;
            SparseHashBucket* bucket = free_set[i + j];
            while (bucket) {
                bucket_size += 1;
                bucket = bucket->next;
            }
            printf(" %*d", 2, bucket_size);
            total_size += bucket_size;
        }
        printf("\n");
    }
    printf("size: %d\n", total_size);
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
    printf(">>>>>>>>>> DENSE TABLE (Rebalanced)\n");
    unsigned rebalanced[NUM_BUCKET];
    unsigned total_size = 0;
    for (int i = 0; i < NUM_BUCKET; i += 1) {
        rebalanced[i] = 0;
    }
    for (auto pool : pools) {
        if (pool) {
            for (auto superblock : pool->m_superblocks) {
                uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (pool->m_bitset_entries * 8);
                uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + block_offset;
                uintptr_t superblock_end = superblock_start + pool->m_block_size * pool->m_superblock_size;
                unsigned hash = (superblock_start >> low_bits) & (NUM_BUCKET - 1);
                uintptr_t hash_end = ((superblock_end - 1) >> low_bits) & (NUM_BUCKET - 1);
                while (true) {
                    rebalanced[hash] += 1;
                    if (hash == hash_end) {
                        break;
                    }
                    hash = (hash + 1) % NUM_BUCKET;
                }
            }
        }
    }
    for (int i = 0; i < NUM_BUCKET; i += 20) {
        printf("%03d:", i);
        for (int j = 0; j < 20 && j + i < NUM_BUCKET; ++j) {
            printf(" %*d", 2, rebalanced[i + j]);
            total_size += rebalanced[i + j];
        }
        printf("\n");
    }
    printf("size: %d\n", total_size);
    printf(">>>>>>>>>> SPARSE TABLE (Rebalanced)\n");
    total_size = 0;
    for (int i = 0; i < NUM_BUCKET; i += 1) {
        rebalanced[i] = 0;
    }
    for (int i = 0; i < NUM_BUCKET; i += 1) {
        SparseHashBucket* bucket = sparse_buckets[i];
        while (bucket) {
            unsigned hash = (reinterpret_cast<uintptr_t>(bucket->data) >> low_bits) & (NUM_BUCKET - 1);
            rebalanced[hash] += 1;
            bucket = bucket->next;
        }
    }
    for (int i = 0; i < NUM_BUCKET; i += 20) {
        printf("%03d:", i);
        for (int j = 0; j < 20 && j + i < NUM_BUCKET; ++j) {
            printf(" %*d", 2, rebalanced[i + j]);
            total_size += rebalanced[i + j];
        }
        printf("\n");
    }
    printf("size: %d\n", total_size);
}

