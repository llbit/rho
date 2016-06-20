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

void* heap_start = reinterpret_cast<void*>(UINTPTR_MAX);
void* heap_end = reinterpret_cast<void*>(0);

struct HashBucket {
    HashBucket* next;
    BlockPool* pool;
    uintptr_t start;
    uintptr_t end;
    unsigned superblock_index;
};

struct SparseHashBucket {
    SparseHashBucket* next;
    void* data;
    size_t size;
};

#define LOW_BITS (14)
#define BUCKET_BITS (9)
#define NUM_BUCKET (1 << 9)

#define SUPERBLOCK_BITS (10)
#define SUPERBLOCK_MASK ((1 << SUPERBLOCK_BITS) - 1)

// NUM_POOLS = (1 + max block size) / 8.
#define NUM_POOLS (513)

static BlockPool* pools[NUM_POOLS];

uintptr_t hash_ptr(uintptr_t ptr)
{
    return (ptr >> LOW_BITS) & (NUM_BUCKET - 1);
}

HashBucket* buckets[NUM_BUCKET];
SparseHashBucket* sparse_buckets[NUM_BUCKET];
SparseHashBucket* free_set[NUM_BUCKET];

inline int first_free(u64 bitset)
{
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

/** Compute next 2-log. */
static int next_log2_16(int size)
{
    if (size == 0) {
        return 0;
    }
    int log2 = 0;
    int temp = size;
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

static HashBucket* bucket_from_pointer(void* p);

#ifndef NO_LOG_ALLOCS
FILE* logfile;
#endif

void BlockPool::initialize()
{
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

static void update_heap_bounds(void* allocation, size_t size)
{
    if (allocation < heap_start) {
        heap_start = allocation;
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    if (allocation_end > heap_end) {
        heap_end = allocation_end;
    }
}

// Find next power of two < 4096.
void* BlockPool::pool_alloc(size_t bytes)
{
    if (bytes > 4096) {
        // Default to separate allocation if block size is larger than page size.
        void* result = separate_alloc(bytes);
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
        lookup(result); // Check lookup table consistency.
#endif
        return result;
    }

    int pool_index = bytes / 8;
    int block_bytes = (pool_index + 1) * 8;
    BlockPool* pool = pools[pool_index];
    if (!pool) {
        size_t superblock_size;
        if (block_bytes < 4096) {
            superblock_size = (4 * 4096) / block_bytes;
        } else {
            superblock_size = 10;
        }
        pool = new BlockPool(block_bytes, superblock_size);
        pools[pool_index] = pool;
    }
    void* result = pool->alloc();
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
    lookup(result); // Check lookup table consistency.
#endif
    return result;
}

void* BlockPool::separate_alloc(size_t bytes)
{
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

void BlockPool::pool_free(void* p)
{
#ifndef NO_LOG_ALLOCS
    fprintf(logfile, "free %p\n", p);
#endif
#ifdef ALLOCATION_CHECK
    if (!lookup(p)) {
        allocerr("can not free unknown/already-freed pointer");
    }
    remove_from_allocation_map(p);
#endif
    HashBucket* bucket = bucket_from_pointer(p);
    if (bucket) {
        bucket->pool->free(p, bucket->superblock_index);
    } else if (remove_sparse_block(p)) {
        // delete[] static_cast<double*>(p); // TODO: cleanup.
    } else {
        allocerr("failed to free pointer - unallocated or double-free error");
    }
}
void* BlockPool::alloc()
{
    Superblock* superblock;
    void* block;
    if (m_num_victims > 0) {
        unsigned block = m_victim[m_last_victim];
        m_num_victims -= 1;
        m_last_victim = (m_last_victim + 63) & 63;
        unsigned index = block >> SUPERBLOCK_BITS;
        unsigned bitset = index / 64;
        superblock = m_superblocks[block & SUPERBLOCK_MASK];
        superblock->num_free -= 1;
        superblock->free[bitset] &= ~(1ull << (index & 63));
        // Update next superblock.
        while (m_next_superblock >= 0) {
            if (m_superblocks[m_next_superblock]->num_free > 0) {
                break;
            } else {
                m_next_superblock += 1;
                if (m_next_superblock == m_superblocks.size()) {
                    m_next_superblock = -1;
                    break;
                }
            }
        }
        return reinterpret_cast<char*>(superblock)
            + std::max(sizeof(int), alignof(u64*))
            + (m_bitset_entries * 8) + index * m_block_size;
    } else if (m_next_superblock >= 0) {
        superblock = m_superblocks[m_next_superblock];
        block = get_next_block(superblock);
        update_next_superblock();
    } else {
        // Allocate new superblock.
        superblock = add_superblock();
        m_next_superblock = m_superblocks.size() - 1;
        block = get_next_block(superblock);
    }
    return block;
}

// Free a pointer inside a given superblock. The block MUST be in the given superblock.
void BlockPool::free(void* pointer, unsigned superblock_id)
{
    size_t block = reinterpret_cast<size_t>(pointer);
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    Superblock* superblock = m_superblocks[superblock_id];
    size_t superblock_start = reinterpret_cast<size_t>(superblock) + block_offset;
    size_t index = (block - superblock_start) / m_block_size;
    size_t bitset = index / 64;
    superblock->num_free += 1;
    superblock->free[bitset] |= 1ull << (index & 63);
    if (superblock_id < m_next_superblock) {
        m_next_superblock = superblock_id;
    }
    // Add to victim buffer only if superblock id fits in the superblock bits.
    if (superblock_id <= SUPERBLOCK_MASK) {
        m_last_victim = (m_last_victim + 1) & 63;
        m_victim[m_last_victim] = (index << SUPERBLOCK_BITS) | superblock_id;
        if (m_num_victims < 64) {
            m_num_victims += 1;
        }
    }
}

/**
 * Inserts new superblock in hash table.
 */
void BlockPool::registerSuperblock(int id) {
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

HashBucket* bucket_from_pointer(void* p)
{
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

BlockPool::Superblock* BlockPool::add_superblock()
{
    Superblock* superblock = (Superblock*) new char[std::max(sizeof(int), alignof(u64*))
        + m_bitset_entries * 8 + m_superblock_size * m_block_size];
    superblock->num_free = m_superblock_size;
    for (int i = 0; i < m_bitset_entries; ++i) {
        superblock->free[i] = ~0ull;
    }
    uintptr_t superblock_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + superblock_offset;
    uintptr_t superblock_end = superblock_start + m_block_size * m_superblock_size;
    m_superblocks.push_back(superblock);
    registerSuperblock(m_superblocks.size() - 1);
    return superblock;
}

void BlockPool::update_next_superblock()
{
    int next = m_next_superblock;
    do {
        if (m_superblocks[next]->num_free > 0) {
            m_next_superblock = next;
            return;
        }
        next += 1;
    } while (next < m_superblocks.size());
    // No free superblock found.
    m_next_superblock = -1;
}

// Tag a block as allocated.
void BlockPool::allocate_block(Superblock* superblock, int block)
{
    int bitset = block / 64;
    superblock->num_free -= 1;
    superblock->free[bitset] &= ~(1ull << (block & 63));
}

void* BlockPool::get_next_block(Superblock* superblock)
{
    if (superblock->num_free > 0) {
        int block;
        unsigned bitset = 0;
        do {
            block = first_free(superblock->free[bitset]);
            if (block >= 0) {
                block += bitset * 64;
                break;
            }
            bitset += 1;
        } while (bitset < m_bitset_entries);
        allocate_block(superblock, block);
        return reinterpret_cast<char*>(superblock) + std::max(sizeof(int), alignof(u64*))
            + (m_bitset_entries * 8) + block * m_block_size;
    }
    return nullptr;
}

void BlockPool::apply_to_blocks(std::function<void(void*)> fun)
{
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

void BlockPool::applyToAllBlocks(std::function<void(void*)> fun)
{
#ifdef ALLOCATION_CHECK
    iterating = true;
#endif
    for (BlockPool* pool : pools) {
        if (pool) {
            pool->apply_to_blocks(fun);
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

void* BlockPool::lookup(void* candidate)
{
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
void add_to_allocation_map(void* allocation, size_t size)
{
    if (iterating) {
        allocerr("allocating node during GC pass");
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    allocations[allocation] = allocation_end;
}

void remove_from_allocation_map(void* allocation)
{
    if (iterating) {
        allocerr("freeing node during GC pass");
    }
    allocations.erase(allocation);
}

void* lookup_in_allocation_map(void* tentative_pointer)
{
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

