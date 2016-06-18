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

struct SLNode {
    SLNode* pred;
    SLNode* succ;
    size_t size;
    void* data;
};

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

// NUM_POOLS is 1 + the 2 log of the max pool allocation size.
#define NUM_POOLS (13)

static BlockPool* pools[NUM_POOLS];

uintptr_t hash_ptr(uintptr_t ptr)
{
    return (ptr >> LOW_BITS) & (NUM_BUCKET - 1);
}

HashBucket* buckets[NUM_BUCKET];
SparseHashBucket* sparse_buckets[NUM_BUCKET];

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

static SLNode* free_list;

static HashBucket* bucket_from_pointer(void* p);

#ifndef NO_LOG_ALLOCS
FILE* logfile;
#endif

void BlockPool::initialize()
{
#ifndef NO_LOG_ALLOCS
    logfile = fopen("/usr/local/google/home/joqvist/foo.log", "w");
#endif
    free_list = new SLNode();
    free_list->pred = free_list;
    free_list->succ = free_list;
    free_list->data = nullptr;
    free_list->size = 0;

    for (int i = 0; i < NUM_POOLS; ++i) {
        pools[i] = nullptr;
    }

    for (int i = 0; i < NUM_BUCKET; ++i) {
        buckets[i] = nullptr;
        sparse_buckets[i] = nullptr;
    }
}

void add_sparse_block(void* data, size_t size) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(data);
    uintptr_t hash = hash_ptr(pointer);

    SparseHashBucket* bucket = sparse_buckets[hash];
    SparseHashBucket* pred = nullptr;
    while (bucket && data > bucket->data) {
        pred = bucket;
        bucket = bucket->next;
    }
    SparseHashBucket* new_bucket = new SparseHashBucket();
    new_bucket->data = data;
    new_bucket->size = size;
    new_bucket->next = bucket;
    if (pred) {
        pred->next = new_bucket;
    } else {
        sparse_buckets[hash] = new_bucket;
    }
}

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
        delete bucket;
        return true;
    }
    return false;
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

    int log2 = next_log2_16(bytes);
    BlockPool* pool = pools[log2];
    if (!pool) {
        int block_bytes = 1 << log2;
        size_t superblock_size;
        if (block_bytes < 4096) {
            superblock_size = 4096 / block_bytes;
        } else {
            superblock_size = 10;
        }
        pool = new BlockPool(block_bytes, superblock_size);
        pools[log2] = pool;
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
    add_to_allocation_map(result, 1 << log2);
#endif
    update_heap_bounds(result, 1 << log2);
#ifdef ALLOCATION_CHECK
    lookup(result); // Check lookup table consistency.
#endif
    return result;
}

void* BlockPool::separate_alloc(size_t bytes)
{
    void* ptr = new double[(bytes + 7) / 8];
    add_sparse_block(ptr, bytes);
    return ptr;
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
        bucket->pool->free(p);
    } else if (remove_sparse_block(p)) {
        delete[] static_cast<double*>(p);
    } else {
        allocerr("failed to free pointer - unallocated or double-free error");
    }
}
void* BlockPool::alloc()
{
    Superblock* superblock;
    void* block;
    if (m_next_superblock >= 0) {
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

void BlockPool::free(void* pointer)
{
    size_t block = reinterpret_cast<size_t>(pointer);
    int i = 0;
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    for (auto superblock : m_superblocks) {
        size_t superblock_start = reinterpret_cast<size_t>(superblock) + block_offset;
        size_t superblock_end = superblock_start + m_block_size * m_superblock_size;
        if (block >= superblock_start && block < superblock_end) {
            size_t index = (block - superblock_start) / m_block_size;
            size_t bitset = index / 64;
            superblock->num_free += 1;
            superblock->free[bitset] |= 1ull << (index & 63);
            if (i < m_next_superblock) {
                m_next_superblock = i;
            }
            // TODO add to victim buffer.
            return;
        }
        i += 1;
    }
#ifndef NO_LOG_ALLOCS
    fprintf(logfile, "Error: can not free unknown block.\n");
    fflush(logfile);
#endif
    allocerr("can not free unknown block");
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
    if (superblock_start < m_block_start) {
        m_block_start = superblock_start;
    }
    if (superblock_end > m_block_end) {
        m_block_end = superblock_end;
    }
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

void* BlockPool::get_block_pointer(void* pointer)
{
    uintptr_t block = reinterpret_cast<uintptr_t>(pointer);
    if (block < m_block_start || block >= m_block_end) {
        return nullptr;
    }
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    for (auto superblock : m_superblocks) {
        uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + block_offset;
        uintptr_t superblock_end = superblock_start + m_block_size * m_superblock_size;
        if (block >= superblock_start && block < superblock_end) {
            unsigned index = (block - superblock_start) / m_block_size;
            unsigned bitset = index / 64;
            if (!(superblock->free[bitset] & (1ull << (index & 63)))) {
                return reinterpret_cast<char*>(superblock_start) + index * m_block_size;
            } else {
                // The block is not allocated.
                return nullptr;
            }
        }
    }
    return nullptr;
}

void BlockPool::apply_to_blocks(std::function<void(void*)> fun)
{
    uintptr_t block_offset = std::max(sizeof(int), alignof(u64*)) + (m_bitset_entries * 8);
    for (auto superblock : m_superblocks) {
        unsigned blockid = 0;
        uintptr_t block = reinterpret_cast<uintptr_t>(superblock) + block_offset;
        for (int i = 0; i < m_bitset_entries; ++i) {
            u64 bitset = superblock->free[i];
            if (bitset != ~0ull) {
                for (int index = 0; index < 64 && blockid < m_superblock_size; ++index) {
                    if (!(bitset & (1ull << index))) {
#ifdef ALLOCATION_CHECK
                        if (!lookup_in_allocation_map(reinterpret_cast<void*>(block))) {
                            allocerr("apply to all blocks iterating over non-alloc'd pointer");
                        }
#endif
                        fun(reinterpret_cast<void*>(block));
                    }
                    block += m_block_size;
                    blockid += 1;
                }
            } else {
                block += m_block_size * 64;
                blockid += 64;
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

