#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

using namespace std;

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

void* BlockPool::alloc()
{
    if (m_next_free < m_superblock_size) {
        unsigned int block = m_next_free;
        allocate_block(block);
        return m_superblock + block * m_block_size;
    } else {
        // No space left.
        return nullptr;
    }
}

bool BlockPool::try_free(void* pointer)
{
    uintptr_t block = (uintptr_t) pointer;
    if (block >= m_block_start && block < m_block_end) {
        uintptr_t index = (block - m_block_start) / m_block_size;
        uintptr_t bitset = index / 64;
        // TODO: test for double free?
        m_free[bitset] |= 1ull << (index & 63);
        if (index < m_next_free) {
            m_next_free = index;
        }
        return true;
    }
    return false;
}

// Tag a block as allocated and update m_next_free.
void BlockPool::allocate_block(unsigned int block)
{
    int bitset = block / 64;
    m_free[bitset] &= ~(1ull << (block & 63));
    do {
        int first = first_free(m_free[bitset]);
        if (first >= 0) {
            m_next_free = bitset * 64 + first;
            return;
        }
        bitset += 1;
    } while (bitset < m_bitset_entries);
    m_next_free = m_superblock_size;
}

void* BlockPool::get_block_pointer(void* pointer)
{
    uintptr_t block = (uintptr_t) pointer;
    if (block < m_block_start || block >= m_block_end) {
        return nullptr;
    }
    uintptr_t index = (block - m_block_start) / m_block_size;
    uintptr_t bitset = index / 64;
    if (!(m_free[bitset] & (1ull << (index & 63)))) {
        return m_superblock + index * m_block_size;
    } else {
        // The block is not currently allocated.
        return nullptr;
    }
}

void* BlockPool::apply_to_blocks(std::function<void(void*)> f)
{
    char* block = m_superblock;
    for (int i = 0; i < m_bitset_entries; ++i) {
        if (m_free[i] != ~0ull) {
            for (int index = 0; index < 64; ++index) {
                if (!(m_free[i] & (1ull << index))) {
                    f(block);
                }
                block += m_block_size;
            }
        } else {
            // Whole block is free -- skip past it.
            block += 64 * m_block_size;
        }
    }
    return nullptr;
}

void BlockPool::print_alloc_stats() {
    printf(">>>>>>>>>> SUPERBLOCK\n");
    for (int i = 0; i < m_bitset_entries; ++i) {
        int num_free = 0;
        for (int index = 0; index < 64; ++index) {
            if (m_free[i] & (1ull << index)) {
                num_free += 1;
            }
        }
        switch (num_free) {
            case 0:
                printf("##"); // Full.
                break;
            case 64:
                printf(".."); // Empty
                break;
            default:
                printf("%02d", 64 - num_free);
                break;
        }
        if (i + 1 < m_bitset_entries) {
            printf(" ");
        }
    }
    printf("\n");
}
