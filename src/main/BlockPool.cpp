#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>

#include "rho/BlockPool.hpp"
#include "rho/errors.hpp"

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
#if defined __GNUC__ || __has_builtin(__builtin_ctzll)
    return __builtin_ctzll(bitset);
#else
    // This is an inefficient implementation of finding the number of
    // trailing zeroes.  We rely on __builtin_ctzll above for performance.
    // This version can be replaced to improve the case where __builtin_ctzll
    // is not available.
    int log2 = 0;
    if (!(bitset & 0xFFFFFFFF)) {
        log2 = 32;
        bitset >>= 32;
    }
    if (!(bitset & 0xFFFFull)) {
        log2 += 16;
        bitset >>= 16;
    }
    if (!(bitset & 0xFFull)) {
        log2 += 8;
        bitset >>= 8;
    }
    if (!(bitset & 0xFull)) {
        log2 += 4;
        bitset >>= 4;
    }
    if (!(bitset & 0x3ull)) {
        log2 += 2;
        bitset >>= 2;
    }
    return log2 + 1 - (bitset & 1ull);
#endif
}

BlockPool::BlockPool(size_t block_size, size_t superblock_size)
    : m_block_size(block_size) {
        m_bitset_entries = (superblock_size + 63) / 64;
        m_superblock_size = m_bitset_entries * 64;
        m_free = new u64[m_bitset_entries];
        for (int i = 0; i < m_bitset_entries; ++i) {
            m_free[i] = ~0ull;
        }
        // Initialize victim buffer to point to the first part of the superblock.
        m_num_victims = 2048;
        m_last_victim = 2047;
        for (int i = 0; i < 2048; ++i) {
            m_victim[i] = 2047 - i;
        }
        m_storage = new char[m_superblock_size * m_block_size + 15];
        // 16-byte align the superblock pointer:
        m_superblock = reinterpret_cast<char*>(((uintptr_t) m_storage + 15) & (~((uintptr_t) 15)));
        m_block_start = (uintptr_t) m_superblock;
        m_block_end = m_block_start + m_superblock_size * m_block_size;
}

void* BlockPool::alloc()
{
    if (m_num_victims > 0) {
        unsigned int block = m_victim[m_last_victim];
        m_num_victims -= 1;
        m_last_victim = (m_last_victim + 2047) & 2047;
        int bitset = block / 64;
        m_free[bitset] &= ~(1ull << (block & 63));
        return m_superblock + block * m_block_size;
    } else {
        while (m_first_bitset < m_bitset_entries) {
            int block = first_free(m_free[m_first_bitset]);
            if (block >= 0) {
                m_free[m_first_bitset] &= ~(1ull << block);
                return m_superblock + (m_first_bitset * 64 + block) * m_block_size;
            }
            m_first_bitset += 1;
        }
    }
    // No space left.
    return nullptr;
}

bool BlockPool::try_free(void* pointer)
{
    uintptr_t block = (uintptr_t) pointer;
    if (block >= m_block_start && block < m_block_end) {
        uintptr_t index = (block - m_block_start) / m_block_size;
        uintptr_t bitset = index / 64;
#if 0
        if (m_free[bitset] & (1ull << (index & 63))) {
            // Double free!
            error("Double free in block pool.");
        }
#endif
        m_free[bitset] |= 1ull << (index & 63);
        m_last_victim = (m_last_victim + 1) & 2047;
        m_victim[m_last_victim] = index;
        if (m_num_victims < 2048) {
            m_num_victims += 1;
        }
        if (bitset < m_first_bitset) {
            m_first_bitset = bitset;
        }
        return true;
    }
    return false;
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

void* BlockPool::apply_to_blocks(std::function<void(rho::GCNode*)> f)
{
    char* block = m_superblock;
    for (int i = 0; i < m_bitset_entries; ++i) {
        if (m_free[i] != ~0ull) {
            for (int index = 0; index < 64; ++index) {
                if (!(m_free[i] & (1ull << index))) {
                    f(reinterpret_cast<rho::GCNode*>(block));
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
    printf(">>>>>>>>>> SUPERBLOCK (blocksize=%zu, num block=%zu)\n", m_block_size, m_superblock_size);
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
