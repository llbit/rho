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
	if (m_next_superblock >= 0) {
		Superblock* superblock = m_superblocks[m_next_superblock];
		return get_next_block(superblock);
	} else {
		// Allocate new superblock.
		return get_next_block(add_superblock());
	}
}

void BlockPool::free(void* pointer)
{
	size_t block = reinterpret_cast<size_t>(pointer);
	int i = 0;
	for (auto superblock : m_superblocks) {
		size_t superblock_start = reinterpret_cast<size_t>(superblock) + block_offset;
		size_t superblock_end = superblock_start + m_block_size * m_superblock_size;
		if (block >= superblock_start && block < superblock_end) {
			size_t index = (block - superblock_start) / m_block_size;
			size_t bitset = index / 64;
			superblock->free[bitset] |= 1ull << (index & 63);
			if (i < m_next_superblock) {
				m_next_superblock = i;
				superblock->next_free = index;
			}
                        if (index < superblock->next_free) {
                            superblock->next_free = index;
                        }
			return;
		}
		i += 1;
	}
	printf("???");
}

// Allocates a new superblock and updates m_next_superblock.
BlockPool::Superblock* BlockPool::add_superblock()
{
	Superblock* superblock = (Superblock*) new char[alignof(u64*)
		+ m_bitset_entries * 8 + m_superblock_size * m_block_size];
	superblock->next_free = 0;
	for (int i = 0; i < m_bitset_entries; ++i) {
		superblock->free[i] = ~0ull;
	}
	size_t superblock_offset = alignof(u64*) + (m_bitset_entries * 8);
	size_t superblock_start = reinterpret_cast<size_t>(superblock) + superblock_offset;
	size_t superblock_end = superblock_start + m_block_size * m_superblock_size;
	if (superblock_start < m_block_start) {
		m_block_start = superblock_start;
	}
	if (superblock_end > m_block_end) {
		m_block_end = superblock_end;
	}
	m_superblocks.push_back(superblock);
        m_next_superblock = m_superblocks.size() - 1;
	return superblock;
}

// Select next free superblock.
void BlockPool::update_next_superblock()
{
	int next = m_next_superblock + 1;
	while (next < m_superblocks.size()) {
		if (m_superblocks[next]->next_free < m_superblock_size) {
			m_next_superblock = next;
			return;
		}
		next += 1;
	}
	// No free superblock found.
	m_next_superblock = -1;
}

// Tag a block as allocated.
void BlockPool::allocate_block(Superblock* superblock, int block)
{
	int bitset = block / 64;
	superblock->free[bitset] &= ~(1ull << (block & 63));
	int next_free;
	do {
		next_free = first_free(superblock->free[bitset]);
		if (next_free >= 0) {
			superblock->next_free = bitset * 64 + next_free;
			return;
		}
		bitset += 1;
	} while (bitset < m_bitset_entries);
	superblock->next_free = m_superblock_size;
        update_next_superblock();
}

void* BlockPool::get_next_block(Superblock* superblock)
{
	size_t block = superblock->next_free;
        allocate_block(superblock, block);
        return reinterpret_cast<char*>(superblock) + alignof(u64*)
            + (m_bitset_entries * 8) + block * m_block_size;
}

void* BlockPool::get_block_pointer(void* pointer)
{
	size_t block = reinterpret_cast<size_t>(pointer);
	if (block < m_block_start || block >= m_block_end) {
		return nullptr;
	}
	for (auto superblock : m_superblocks) {
		size_t superblock_start = reinterpret_cast<size_t>(superblock) + block_offset;
		size_t superblock_end = superblock_start + m_block_size * m_superblock_size;
		if (block >= superblock_start && block < superblock_end) {
			size_t index = (block - superblock_start) / m_block_size;
			size_t bitset = index / 64;
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

void* BlockPool::apply_to_blocks(std::function<void(void*)> f)
{
	for (auto superblock : m_superblocks) {
		size_t block = reinterpret_cast<size_t>(superblock) + block_offset;
		for (int i = 0; i < m_bitset_entries; ++i) {
                    if (superblock->free[i] != ~0ull) {
			for (int index = 0; index < 64; ++index) {
				if (!(superblock->free[i] & (1ull << index))) {
					f(reinterpret_cast<void*>(block));
				}
				block += m_block_size;
			}
                    } else {
                        // Whole block is free -- skip past it.
                        block += 64 * m_block_size;
                    }
		}
	}
	return nullptr;
}

void BlockPool::print_alloc_stats() {
    int superblock_index = 0;
    for (auto superblock : m_superblocks) {
        superblock_index += 1;
        printf(">>>>>>>>>> SUPERBLOCK %d\n", superblock_index);
        for (int i = 0; i < m_bitset_entries; ++i) {
            u64 bitset = superblock->free[i];
            for (int index = 0; index < 64; ++index) {
                if (bitset & (1ull << index)) {
                    printf(".");
                } else {
                    printf("#");
                }
            }
            if (i + 1 < m_bitset_entries) {
                printf(" ");
            }
        }
        printf("\n");
    }
}
