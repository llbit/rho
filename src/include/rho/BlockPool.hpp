#include <algorithm>
#include <cstdint>
#include <vector>
#include <limits>

#include "rho/GCNode.hpp"

using namespace std;
typedef std::uint64_t u64;

class BlockPool {
    public:
        BlockPool(size_t block_size, size_t superblock_size)
            : m_block_size(block_size) {
                m_bitset_entries = (superblock_size + 63) / 64;
                m_superblock_size = m_bitset_entries * 64;
                m_free = new u64[m_bitset_entries];
                for (int i = 0; i < m_bitset_entries; ++i) {
                    m_free[i] = ~0ull;
                }
                m_storage = new char[m_superblock_size * m_block_size + 15];
                // 16-byte align the superblock pointer:
                m_superblock = reinterpret_cast<char*>(((uintptr_t) m_storage + 15) & (~((uintptr_t) 15)));
                m_block_start = (uintptr_t) m_superblock;
                m_block_end = m_block_start + m_superblock_size * m_block_size;
        }

        ~BlockPool() {
            delete[] m_free;
            delete[] m_storage;
        }

        void* alloc();
        bool try_free(void* p);

        void* get_block_pointer(void* pointer);
        void* apply_to_blocks(std::function<void(rho::GCNode*)> f);
        void print_alloc_stats();

    private:
        uintptr_t m_block_size; // Number of bytes per block.
        uintptr_t m_superblock_size; // Number of blocks in the superblock.
        uintptr_t m_bitset_entries; // Number of entries in the free bitset.
        uintptr_t m_block_start; // Address of first block.
        uintptr_t m_block_end; // Address of one-past last block end.

        unsigned int m_first_bitset = 0;
        unsigned int m_num_victims = 0;
        unsigned int m_last_victim = 0;
        unsigned int m_victim[128];

        u64 volatile* m_free; // Free bitset.
        char* m_storage;
        char* m_superblock;
};

