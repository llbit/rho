#include <algorithm>
#include <cstdint>
#include <vector>
#include <limits>

#include "rho/GCNode.hpp"

using namespace std;
typedef std::uint64_t u64;

class BlockPool {
    public:
        BlockPool(size_t block_size, size_t superblock_size);

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
        unsigned int m_num_victims;
        unsigned int m_last_victim;
        unsigned int m_victim[2048];

        u64 volatile* m_free; // Free bitset.
        char* m_storage;
        char* m_superblock;
};

