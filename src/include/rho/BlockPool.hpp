#include <algorithm>
#include <cstdint>
#include <vector>
#include <limits>

using namespace std;
typedef std::uint64_t u64;

class BlockPool {
    public:
        BlockPool(size_t block_size, size_t superblock_size)
            : m_block_size(block_size),
            m_num_victims(0),
            m_last_victim(0),
            m_next_superblock(-1) {
                m_bitset_entries = (superblock_size + 63) / 64;
                m_superblock_size = m_bitset_entries * 64; // TODO shrink to page size?
        }

        ~BlockPool() {
            for (auto superblock : m_superblocks) {
                delete[] superblock;
            }
        }

        static void initialize();

        static void* pool_alloc(size_t bytes);
        static void* separate_alloc(size_t bytes);
        static void pool_free(void* p);
        static void applyToAllBlocks(std::function<void(void*)> f);

        /** Find heap allocation start pointer. */
        static void* lookup(void* candidate);

        void apply_to_blocks(std::function<void(void*)> f);

    private:
        struct Superblock {
            unsigned num_free;
            u64 free[]; // Free bitset.
        };
        size_t m_block_size; // Number of bytes per block.
        size_t m_superblock_size; // Number of blocks in the superblock.
        size_t m_bitset_entries; // Number of entries in the free bitset.

        unsigned m_num_victims;
        unsigned m_last_victim;
        unsigned m_victim[64];

        int m_next_superblock;
        vector<Superblock*> m_superblocks;

        void* get_next_block(Superblock* superblock);
        void update_next_superblock();
        Superblock* add_superblock();
        void allocate_block(Superblock* superblock, int block);
        void registerSuperblock(int id);

        void* alloc();
        void free(void* p, unsigned superblock_id);
};

