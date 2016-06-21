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

        static void Initialize();

        static void* AllocBlock(size_t bytes);
        static void FreeBlock(void* p);

        /** Apply function to all allocated blocks (small + large). */
        static void ApplyToAllBlocks(std::function<void(void*)> f);

        /** Find heap allocation start pointer. */
        static void* Lookup(void* candidate);

    private:
        /**
         * Used as an arena for allocating blocks. The Superblock members are
         * only the header part of the structure, the rest is block-size dependent.
         */
        struct Superblock {
            unsigned num_free;
            u64 free[]; // Free bitset.
        };
        size_t m_block_size; // Number of bytes per block.
        size_t m_superblock_size; // Number of blocks in each superblock.
        size_t m_bitset_entries; // Number of entries in free bitsets.

        unsigned m_num_victims; // Current number of blocks in the victim buffer.
        unsigned m_last_victim; // Index to latest victim in victim buffer.
        unsigned m_victim[64]; // Fixed size FIFO victim buffer.
        // TODO: adjust victim buffer size based on profile/benchmark results.

        int m_next_superblock; // Index to next superblock with a free block.
        vector<Superblock*> m_superblocks;

        /** Allocate a block from this block pool. */
        void* AllocSmall();

        /** Free a block in this block pool. */
        void FreeSmall(void* p, unsigned superblock_id);

        /** Apply a function to all blocks in this pool. */
        void ApplyToPoolBlocks(std::function<void(void*)> f);

        /** Allocate a large block in the sparse block table. */
        static void* AllocLarge(size_t bytes);

        Superblock* AddSuperblock();

        /** Get next free block from a superblock. */
        void* GetNextBlock(Superblock* superblock);

        /** Tag a block as allocated. */
        void AllocateBlock(Superblock* superblock, int block);

        /** Update the index to the next superblock with a free block. */
        void UpdateNextSuperblock();

        /** Add a new superblock to the block pool hash table. */
        void RegisterSuperblock(int id);
};

