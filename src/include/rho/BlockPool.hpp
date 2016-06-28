#include <algorithm>
#include <cstdint>
#include <vector>
#include <limits>

using namespace std;
typedef std::uint64_t u64;
typedef std::uint32_t u32;

class BlockPool {
    public:
        BlockPool(size_t block_size, size_t superblock_bytes);

        static void Initialize();

        static void* AllocBlock(size_t bytes);
        static void FreeBlock(void* p);

        /** Apply function to all allocated blocks (small + large). */
        static void ApplyToAllBlocks(std::function<void(void*)> f);

        /** Find heap allocation start pointer. */
        static void* Lookup(void* candidate);

        /** Print allocation overview for debugging. */
        static void DebugPrint();

        /** Visualize this block pool. */
        void DebugPrintPool();

    private:
        /**
         * Used as an arena for allocating blocks. The Superblock members are
         * only the header part of the structure, the rest is block-size dependent.
         */
        struct Superblock {
            u32 block_size;
            u32 index;
            BlockPool* pool;
            u64 free[]; // Free bitset.
        };

        struct FreeListEntry {
            FreeListEntry* next;
            u32 block;
            u32 superblock_index;
        };

        size_t m_block_size; // Number of bytes per block.
        size_t m_superblock_size; // Number of blocks in each superblock.
        size_t m_bitset_entries; // Number of entries in free bitsets.

        FreeListEntry* m_last_freed; // Free list pointer.

        // Index to next untouched block. All other blocks are either allocated or
        // linked up to the free list. This index is used when allocating and the
        // free list is empyt.
        u32 m_next_untouched;
        u32 m_next_superblock;
        vector<Superblock*> m_superblocks;

        /** Allocate a block from this block pool. */
        void* AllocSmall();

        /** Free a block in this block pool. */
        void FreeSmall(void* p, unsigned superblock_id);

        /** Apply a function to all blocks in this pool. */
        void ApplyToPoolBlocks(std::function<void(void*)> f);

        /** Allocate a large block in the sparse block table. */
        static void* AllocLarge(unsigned bytes);

        Superblock* AddSuperblock();

        /** Tag a block as allocated. */
        void* AllocateBlock(Superblock* superblock, int block);

        /** Tag a block as allocated. */
        static Superblock* SuperblockFromPointer(void* pointer);
};

