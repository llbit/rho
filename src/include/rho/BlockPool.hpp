#include <algorithm>
#include <cstdint>
#include <vector>
#include <limits>

using namespace std;
typedef std::uint64_t u64;

class BlockPool {
	private:
		struct Superblock {
			unsigned int next_free;
			volatile u64 free[]; // Free bitset.
		};
		size_t m_block_size; // Number of bytes per block.
		size_t m_superblock_size; // Number of blocks in the superblock.
		size_t m_bitset_entries; // Number of entries in the free bitset.
		size_t m_block_start; // Address of first block.
		size_t m_block_end; // Address of one-past last block end.
		int m_next_superblock;
		vector<Superblock*> m_superblocks;

                // Offset from start of superblock to first block.
                const size_t block_offset;

		void* get_next_block(Superblock* superblock);
		void update_next_superblock();
		Superblock* add_superblock();
		void allocate_block(Superblock* superblock, int block);
	public:
		BlockPool(size_t block_size, size_t superblock_size)
			: m_block_size(block_size),
			m_block_start(std::numeric_limits<size_t>::max()),
			m_block_end(0),
			m_next_superblock(-1),
                        block_offset(std::max(sizeof(unsigned int), alignof(u64*)) + (((superblock_size + 63) / 64) * sizeof(u64))) {
				m_bitset_entries = (superblock_size + 63) / 64;
				m_superblock_size = m_bitset_entries * 64;
                }

		~BlockPool() {
			for (auto superblock : m_superblocks) {
				delete[] superblock;
			}
		}

		void* alloc();
		void free(void* p);

		void* get_block_pointer(void* pointer);
		void* apply_to_blocks(std::function<void(void*)> f);
};

