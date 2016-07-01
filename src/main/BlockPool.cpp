#include "rho/BlockPool.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

#define USE_STD_MAP

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

// Tracks heap bounds for fast pointer rejection during lookup.
void* heap_start = reinterpret_cast<void*>(UINTPTR_MAX);
void* heap_end = reinterpret_cast<void*>(0);

#define RED (1)
#define BLACK (0)

// Red-black tree node. Null pointer represents a black leaf node.
struct TreeNode {
    TreeNode* parent;
    TreeNode* left;
    TreeNode* right;
    uintptr_t data;
    unsigned size;
    unsigned color;
};

static std::map<uintptr_t, unsigned> sparse_allocs;

TreeNode* root = nullptr;

TreeNode* new_tree_node(uintptr_t data, unsigned size) {
    TreeNode* node = new TreeNode();
    node->data = data;
    node->size = size;
    node->left = nullptr;
    node->right = nullptr;
    node->parent = nullptr;
    node->color = RED;
    return node;
}

struct FreeNode {
    FreeNode* next;
};

// NUM_POOLS = 1 + (max block size / 8) = 1 + 256 / 8 = 33.
#define NUM_POOLS (33)

static BlockPool* pools[NUM_POOLS];

FreeNode* free_set[64];

inline int first_free(u64 bitset) {
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

// Find next power of two greater than or equal to the input size.
static unsigned next_log2_32(int size) {
    if (size == 0) {
        return 0;
    }
    int log2 = 0;
    int temp = size;
    if (temp & 0xFFFF0000) {
        log2 += 16;
        temp >>= 16;
    }
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

// Megaarena is 2 gig = 31 bits.
#define MEGASIZE (1 << 30)

#define SB_BITS (18)
#define SB_SIZE (1 << SB_BITS)

// NOTE: We use a fixed superblock header size, regardless of block size. This
// leaves some unused bitset entries for larger block sizes.
#define SB_HEADER_SIZE (1040)

void add_free_block(uintptr_t data, unsigned size);

void* megaarena = nullptr;
uintptr_t sbstart;
uintptr_t sbend;
uintptr_t sbnext;

BlockPool::Superblock* BlockPool::SuperblockFromPointer(void* p) {
    uintptr_t pointer = reinterpret_cast<uintptr_t>(p);
    if (pointer >= sbstart && pointer < sbnext) {
        return reinterpret_cast<Superblock*>(pointer & (~static_cast<uintptr_t>(SB_SIZE - 1)));
    } else {
        return nullptr;
    }
}

BlockPool::BlockPool(size_t block_size, size_t superblock_bytes)
    : m_block_size(block_size),
    m_next_untouched(0),
    m_next_superblock(0),
    m_last_freed(nullptr) {
        m_superblock_size = (superblock_bytes - 1040) / block_size;
        m_bitset_entries = (m_superblock_size + 63) / 64;
}

void BlockPool::Initialize() {
    megaarena = sbrk(MEGASIZE);
    uintptr_t start = (uintptr_t) megaarena;
    uintptr_t end = start + MEGASIZE;
    uintptr_t pad = SB_SIZE - (start & (SB_SIZE - 1));
    uintptr_t num_sb = (end - sbstart) >> SB_BITS;
    sbstart = start + pad;
    sbend = sbstart + num_sb * SB_SIZE;
    sbnext = sbstart;

    for (int i = 0; i < NUM_POOLS; ++i) {
        pools[i] = nullptr;
    }

    for (int i = 0; i < 64; ++i) {
        free_set[i] = nullptr;
    }
}

TreeNode* search_tree(uintptr_t data) {
    TreeNode* node = root;
    while (node) {
        if (node->data <= data && (node->data + (1L << node->size)) > data) {
            return node;
        } else if (data < node->data) {
            node = node->left;
        } else {
            node = node->right;
        }
    }
    return nullptr;
}

void left_rotate(TreeNode* x) {
    TreeNode* y = x->right;
    if (y == nullptr) {
        return;
    }
    y->parent = x->parent;
    if (x->parent) {
        if (x->parent->left == x) {
            x->parent->left = y;
        } else {
            x->parent->right = y;
        }
    } else {
        root = y;
    }
    if (y->left) {
        y->left->parent = x;
    }
    x->right = y->left;
    y->left = x;
    x->parent = y;
}

void right_rotate(TreeNode* y) {
    TreeNode* x = y->left;
    if (x == nullptr) {
        return;
    }
    x->parent = y->parent;
    if (y->parent != nullptr) {
        if (y->parent->left == y) {
            y->parent->left = x;
        } else {
            y->parent->right = x;
        }
    } else {
        root = x;
    }
    if (x->right != nullptr) {
        x->right->parent = y;
    }
    y->left = x->right;
    x->right = y;
    y->parent = x;
}

void insert_fixup(TreeNode* p) {
    while (p->parent != nullptr && p->parent->color == RED) {
	TreeNode* pp = p->parent;
	TreeNode* gp = pp->parent;
	TreeNode* y;
        if (gp == nullptr) {
            break;
        }
	if (pp == gp->left) {
	    y = gp->right;
	    if (y && y->color == RED) {
		// Case 1
		pp->color = BLACK;
		y->color = BLACK;
		gp->color = RED;
		p = gp;
	    } else {
		if (pp->right == p) {
		    // Case 2
		    left_rotate(pp);
		    p = pp;
		    pp = p->parent;
		}
		pp->color = BLACK;
		gp->color = RED;
		right_rotate(gp);
	    }
	} else {
	    y = gp->left;
	    if (y && y->color == RED) {
		// Case 1
		pp->color = BLACK;
		y->color = BLACK;
		gp->color = RED;
		p = gp;
	    } else {
		if (pp->left == p) {
		    // Case 2
		    right_rotate(pp);
		    p = pp;
		    pp = p->parent;
		}
		pp->color = BLACK;
		gp->color = RED;
		left_rotate(gp);
	    }
	}
    }
}

void add_sparse_block(uintptr_t data, unsigned size) {
#ifdef USE_STD_MAP
    sparse_allocs[data] = size;
#else
    if (!root) {
        root = new_tree_node(data, size);
        root->color = BLACK;
        return;
    }

    TreeNode* node = root;
    while (true) {
        if (data < node->data) {
            if (node->left) {
                node = node->left;
            } else {
                node->left = new_tree_node(data, size);
                node->left->parent = node;
                insert_fixup(node->left);
                return;
            }
        } else {
            if (node->right) {
                node = node->right;
            } else {
                node->right = new_tree_node(data, size);
                node->right->parent = node;
                insert_fixup(node->right);
                return;
            }
        }
    }
#endif
}

void delete_fixup(TreeNode* x, TreeNode* pp) {
    while (x != root && (x == nullptr || x->color == BLACK)) {
        if (x == pp->left) {
            TreeNode* w = pp->right;
            if (w != nullptr && w->color == RED) {
                w->color = BLACK;
                pp->color = RED;
                left_rotate(pp);
                w = pp->right;
            }
            if ((w->left == nullptr || w->left->color == BLACK)
                    && (w->right == nullptr || w->right->color == BLACK)) {
                w->color = RED;
                x = pp;
            } else {
                if (w->right == nullptr || w->right->color == BLACK) {
                    w->color = RED;
                    right_rotate(w);
                    w = pp->right;
                }
                w->color = pp->color;
                pp->color = BLACK;
                w->right->color = BLACK;
                left_rotate(pp);
                x = root;
            }
        } else {
            TreeNode* w = pp->left;
            if (w->color == RED) {
                w->color = BLACK;
                pp->color = RED;
                right_rotate(pp);
                w = pp->left;
            }
            if ((w->left == nullptr || w->left->color == BLACK)
                    && (w->right == nullptr || w->right->color == BLACK)) {
                w->color = RED;
                x = pp;
            } else {
                if (w->left == nullptr || w->left->color == BLACK) {
                    w->color = RED;
                    left_rotate(w);
                    w = pp->left;
                }
                w->color = pp->color;
                pp->color = BLACK;
                w->left->color = BLACK;
                right_rotate(pp);
                x = root;
            }
        }
        pp = x->parent;
    }
    if (x) {
        x->color = BLACK;
    }
}

TreeNode* minimum(TreeNode* x) {
    while (x->left != nullptr) {
        x = x->left;
    }
    return x;
}

TreeNode* successor(TreeNode* x) {
    if (x->right != nullptr) {
        return minimum(x->right);
    }
    TreeNode* y = x->parent;
    while (y != nullptr && x == y->right) {
        x = y;
        y = y->parent;
    }
    return y;
}

void delete_rb_node(TreeNode* p) {
    TreeNode* y;
    TreeNode* x;
    if (p->left == nullptr || p->right == nullptr) {
        y = p;
    } else {
        y = successor(p);
    }
    if (y->left != nullptr) {
        x = y->left;
    } else {
        x = y->right;
    }
    TreeNode* parentx = y->parent;
    if (x != nullptr) {
        x->parent = y->parent;
    }
    if (y->parent == nullptr) {
        root = x;
    } else {
        if (y == y->parent->left) {
            y->parent->left = x;
        } else {
            y->parent->right = x;
        }
    }
    if (y != p) {
        p->data = y->data;
        p->size = y->size;
    }
    if (y->color == BLACK) {
        delete_fixup(x, parentx);
    }
    delete y;
}

bool remove_sparse_block(uintptr_t data) {
#ifdef USE_STD_MAP
    std::map<uintptr_t, unsigned>::const_iterator next_allocation =
        sparse_allocs.upper_bound(data);
    if (next_allocation != sparse_allocs.begin()) {
        auto allocation = std::prev(next_allocation);
        add_free_block(allocation->first, allocation->second);
        sparse_allocs.erase(allocation);
        return true;
    }
    return false;
#else
    TreeNode* node = root;
    while (node) {
        if (node->data == data) {
            add_free_block(data, node->size);
            delete_rb_node(node);
            return true;
        } else if (data < node->data) {
            node = node->left;
        } else {
            node = node->right;
        }
    }
    return false;
#endif
}

void add_free_block(uintptr_t data, unsigned log2) {
    FreeNode* new_node = reinterpret_cast<FreeNode*>(data);
    new_node->next = free_set[log2];
    free_set[log2] = new_node;
}

void* remove_free_block(unsigned log2) {
    FreeNode* node = reinterpret_cast<FreeNode*>(free_set[log2]);
    if (node) {
        free_set[log2] = node->next;
    }
    return static_cast<void*>(node);
}

static void update_heap_bounds(void* allocation, size_t size) {
    if (allocation < heap_start) {
        heap_start = allocation;
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    if (allocation_end > heap_end) {
        heap_end = allocation_end;
    }
}

void* BlockPool::AllocBlock(size_t bytes) {
    void* result;
    unsigned block_bytes;
    if (bytes > 256) {
        // Default to separate allocation if block size is larger than page size.
        unsigned log2 = next_log2_32(bytes);
        block_bytes = 1 << log2;
        result = AllocLarge(log2);
        update_heap_bounds(result, block_bytes);
    } else {
        int pool_index = (bytes + 7) / 8;
        if (pool_index < 4) {
            // Ensure at least 32-byte blocks. This is required both to to fit a
            // FreeListEntry (16 bytes), and to reduce the number of bytes needed
            // for the fixed size bitset (as part of the constant SB_HEADER_SIZE).
            pool_index = 4;
        }
        block_bytes = pool_index * 8;
        BlockPool* pool = pools[pool_index];
        if (!pool) {
            pool = new BlockPool(block_bytes, SB_SIZE);
            pools[pool_index] = pool;
        }
        result = pool->AllocSmall();
    }
    if (!result) {
        allocerr("returning null from alloc");
    }

#ifdef ALLOCATION_CHECK
    if (lookup_in_allocation_map(result)) {
        allocerr("reusing live allocation");
    }
    add_to_allocation_map(result, block_bytes);
#endif
#ifdef ALLOCATION_CHECK
    Lookup(result); // Check lookup table consistency.
#endif
    return result;
}

void* BlockPool::AllocLarge(unsigned log2) {
    void* result = remove_free_block(log2);
    if (result) {
        add_sparse_block(reinterpret_cast<uintptr_t>(result), log2);
    } else {
        result = new double[1 << log2];
        add_sparse_block(reinterpret_cast<uintptr_t>(result), log2);
    }
    return result;
}

void BlockPool::FreeBlock(void* p) {
#ifdef ALLOCATION_CHECK
    if (!Lookup(p)) {
        allocerr("can not free unknown/already-freed pointer");
    }
    remove_from_allocation_map(p);
#endif
    Superblock* superblock = SuperblockFromPointer(p);
    if (superblock) {
        superblock->pool->FreeSmall(p, superblock);
    } else if (!remove_sparse_block(reinterpret_cast<uintptr_t>(p))) {
        allocerr("failed to free pointer - unallocated or double-free problem");
    }
}

void* BlockPool::AllocSmall() {
    if (m_last_freed) {
        unsigned index = m_last_freed->block;
        Superblock* superblock = m_last_freed->superblock;
        m_last_freed = m_last_freed->next;
        return AllocateBlock(superblock, index);
    } else {
        unsigned index = m_next_untouched;
        unsigned superblock_index = m_next_superblock;
        if (index + 1 < m_superblock_size) {
            m_next_untouched += 1;
        } else {
            m_next_untouched = 0;
            m_next_superblock += 1;
        }
        if (superblock_index < m_superblocks.size()) {
            return AllocateBlock(m_superblocks[superblock_index], index);
        } else {
            // Allocate new superblock.
            Superblock* superblock = AddSuperblock();
            return AllocateBlock(superblock, index);
        }
    }
}

// Free a pointer inside a given superblock. The block MUST be in the given superblock.
void BlockPool::FreeSmall(void* pointer, Superblock* superblock) {
    uintptr_t block = reinterpret_cast<uintptr_t>(pointer);
    uintptr_t superblock_start = reinterpret_cast<uintptr_t>(superblock) + SB_HEADER_SIZE;
    unsigned index = (block - superblock_start) / m_block_size;
    unsigned bitset = index / 64;
    superblock->free[bitset] |= 1ull << (index & 63);

    // Use the block as a free list node and prepend to the free list.
    FreeListEntry* node = reinterpret_cast<FreeListEntry*>(pointer);
    node->next = m_last_freed;
    node->block = index;
    node->superblock = superblock;
    m_last_freed = node;
}

BlockPool::Superblock* BlockPool::AddSuperblock() {
    if (sbnext >= sbend) {
        allocerr("out of superblock space");
    }
    Superblock* superblock = reinterpret_cast<Superblock*>(sbnext);
    superblock->block_size = m_block_size;
    superblock->index = m_superblocks.size();
    superblock->pool = this;
    sbnext += SB_SIZE;
    // Here we mark all bitset entries as free, later we don't have to do
    // precise range checking while iterating allocated blocks.
    for (int i = 0; i < m_bitset_entries; ++i) {
        superblock->free[i] = ~0ull;
    }
    m_superblocks.push_back(superblock);
    return superblock;
}

// Tag a block as allocated.
void* BlockPool::AllocateBlock(Superblock* superblock, unsigned block) {
    unsigned bitset = block / 64;
    superblock->free[bitset] &= ~(1ull << (block & 63));
    return reinterpret_cast<char*>(superblock)
        + SB_HEADER_SIZE + block * m_block_size;
}

void BlockPool::ApplyToPoolBlocks(std::function<void(void*)> fun) {
    for (auto superblock : m_superblocks) {
        uintptr_t block = reinterpret_cast<uintptr_t>(superblock) + SB_HEADER_SIZE;
        uintptr_t block_end = reinterpret_cast<uintptr_t>(superblock) + SB_SIZE;
        for (int i = 0; i < m_bitset_entries; ++i) {
            u64 bitset = superblock->free[i];
            if (bitset != ~0ull) {
                for (int index = 0; index < 64; ++index) {
                    if (!(bitset & (1ull << index))) {
#ifdef ALLOCATION_CHECK
                        if (!lookup_in_allocation_map(reinterpret_cast<void*>(block))) {
                            allocerr("apply to all blocks iterating over non-alloc'd pointer");
                        }
#endif
                        fun(reinterpret_cast<void*>(block));
                    }
                    block += m_block_size;
                }
            } else {
                block += m_block_size * 64;
            }
        }
    }
}

void apply_to_tree(TreeNode* node, std::function<void(void*)> fun) {
    if (node) {
        fun(reinterpret_cast<void*>(node->data));
        apply_to_tree(node->left, fun);
        apply_to_tree(node->right, fun);
    }
}

void BlockPool::ApplyToAllBlocks(std::function<void(void*)> fun) {
#ifdef ALLOCATION_CHECK
    iterating = true;
#endif
    for (BlockPool* pool : pools) {
        if (pool) {
            pool->ApplyToPoolBlocks(fun);
        }
    }
    apply_to_tree(root, fun);
#ifdef ALLOCATION_CHECK
    iterating = false;
#endif
}

void* BlockPool::Lookup(void* candidate) {
    void* result = nullptr;
    uintptr_t candidate_uint = reinterpret_cast<uintptr_t>(candidate);
    Superblock* superblock = SuperblockFromPointer(candidate);
    if (superblock) {
        uintptr_t first_block = reinterpret_cast<uintptr_t>(superblock) + SB_HEADER_SIZE;
        if (candidate_uint < first_block) {
            // The candidate pointer points inside the superblock header.
            return nullptr;
        }
        unsigned index = (candidate_uint - first_block) / superblock->block_size;
        unsigned bitset = index / 64;
        if (!(superblock->free[bitset] & (1ull << (index & 63)))) {
            result = reinterpret_cast<char*>(first_block) + index * superblock->block_size;
        } else {
            // The block is not allocated.
            result = nullptr;
        }
    } else if (candidate >= heap_start && candidate < heap_end) {
#ifdef USE_STD_MAP
        std::map<uintptr_t, unsigned>::const_iterator next_allocation =
            sparse_allocs.upper_bound(candidate_uint);
        if (next_allocation != sparse_allocs.begin()) {
            auto allocation = std::prev(next_allocation);

            // Check that tentative_pointer is before the end of the allocation.
            unsigned size = allocation->second;
            if (candidate_uint < (allocation->first + (1L << size))) { // Less-than-or-equal handles one-past-end pointers.
                result = reinterpret_cast<void*>(allocation->first);
            }
        }
#else
        TreeNode* node = search_tree(candidate_uint);
        if (node) {
            result = reinterpret_cast<void*>(node->data);
        }
#endif
    }
#ifdef ALLOCATION_CHECK
    if (result != lookup_in_allocation_map(candidate)) {
        allocerr("allocation map mismatch");
    }
#endif
    return result;
}

#ifdef ALLOCATION_CHECK
void add_to_allocation_map(void* allocation, size_t size) {
    if (iterating) {
        allocerr("allocating node during GC pass");
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    allocations[allocation] = allocation_end;
}

void remove_from_allocation_map(void* allocation) {
    if (iterating) {
        allocerr("freeing node during GC pass");
    }
    allocations.erase(allocation);
}

void* lookup_in_allocation_map(void* tentative_pointer) {
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

static void tree_print(TreeNode* node, int indent) {
    printf("%*c", indent, ' ');
    if (node) {
        if (node->color == RED) {
            printf("red: %p\n", reinterpret_cast<void*>(node->data));
        } else {
            printf("blk: %p\n", reinterpret_cast<void*>(node->data));
        }
        tree_print(node->left, indent + 2);
        tree_print(node->right, indent + 2);
    } else {
        printf("nil\n");
    }
}

static void rb_tree_print() {
    tree_print(root, 0);
}

void BlockPool::DebugPrint() {
    for (auto pool : pools) {
        if (pool) {
            pool->DebugPrintPool();
        }
    }

#ifndef USE_STD_MAP
    rb_tree_print();
#endif

    printf("#### FREE SET\n");
    for (int i = 7; i < 64; ++i) {
        unsigned count = 0;
        FreeNode* node = free_set[i];
        while (node) {
            node = node->next;
            count += 1;
        }
        if (count > 0) {
            printf("%d: %d\n", i, count);
        }
    }
}

void BlockPool::DebugPrintPool() {
    printf(">>>>>>>>>> POOL (blocksize=%zu, num block=%zu)\n", m_block_size, m_superblock_size);
    for (int i = 0; i < m_superblocks.size(); ++i) {
        Superblock* superblock = m_superblocks[i];
        printf(">>>>> SUPERBLOCK %d\n", i);
        for (int i = 0; i < m_bitset_entries; ++i) {
            int num_free = 0;
            for (int index = 0; index < 64; ++index) {
                if (superblock->free[i] & (1ull << index)) {
                    num_free += 1;
                }
            }
            switch (num_free) {
                case 0:
                    printf("#"); // Full.
                    break;
                case 64:
                    printf("."); // Empty
                    break;
                default:
                    printf("%d", (64 - num_free) / 10);
                    break;
            }
            if (i + 1 < m_bitset_entries) {
                printf(" ");
            }
        }
        printf("\n");
    }
}

