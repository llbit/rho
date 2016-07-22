/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2007   The R Development Core Team.
 *  Copyright (C) 2008-2014  Andrew R. Runnalls.
 *  Copyright (C) 2014 and onwards the Rho Project Authors.
 *
 *  Rho is not part of the R project, and bugs and other issues should
 *  not be reported via r-bugs or other R project channels; instead refer
 *  to the Rho website.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <assert.h>
#include <stdio.h>

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>

#include "rho/AddressSanitizer.hpp"
#include "rho/AllocationTable.hpp"
#include "rho/GCNodeAllocator.hpp"

#ifdef HAVE_ADDRESS_SANITIZER
// Quarantine free lists are used to store freed objects for a while before
// they can be reused. Allocations are poisoned while in the quarantine. This
// makes it more likely that Address Sanitizer can detect use-after-free
// errors.
rho::FreeListNode* rho::GCNodeAllocator::s_quarantine[s_num_small_pools + 64];
#endif

// The allocation table tracks medium and large allocations.
rho::AllocationTable* rho::GCNodeAllocator::s_alloctable = nullptr;

// Statically declared superblock freelists.
rho::AllocatorSuperblock* rho::GCNodeAllocator::s_superblocks[s_num_small_pools + s_num_medium_pools];

// Tracking heap bounds for fast rejection in pointer lookup.
uintptr_t rho::GCNodeAllocator::s_heap_start = UINTPTR_MAX;
uintptr_t rho::GCNodeAllocator::s_heap_end = 0;

// Free lists head pointers.
rho::FreeListNode* rho::GCNodeAllocator::s_freelists[s_num_freelists];

#ifdef ALLOCATION_CHECK
// Helper function for allocator consistency checking.
// An additional allocation map is added which shadows the state of the
// allocator. The extra allocation map is checked to verify each operation
// on the allocator.

static void add_to_allocation_map(void* allocation, std::size_t size);
static void remove_from_allocation_map(void* allocation);

/**
 * This should always return the same result as GCNodeAllocator::lookupPointer().
 * If it does not, something is wrong.
 * When allocation checking is enabled, each call to lookupPointer() is checked
 * against the result of this function.
 * This is also used to test if we try to free something
 * unalloated or allocate something which is already allocated.
 */
static void* lookup_in_allocation_map(void* tentative_pointer);
#endif

void allocerr(const char* message) {
  fprintf(stderr, "ERROR: %s\n", message);
  abort();
}

/**
 * Helper function for finding the the 2-log of the next power of two greater
 * than or equal to size, such that (1 << result) >= size.
 *
 * Handles 32-bit integers.
 */
static unsigned next_log2_32(unsigned size) {
  if (size == 0) {
    return 0;
  }

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

  // First we figure out the 2-log of the highest set bit.
  // Use builtin to count leading zeroes if available.
#if __has_builtin(__builtin_clz)
  unsigned log2 = 31 - __builtin_clz(size);
#else
  // Use manual log-2 finding implementation. The builtin version is available
  // for most compilers, but if not we fall back on this version.
  unsigned log2 = 0;
  // This holds remaining bits after shifting out uninteresting bits:
  unsigned remaining = size;
  if (remaining & 0xFFFF0000) {
    // The high 16 bits have at least one bit set.
    log2 += 16;
    remaining >>= 16;
  }
  if (remaining & 0xFF00) {
    // The high 8 bits of remaining bits have at least one bit set.
    log2 += 8;
    remaining >>= 8;
  }
  if (remaining & 0xF0) {
    // The high 4 bits of remaining bits have at least one bit set.
    log2 += 4;
    remaining >>= 4;
  }
  // Test the upper 2 bits of the remaining bits.
  if (remaining & 0xC) {
    log2 += 2;
    remaining >>= 2;
  }
  // Adjust for the final bits.
  if (remaining > 0) {
    log2 += remaining - 1;
  }
#endif // __has_builtin(__builtin_clz)

  // Now make sure we go one above the highest set bit if the input is larger
  // than that power of two:
  if ((size & (1 << log2)) && (size & ((1 << log2) - 1))) {
    log2 += 1;
  }
  return log2;
}

/**
 * Allocate the small object arena, null out freelist head pointers,
 * and initialize the hashtable.
 */
void rho::GCNodeAllocator::initialize() {
  AllocatorSuperblock::allocateArena();

#ifdef HAVE_ADDRESS_SANITIZER
  // Initialize quarantine freelists.
  for (int i = 0; i < s_num_freelists; ++i) {
    s_quarantine[i] = nullptr;
  }
#endif

  for (int i = 0; i < s_num_small_pools + s_num_medium_pools; ++i) {
    s_superblocks[i] = nullptr;
  }

  for (int i = 0; i < s_num_freelists; ++i) {
    s_freelists[i] = nullptr;
  }

  s_alloctable = new rho::AllocationTable(16);
}

void* rho::GCNodeAllocator::allocate(size_t bytes) {
  void* result = nullptr;
#ifdef HAVE_ADDRESS_SANITIZER
  // Increase size to include redzones.
  bytes += 2 * s_redzone_size;
#endif
  unsigned block_bytes;  // The actual allocation size (may be more than requested).
  if (bytes <= 256) {
    int size_class = (bytes + 7) / 8; // Computes ceil(bytes / 8).
    if (size_class < 4) {
      // Ensure at least 32-byte blocks for the small-object arena. This is
      // required both to fit a FreeListNode (20 bytes), and to reduce the
      // number of bytes needed for the fixed size bitset (as part of the
      // constant SUPERBLOCK_HEADER_SIZE).
      size_class = 4;
    }
    block_bytes = size_class * 8;
    result = removeFromFreelist(size_class);
    if (!result) {
      result = AllocatorSuperblock::allocateBlock(block_bytes);
    }
    if (!result) {
      // Allocating in the small object arena failed, so we have to fall
      // back on using the medium block allocator. To make this work
      // the object size must be at least 64 bytes:
      if (bytes < 64) {
        bytes = 64;
      }
    }
  }
  if (!result) {
    // Default to separate allocation if block size is larger than small block
    // threshold.
    // These allocations are rounded up to the next power of two size.
    unsigned size_log2 = next_log2_32(bytes);
    block_bytes = 1 << size_log2;
    unsigned size_class = size_log2 + s_num_small_pools;
    result = removeFromFreelist(size_class);
    if (!result) {
      if (size_log2 < s_num_medium_pools) {
        result = AllocatorSuperblock::allocateLarge(size_log2);
      } else {
        result = new double[(1L << size_log2) / sizeof(double)];
        GCNodeAllocator::s_alloctable->insert(result, size_log2);
      }
      // Only update heap bounds if allocating a new block.
      updateHeapBounds(result, block_bytes);
    }
  }
  if (!result) {
    allocerr("failed to allocate object");
  }

#ifdef ALLOCATION_CHECK
  if (lookup_in_allocation_map(result)) {
    allocerr("reusing live allocation");
  }
  add_to_allocation_map(result, block_bytes);
  lookupPointer(result);  // Check lookup table consistency.
#endif

#ifdef HAVE_ADDRESS_SANITIZER
  // Poison the redzones.
  // The end redzone can be larger than s_redzone_size if the allocator
  // grew the allocation larger than the requested size.
  // The additional size (block_bytes - bytes) is added to the default
  // redzone size.
  void* end_redzone = offsetPointer(result, bytes - s_redzone_size);
  unsigned end_redzone_size = s_redzone_size + (block_bytes - bytes);
  ASAN_POISON_MEMORY_REGION(result, s_redzone_size);
  ASAN_POISON_MEMORY_REGION(end_redzone, end_redzone_size);

  // Offset the result pointer past first redzone.
  result = offsetPointer(result, s_redzone_size);
#endif // HAVE_ADDRESS_SANITIZER
  return result;
}

void rho::GCNodeAllocator::free(void* pointer) {
#ifdef HAVE_ADDRESS_SANITIZER
  // Adjust for redzone.
  pointer = offsetPointer(pointer, -s_redzone_size);

  // Unpoison the start of the allocation so we can link it into a freelist.
  ASAN_UNPOISON_MEMORY_REGION(pointer, sizeof(FreeListNode));
#endif // HAVE_ADDRESS_SANITIZER
#ifdef ALLOCATION_CHECK
  if (!lookupPointer(pointer)) {
    allocerr("can not free unknown/already-freed pointer");
  }
  remove_from_allocation_map(pointer);
#endif
  uintptr_t pointer_uint = reinterpret_cast<uintptr_t>(pointer);
  AllocatorSuperblock* superblock =
      AllocatorSuperblock::arenaSuperblockFromPointer(pointer_uint);
  if (superblock) {
    superblock->freeBlock(pointer);
  } else if (!s_alloctable->freeAllocation(pointer_uint)) {
    allocerr("failed to free pointer - unallocated or double-free problem");
  }
}

void rho::GCNodeAllocator::applyToAllAllocations(std::function<void(void*)> fun) {
  AllocatorSuperblock::applyToArenaAllocations(fun);
  s_alloctable->applyToAllAllocations(fun);
}

rho::GCNode* rho::GCNodeAllocator::lookupPointer(void* candidate) {
  uintptr_t candidate_uint = reinterpret_cast<uintptr_t>(candidate);
  void* result = AllocatorSuperblock::lookupAllocation(candidate_uint);
  if (!result && (candidate_uint >= s_heap_start && candidate_uint < s_heap_end)) {
    result = s_alloctable->lookup_pointer(candidate_uint);
  }
#ifdef ALLOCATION_CHECK
  if (result != lookup_in_allocation_map(candidate)) {
    allocerr("allocation map mismatch");
  }
#endif
#ifdef HAVE_ADDRESS_SANITIZER
  // Adjust for redzone.
  if (result) {
    result = offsetPointer(result, s_redzone_size);
  }
#endif
  return static_cast<GCNode*>(result);
}

void rho::GCNodeAllocator::addToFreelist(uintptr_t data, unsigned size_class) {
  FreeListNode* free_node = reinterpret_cast<FreeListNode*>(data);
  free_node->m_superblock = nullptr;
  addToFreelist(free_node, size_class);
}

void rho::GCNodeAllocator::addToFreelist(FreeListNode* free_node, unsigned size_class) {
#ifdef HAVE_ADDRESS_SANITIZER
  addToQuarantine(free_node, size_class);
#else
  free_node->m_next = s_freelists[size_class];
  s_freelists[size_class] = free_node;
#endif
}

void* rho::GCNodeAllocator::removeFromFreelist(unsigned size_class) {
  FreeListNode* node = s_freelists[size_class];
  if (node) {
#ifdef HAVE_ADDRESS_SANITIZER
    ASAN_UNPOISON_MEMORY_REGION(node, bytesFromSizeClass(size_class));
#endif
    s_freelists[size_class] = node->m_next;

    if (node->m_superblock) {
      // Tag the resued allocation as allocated since it is part of a superblock.
      node->m_superblock->tagBlockAllocated(node->m_block);
    } else {
      // This is a non-superblock allocation. Insert it back into the allocation table.
      s_alloctable->insert(static_cast<void*>(node),
          size_class - s_num_small_pools);
    }
  }
  return static_cast<void*>(node);
}

void rho::GCNodeAllocator::updateHeapBounds(void* pointer, size_t size) {
  uintptr_t allocation = reinterpret_cast<uintptr_t>(pointer);
  if (allocation < s_heap_start) {
    s_heap_start = allocation;
  }
  uintptr_t allocation_end = allocation + size;
  if (allocation_end > s_heap_end) {
    s_heap_end = allocation_end;
  }
}

void rho::GCNodeAllocator::printSummary() {
  AllocatorSuperblock::debugPrintSmallSuperblocks();
  s_alloctable->printSummary();
}

#ifdef ALLOCATION_CHECK
/*
 * Defining ALLOCATION_CHECK adds a separate allocation-map to
 * track all current allocations. This map is used to double-check
 * each operation the GCNodeAllocator does to ensure it is consistent
 * with the state recorded in the ALLOCATION_CHECK map.
 */

typedef std::map<void*, void*> allocation_map;

/** Extra map for allocator sanity checking. */
static allocation_map allocations;

/** Adds an allocation to the extra allocation map. */
static void add_to_allocation_map(void* allocation, std::size_t size) {
  void* allocation_end = static_cast<char*>(allocation) + size;
  allocations[allocation] = allocation_end;
}

/** Removes an allocation from the extra allocation map. */
static void remove_from_allocation_map(void* allocation) {
  allocations.erase(allocation);
}

/**
 * Check if a pointer is in the extra allocation map.
 * Returns the start pointer of the allocation if the pointer
 * is in the map, otherwise returns null.
 */
static void* lookup_in_allocation_map(void* tentative_pointer) {
  // Find the largest key less than or equal to tentative_pointer.
  allocation_map::const_iterator next_allocation =
      allocations.upper_bound(tentative_pointer);
  if (next_allocation != allocations.begin()) {
    allocation_map::const_iterator allocation = std::prev(next_allocation);

    // Check that tentative_pointer is before the end of the allocation.
    void* allocation_end = allocation->second;
    if (tentative_pointer < allocation_end) {
      return allocation->first;
    }
  }
  return nullptr;
}
#endif

#ifdef HAVE_ADDRESS_SANITIZER
void* rho::GCNodeAllocator::offsetPointer(void* pointer, std::size_t bytes) {
    return static_cast<char*>(pointer) + bytes;
}

void rho::GCNodeAllocator::addToQuarantine(rho::FreeListNode* free_node, unsigned size_class) {
  // Add an allocation to the quarantine and poison it.
  static unsigned quarantine_size = 0;  // Separate from small object quarantine size.

  size_t block_size = bytesFromSizeClass(size_class);
  quarantine_size += block_size;
  free_node->m_next = s_quarantine[size_class];
  s_quarantine[size_class] = free_node;
  ASAN_POISON_MEMORY_REGION(free_node, block_size);

  // Clear the quarantine if it has grown too large.
  if (quarantine_size > s_max_quarantine_size) {
    quarantine_size = 0;
    // Now we clear the quarantine.
    // The quarantine contains both blocks in superblocks and large allocations.
    // Remove all free nodes from quarantine and insert in corresponding freelist.
    for (int i = 0; i < s_num_freelists; ++i) {
      while (s_quarantine[i]) {
        FreeListNode* quarantined_node = s_quarantine[i];

        // Unpoison the quarantined allocation so we can link it into a freelist.
        ASAN_UNPOISON_MEMORY_REGION(quarantined_node, sizeof(FreeListNode));
        s_quarantine[i] = quarantined_node->m_next;
        quarantined_node->m_next = s_freelists[size_class];
        s_freelists[size_class] = quarantined_node;

        // Re-poison so the freelist nodes stay poisoned while free.
        ASAN_POISON_MEMORY_REGION(quarantined_node, sizeof(FreeListNode));
      }
    }
  }
}

#endif // HAVE_ADDRESS_SANITIZER

