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

#ifndef RHO_GCNODEALLOCATOR_HPP
#define RHO_GCNODEALLOCATOR_HPP

#include <cstdint>
#include <functional>

#include "rho/AddressSanitizer.hpp"

/** Reports an allocation error and calls abort(). */
void allocerr(const char* message);

namespace rho {
  // Forward declarations.
  class AllocationTable;
  class AllocatorSuperblock;
  class GCNode;

  /**
   * Freelists are used to store superblock nodes and large allocations.
   * For large allocations the block and superblock members are unused.
   */
  struct FreeListNode {
    FreeListNode* m_next;
    std::uint32_t m_block;  // Used only for superblock free nodes.
    AllocatorSuperblock* m_superblock;  // Used only for superblock free nodes.

    /**
     * Insert this allocation into the superblock freelist for the given
     * allocation size.
     *
     * @param size_class the size class for the allocation.
     */
    void addToSuperblockFreelist(unsigned size_class);
  };

  /**
   * \brief An allocator for allocating and iterating GCNode objects of
   * different sizes.
   *
   * This allocator should only be used to allocate GCNodes, otherwise
   * you will trick the garbage collector to treat non-GCNodes as GCNodes
   * during the Sweep phase of garbage collection.
   */
  class GCNodeAllocator {
    public:
      GCNodeAllocator() = delete;

      /**
       * \brief Allocate an object of at least the given size.
       *
       * This should only be used to allocate memory for rho::GCNode objects.
       */
      static void* allocate(std::size_t bytes);

      /** \brief Apply function to all current allocations. */
      static void applyToAllAllocations(std::function<void(void*)> f);

      /** \brief Print allocation overview for debugging. */
      static void debugPrint();

      /** \brief Free a previously allocated object. */
      static void free(void* p);

      /** \brief Must be called before any allocations can be made. */
      static void initialize();

      /**
       * \brief Find heap allocation start pointer.
       *
       * This function finds the corresponding allocation for an internal or
       * exact object pointer.
       *
       * If the pointer is inside the bounds of the small object arena, then
       * the corresponding superblock is found via pointer manipulation.
       *
       * If the pointer is inside the heap bounds, then we iterate over all hash
       * collisions for the pointer to find the corresponding allocation.
       */
      static GCNode* lookupPointer(void* candidate);

    private:
      friend class AllocatorSuperblock;
      friend class AllocationTable;
      friend struct FreeListNode;

      /** Allocation table for medium and large allocations. */
      static AllocationTable* s_alloctable;

      /** Number of different object sizes plus one for the small-object superblocks. */
      static constexpr int s_num_small_pools = (256 / 8) + 1;

      /** Number of different object sizes plus one for the medium-object superblocks. */
      static constexpr int s_num_medium_pools = 18;

      /** Pointers to superblocks with available blocks, indexed by size class. */
      static AllocatorSuperblock* s_superblocks[s_num_small_pools + s_num_medium_pools];

      /**
       * Adds an allocation to a freelist.
       *
       * This should only be used to add medium superblocks and large separate
       * allocations to a freelist. Small allocations are managed in a separate
       * freelist.
       */
      static void addToFreelist(uintptr_t data, unsigned size_log2);

      /**
       * Removes an allocation from a freelist, returns nullptr if there was no
       * available free object.
       */
      static void* removeFromFreelist(unsigned size_log2);

#ifdef HAVE_ADDRESS_SANITIZER

      /*
       * When compiling with the address sanitizer, we modify the allocation routines
       * to detect errors more reliably.  In particular:
       * - A redzone is added on both sides of the object to detect underflows and
       *   overflows.
       * - Memory is poisoned when it is freed to detect use-after-free errors.
       * - Freed objects are held in a quarantine to prevent the memory from
       *   being reallocated quickly.  This helps detect use-after-free errors.
       *
       * Note that in the context of GC, 'use after free' is really premature garbage
       * collection.
       */

      /** Maximum quarantine size = 10Mb. */
      static constexpr int s_max_quarantine_size = 10 << 20;

      /**
       * The quarantine must fit all small superblock size classes, plus
       * all large allocation sizes (= 33 + 64).
       */
      static constexpr int s_num_quarantine_entries = s_num_small_pools + 64;
      static FreeListNode* s_quarantine[s_num_quarantine_entries];

      /** Redzones are added before and after the allocation. */
      static constexpr int s_redzone_size = 16;

      /** Helper function for redzone address calculations. */
      static void* offsetPointer(void* pointer, std::size_t bytes);

      /**
       * Adds an object to the quarantine freelist.
       *
       * The size class is identical to a superlbock size class, or the 2-log
       * of the allocation size plus s_num_small_pools.
       */
      static void addToQuarantine(FreeListNode* free_node, unsigned size_class);

#endif // HAVE_ADDRESS_SANITIZER

  };
}

#endif // RHO_GCNODEALLOCATOR_HPP

