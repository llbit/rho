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

/** @file GCNode.cpp
 *
 * Class GCNode and associated C-callable functions.
 */

#include "rho/GCNode.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>

#include "rho/AddressSanitizer.hpp"
#include "rho/GCManager.hpp"
#include "rho/GCRoot.hpp"
#include "rho/GCStackFrameBoundary.hpp"
#include "rho/ProtectStack.hpp"
#include "rho/RAllocStack.hpp"
#include "rho/WeakRef.hpp"

#include "rho/BlockPool.hpp"

using namespace std;
using namespace rho;

vector<const GCNode*>* GCNode::s_moribund = 0;
unsigned int GCNode::s_num_nodes = 0;
bool GCNode::s_on_stack_bits_correct = false;

static bool iterating = false;
typedef std::map<void*, void*> allocation_map;

#define ALLOCATION_CHECK

#ifdef ALLOCATION_CHECK
static allocation_map allocations;

static void add_to_allocation_map(void* allocation, size_t size);
static void remove_from_allocation_map(void* allocation);
static void* lookup_in_allocation_map(void* tentative_pointer);
#endif

void* heap_start = reinterpret_cast<void*>(UINTPTR_MAX);
void* heap_end = reinterpret_cast<void*>(0);

// Used to update reference count bits of a GCNode. The array element at index
// n+1 is XORed with the current rcmms bits to compute the updated reference
// count n+1.  This does not overflow the reference count bits and preserves
// the other bits in m_rcmms.
const unsigned char GCNode::s_decinc_refcount[]
= {0,    2, 2, 6, 6, 2, 2, 0xe, 0xe, 2, 2, 6, 6, 2, 2, 0x1e,
   0x1e, 2, 2, 6, 6, 2, 2, 0xe, 0xe, 2, 2, 6, 6, 2, 2, 0x3e,
   0x3e, 2, 2, 6, 6, 2, 2, 0xe, 0xe, 2, 2, 6, 6, 2, 2, 0x1e,
   0x1e, 2, 2, 6, 6, 2, 2, 0xe, 0xe, 2, 2, 6, 6, 2, 0,    0};

unsigned char GCNode::s_mark = 0;

HOT_FUNCTION void* GCNode::operator new(size_t bytes)
{
    GCManager::maybeGC();
    MemoryBank::notifyAllocation(bytes);
    void *result;

    result = BlockPool::pool_alloc(bytes);
#ifdef ALLOCATION_CHECK
    add_to_allocation_map(result, bytes);
#endif

    // Because garbage collection may occur between this point and the GCNode's
    // constructor running, we need to ensure that this space is at least
    // minimally initialized.
    // We construct a GCNode with a ficticious reference count of 1, to prevent
    // it from getting marked as moribund or garbage collected (the mark-sweep
    // GC uses the reference counts, so this is always effective).
    // It will be overwritten by the real constructor.
    new (result)GCNode(static_cast<CreateAMinimallyInitializedGCNode*>(nullptr));
    return result;
}

GCNode::GCNode(CreateAMinimallyInitializedGCNode*) : m_rcmms(s_decinc_refcount[1]) {
}

void GCNode::operator delete(void* p, size_t bytes)
{
    MemoryBank::notifyDeallocation(bytes);

#ifdef ALLOCATION_CHECK
    remove_from_allocation_map(p);
#endif
    BlockPool::pool_free(p);
}

bool GCNode::check()
{
    // Check moribund list:
    for (const GCNode* node : *s_moribund) {
	if (!(node->m_rcmms & s_moribund_mask)) {
	    cerr << "GCNode::check() : "
		"Node on moribund list without moribund bit set.\n";
	    abort();
	}
    }

    return true;
}

void GCNode::destruct_aux()
{
    // Erase this node from the moribund list:
    typedef std::vector<const GCNode*>::iterator Iter;
    Iter it = std::find(s_moribund->begin(), s_moribund->end(), this);
    if (it == s_moribund->end())
	abort();  // Should never happen!
    s_moribund->erase(it);
}
    
extern RObject* R_Srcref;

void GCNode::gc(bool markSweep)
{
    if (GCManager::GCInhibitor::active())
	return;
    GCManager::GCInhibitor inhibitor;

    ProtectStack::protectAll();
    incRefCount(R_Srcref);

    if (markSweep) {
	GCStackRootBase::withAllStackNodesProtected(markSweepGC);
    } else {
	GCStackRootBase::withAllStackNodesProtected(gclite);
    }

    decRefCount(R_Srcref);
}

void GCNode::markSweepGC()
{
    // NB: setting this flag implies that the garbage collection will ignore
    // any new stack roots.  To ensure correctness, this function must not call
    // any code that depends on normal operation of the garbage collector.
    s_on_stack_bits_correct = true;

    mark();
    sweep();

    s_on_stack_bits_correct = false;
}

void GCNode::gclite()
{
    s_on_stack_bits_correct = true;

    while (!s_moribund->empty()) {
	// Last in, first out, for cache efficiency:
	const GCNode* node = s_moribund->back();
	s_moribund->pop_back();
	// Clear moribund bit.  Beware ~ promotes to unsigned int.
	node->m_rcmms &= static_cast<unsigned char>(~s_moribund_mask);

	if (node->maybeGarbage()) {
	    delete node;
        }
    }

    s_on_stack_bits_correct = false;
}

void GCNode::initialize()
{
    BlockPool::initialize();
    s_moribund = new vector<const GCNode*>();
}

void GCNode::makeMoribund() const
{
    if (s_on_stack_bits_correct) {
	// In this case, the node can be deleted immediately.
	delete this;
    } else {
        addToMoribundList();
    }
}

void GCNode::addToMoribundList() const
{
    m_rcmms |= s_moribund_mask;
    s_moribund->push_back(this);
}

void GCNode::mark()
{
    // In the first mark-sweep collection, the marking of a node is
    // indicated by the mark bit being set; in the second mark sweep,
    // marking is indicated by the bit being clear, and so on in
    // alternation.  This avoids the need for the sweep phase to
    // iterate through the surviving nodes simply to remove marks.
    s_mark ^= s_mark_mask;
    GCNode::Marker marker;
    GCRootBase::visitRoots(&marker);
    GCStackRootBase::visitRoots(&marker);
    ProtectStack::visitRoots(&marker);
    WeakRef::markThru();
    if (R_Srcref)
	marker(R_Srcref);
}

void GCNode::sweep()
{
    // Detach the referents of nodes that haven't been marked.
    // Once this is done, all of the nodes in the cycle will be unreferenced
    // and they will have been deleted unless their reference count is
    // saturated.
    vector<void*> work;
    vector<GCNode*> to_delete;
    iterating = true;
    BlockPool::applyToAllBlocks([&](void* p) {
#ifdef ALLOCATION_CHECK
            if (!lookup_in_allocation_map(p)) {
                error("pointer not in alloc map!");
            }
#endif
            work.push_back(p);
            });
    iterating = false;
    for (void* pointer : work) {
        if (BlockPool::lookup(pointer)) {
#ifdef ALLOCATION_CHECK
            if (!lookup_in_allocation_map(pointer)) {
                error("pointer not in alloc map 2!");
            }
#endif
            // The pointer is still allocated, so detach referents.
            GCNode* node = static_cast<GCNode*>(pointer);
            if (!node->isMarked()) {
                int ref_count = node->getRefCount();
                incRefCount(node);
                if (node->getRefCount() == ref_count) {
                    // The reference count has saturated.
                    node->detachReferents();
                    to_delete.push_back(node);
                } else {
                    node->detachReferents();
                    decRefCount(node);
                }
            }
        }
    }
    // At this point, the only unmarked objects are GCNodes with saturated
    // reference counts.  Delete them.
    for (GCNode* node : to_delete) {
	delete node;
    }
}

void GCNode::Marker::operator()(const GCNode* node)
{
    if (node->isMarked()) {
	return;
    }
    // Update mark  Beware ~ promotes to unsigned int.
    node->m_rcmms &= static_cast<unsigned char>(~s_mark_mask);
    node->m_rcmms |= s_mark;
    node->visitReferents(this);
}

void GCNode::CountingMarker::operator()(const GCNode* node)
{
    if (node->isMarked()) {
	return;
    }
    Marker::operator()(node);
    m_marks_applied += 1;
}

void rho::initializeMemorySubsystem()
{
    static bool initialized = false;
    if (!initialized) {
	MemoryBank::initialize();
	GCNode::initialize();
	ProtectStack::initialize();
	RAllocStack::initialize();

	initialized = true;
    }
}

// Test if the argument is a GCNode pointer.
GCNode* GCNode::asGCNode(void* candidate_pointer)
{
    if (candidate_pointer < heap_start || candidate_pointer > heap_end) {
	return nullptr;
    }

    void* base_pointer = BlockPool::lookup(candidate_pointer);
#ifdef ALLOCATION_CHECK
    if (base_pointer != lookup_in_allocation_map(candidate_pointer)) {
        error("allocation map mismatch");
    }
#endif
    return static_cast<GCNode*>(base_pointer);
}

GCNode::InternalData GCNode::storeInternalData() const {
    return m_rcmms;
}

void GCNode::restoreInternalData(InternalData data) {
    m_rcmms = data;
}

#ifdef ALLOCATION_CHECK
void add_to_allocation_map(void* allocation, size_t size)
{
    if (iterating) {
        error("allocating node during GC pass");
    }
    if (allocation < heap_start) {
        heap_start = allocation;
    }
    void* allocation_end = static_cast<char*>(allocation) + size;
    if (allocation_end > heap_end) {
        heap_end = allocation_end;
    }
    allocations[allocation] = allocation_end;
}

void remove_from_allocation_map(void* allocation)
{
    if (iterating) {
        error("freeing node during GC pass");
    }
    allocations.erase(allocation);
}

void* lookup_in_allocation_map(void* tentative_pointer)
{
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
