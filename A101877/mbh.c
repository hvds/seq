/*
 * Support for a stack of binary (min-)heaps.
 *
 * This is tailored specifically for use by a recursive algorithm: each new
 * heap sits on top of previous heaps, so only the heap most recently
 * returned by mbh_new() (and not yet freed by mbh_delete()) should ever be
 * modified.
 *
 * Sorting is done by a user-supplied comparison routine:
 *   extern int mbh_compare(void* context, void* left, void* right);
 * which should return a negative integer to represent 'left < right',
 * zero to represent 'left == right' and a positive integer to represent
 * 'left > right'.
 */

#include <stdlib.h>
#undef SAFE_BUT_SLOW /* we need the private definitions here if nowhere else */
#include "mbh.h"

bhsize_t* mbharena;
bhsize_t mbhmaxsize;
bhsize_t mbhsize;
#define MBH_MINARENA 100
#define OVERHEAD (sizeof(bh_heap) / sizeof(bhsize_t))

/*
 * Initialize for use. Must be called before any other functions are called.
 */
void setup_mbh(void) {
	mbhmaxsize = MBH_MINARENA;
	mbharena = (bhsize_t*)malloc(mbhmaxsize * sizeof(bhsize_t));
	mbhsize = 0;
}

/*
 * Clean up after use. After calling this, state should be as it was before
 * setup_mbh() was called.
 */
void teardown_mbh(void) {
	free(mbharena);
}

/*
 * Resize the arena.
 * Input:
 *   bhsize_t minsize: grow the arena to be at least minsize sizeof(bhsize_t).
 * Returns:
 *   Nothing.
 */
void mbh_grow(bhsize_t minsize) {
	if (minsize <= mbhmaxsize)
		return;
	mbhmaxsize = mbhmaxsize * 3 / 2;
	if (mbhmaxsize < minsize)
		mbhmaxsize = minsize;
	mbharena = realloc(mbharena, mbhmaxsize * sizeof(bhsize_t));
}

/*
 * Resize a heap.
 * Input:
 *   bhp h: the heap pointer to resize
 *   bhsize_t size: the minimum size the heap should be
 * Returns:
 *   Nothing.
 */
void mbh_assert(bhp h, bhsize_t size) {
	mbhsize = h + size + OVERHEAD;
	mbh_grow(mbhsize);
}

/*
 * Initialize a new heap.
 * Input:
 *   void* context: an opaque context pointer
 * Returns:
 *   A bhp referring to the new heap.
 * Notes:
 *   The context pointer is accessible as BHP(h)->context.
 *   Any bhp acquired before this should not be modified until after this
 *   one has been deleted with mbh_delete().
 */
bhp mbh_new(void* context, mbh_compare_func* comparator) {
	bhp h = mbhsize;
	mbh_assert(h, 0);
	BHP(h)->size = 0;
	BHP(h)->context = context;
	BHP(h)->comparator = comparator;
	return h;
}

/*
 * Free a heap.
 * Input:
 *   bhp h: the pointer to the heap to be freed.
 * Returns:
 *   Nothing.
 */
void mbh_delete(bhp h) {
	mbhsize = h;
}

/*
 * Swap a pair of nodes in a heap.
 * Input:
 *   bhp h: the pointer to the heap.
 *   bhsize_t left: the index of the first node to be swapped.
 *   bhsize_t right: the index of the second node to be swapped.
 * Returns:
 *   Nothing.
 */
static inline void mbh_swapnode(bhp h, bhsize_t left, bhsize_t right) {
	void* temp = BHP(h)->heap[left];
	BHP(h)->heap[left] = BHP(h)->heap[right];
	BHP(h)->heap[right] = temp;
}

/*
 * Insert a new node into a heap.
 * Input:
 *   bhp h: the pointer to the heap.
 *   void* v: an opaque value object to insert.
 * Returns:
 *   Nothing.
 * Notes:
 *   The value becomes referenced by the heap until either it is shifted
 * out (i.e. returned by a call to mbh_shift(h)), or the heap is deleted
 * (with mbh_delete(h)).
 */
void mbh_insert(bhp h, void* v) {
	bhsize_t node = BHP(h)->size++;
	bhsize_t parent;
	mbh_assert(h, node + 1);
	BHP(h)->heap[node] = v;
	while (node > 0) {
		parent = (node - 1) >> 1;
		if (BHP(h)->comparator(
			BHP(h)->context, BHP(h)->heap[parent], BHP(h)->heap[node]
		) <= 0) {
			return;
		}
		mbh_swapnode(h, parent, node);
		node = parent;
	}
	return;
}

/*
 * Shift the least node out of the heap.
 * Input:
 *   bhp h: the pointer to the heap.
 * Returns:
 *   void* v, the opaque value object shift out of the heap.
 */
void* mbh_shift(bhp h) {
	void* value = BHP(h)->heap[0];
	bhsize_t node = --BHP(h)->size;
	bhsize_t child;
	mbh_assert(h, node);
	BHP(h)->heap[0] = BHP(h)->heap[node];
	node = 0;
	while (1) {
		child = (node << 1) + 1;
		if (child >= BHP(h)->size)
			break;
		if (child + 1 < BHP(h)->size && BHP(h)->comparator(
			BHP(h)->context, BHP(h)->heap[child], BHP(h)->heap[child + 1]
		) > 0) {
			child = child + 1;
		}
		if (BHP(h)->comparator(
			BHP(h)->context, BHP(h)->heap[node], BHP(h)->heap[child]
		) <= 0) {
			break;
		}
		mbh_swapnode(h, node, child);
		node = child;
	}
	return value;
}

/*
 * The current size of the heap.
 * Input:
 *   bhp h: the pointer to the heap.
 * Returns:
 *   bhsize_t size, the number of items currently stored in the heap.
 */
bhsize_t mbh_size(bhp h) {
	return BHP(h)->size;
}
