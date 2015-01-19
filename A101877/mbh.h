#ifndef MBH_H
#define MBH_H

#include <stdint.h>

typedef size_t bhsize_t;    /* must be same size as void* */
typedef bhsize_t bhp;

/*
 * Typedef for function to compare two heap objects
 * Input:
 *   void* context: the context object supplied to mbh_new()
 *   void* left, void* right: the two value objects to compare
 * Returns:
 *   a negative integer to represent 'left < right'
 *   a positive integer to represent 'left > right'
 *   zero to represent 'left == right'
 */
typedef int (mbh_compare_func)(void* context, void* left, void* right);

/*
 * Initialize for use. Must be called before any other functions are called.
 */
extern void setup_mbh(void);

/*
 * Clean up after use. After calling this, state should be as it was before
 * setup_mbh() was called.
 */
extern void teardown_mbh(void);

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
extern bhp mbh_new(void* context, mbh_compare_func* comparator);

/*
 * Free a heap.
 * Input:
 *   bhp h: the pointer to the heap to be freed.
 * Returns:
 *   Nothing.
 */
extern void mbh_delete(bhp h);

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
extern void mbh_insert(bhp h, void* v);

/*
 * Shift the least node out of the heap.
 * Input:
 *   bhp h: the pointer to the heap.
 * Returns:
 *   void* v, the opaque value object shift out of the heap.
 */
extern void* mbh_shift(bhp h);

#ifdef SAFE_BUT_SLOW
extern bhsize_t mbh_size(bhp h);
#else /* ifdef SAFE_BUT_SLOW */

/* Note: mbharena, bh_heap and BHP(h) are not intended for use by callers.
 * They are exposed here purely so that mbh_size() can be inlined.
 */
extern bhsize_t* mbharena;
typedef struct bh_s_heap {
	bhsize_t size;
	void* context;
	mbh_compare_func* comparator;
	void* heap[0];
} bh_heap;
#define BHP(h) ((bh_heap*)&mbharena[h])

/*   
 * The current size of the heap.
 * Input:
 *   bhp h: the pointer to the heap.
 * Returns:
 *   bhsize_t size, the number of items currently stored in the heap.
 */
#ifdef ALL_C
inline bhsize_t mbh_size(bhp h);
#else /* ALL_C */
extern inline bhsize_t mbh_size(bhp h) {
	return BHP(h)->size;
}
#endif /* ALL_C */
#endif /* ifdef SAFE_BUT_SLOW */

#define P2I(x) (int)(intptr_t)(x)
#define I2P(x) (void*)(intptr_t)(x)

#endif /* MBH_H */
