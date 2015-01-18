#include <stdlib.h>
#include "mbh.h"

bhsize_t* mbharena;
bhsize_t mbhmaxsize;
bhsize_t mbhsize;
#define MINARENA 100

#define OVERHEAD (sizeof(bh_heap) / sizeof(bhsize_t))

void mbh_init(void) {
	mbhmaxsize = MINARENA;
	mbharena = (bhsize_t*)malloc(mbhmaxsize * sizeof(bhsize_t));
	mbhsize = 0;
}

void mbh_final(void) {
	free(mbharena);
}

void mbh_grow(bhsize_t minsize) {
	if (minsize <= mbhmaxsize)
		return;
	mbhmaxsize = mbhmaxsize * 3 / 2;
	if (mbhmaxsize < minsize)
		mbhmaxsize = minsize;
	mbharena = realloc(mbharena, mbhmaxsize * sizeof(bhsize_t));
}

void mbh_assert(bhp h, bhsize_t size) {
	mbhsize = h + size + OVERHEAD;
	mbh_grow(mbhsize);
}

bhp mbh_new(void* context) {
	bhp h = mbhsize;
	mbh_assert(h, 0);
	BHP(h)->size = 0;
	BHP(h)->context = context;
	return h;
}

void mbh_delete(bhp h) {
	mbhsize = h;
}

static inline void mbh_swapnode(bhp h, bhsize_t left, bhsize_t right) {
	void* temp = BHP(h)->heap[left];
	BHP(h)->heap[left] = BHP(h)->heap[right];
	BHP(h)->heap[right] = temp;
}

void mbh_insert(bhp h, void* v) {
	bhsize_t node = BHP(h)->size++;
	bhsize_t parent;
	mbh_assert(h, node + 1);
	BHP(h)->heap[node] = v;
	while (node > 0) {
		parent = (node - 1) >> 1;
		if (mbh_compare(h, parent, node) <= 0)
			return;
		mbh_swapnode(h, parent, node);
		node = parent;
	}
	return;
}

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
		if (child + 1 < BHP(h)->size
				&& mbh_compare(h, child, child + 1) > 0)
			child = child + 1;
		if (mbh_compare(h, node, child) <= 0)
			break;
		mbh_swapnode(h, node, child);
		node = child;
	}
	return value;
}

bhsize_t mbh_size(bhp h) {
	return BHP(h)->size;
}
