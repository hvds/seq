#include "walker.h"
#include "pp.h"
#include "mbh.h"
#include <stdlib.h>
#include <stdio.h>

static inline void* wr_abstract(walker* w, walk_result* wr) {
	return (void*)((char*)wr - w->arena);
}

static inline walk_result* wr_concrete(walker* w, void* abstract) {
	return (walk_result*)(w->arena + P2I(abstract));
}

void setup_walker(void) {
	setup_mbh();
}

void teardown_walker(void) {
	teardown_mbh();
}

static inline int wr_size(walker* w) {
	return sizeof(walk_result) + w->vecsize * sizeof(int);
}

static inline walk_result* wr_address(walker* w, int i) {
	return (walk_result*)&w->arena[i * wr_size(w)];
}

walk_result* w_pick_arena(walker* w) {
	walk_result* wr;
	if (!w->arenanext) {
		int newsize = w->arenamax * 3 / 2;
		if (newsize < MINARENA)
			newsize = MINARENA;
		w->arena = realloc(w->arena, newsize * wr_size(w));
		wr = wr_address(w, w->arenamax);
		w->arenanext = wr;
		while (w->arenamax < newsize - 1) {
			wr->next = wr_address(w, ++w->arenamax);
			wr = wr->next;
		}
		++w->arenamax;
		wr->next = (walk_result*)NULL;
	}
	wr = w->arenanext;
	w->arenanext = wr->next;
	return wr;
}

walker* new_walker(pp_pp* pp, int limit, int invsum) {
	walker* w = malloc(sizeof(walker));
	walk_result* wr;

	w->heap = mbh_new((void*)w);
	w->pp = pp;
	w->limit = limit;
	w->invsum = invsum;
	w->vecsize = (pp->valsize + 31) >> 5;
	w->arena = (char*)NULL;
	w->arenanext = (walk_result*)NULL;
	w->arenamax = 0;

	wr = w_pick_arena(w);
	wr->next_discard = 0;
	wr->discard = 0;
	wr->invsum = 0;
	wr->nextbit = pp->valsize;
	memset(&wr->vec[0], 0, w->vecsize * sizeof(int));
	mbh_insert(w->heap, wr_abstract(w, wr));

	return w;
}

void delete_walker(walker* w) {
	mbh_delete(w->heap);
	free(w->arena);
	free(w);
}

void w_free_arena(walker* w, walk_result* wr) {
	wr->next = w->arenanext;
	w->arenanext = wr;
}

static inline void wr_setbit(walk_result* wr, int bit) {
	int byte = bit >> 5;
	int offset = 1 << (bit & 31);
	wr->vec[byte] |= offset;
}

static inline int wr_testbit(walk_result* wr, int bit) {
	int byte = bit >> 5;
	int offset = 1 << (bit & 31);
	return (wr->vec[byte] & offset) ? 1 : 0;
}

int qualify(walker* w, walk_result* wr) {
	wr->next_discard = wr->discard + w->pp->value[wr->nextbit].value;
	if (w->limit == 0 || w->limit >= wr->next_discard)
		return 1;
	return 0;
}

walk_result* walker_find(walker* w) {
	while (1) {
		walk_result* next = walker_next(w);
		if (!next || next->invsum == w->invsum)
			return next;
	}
}

walk_result* walker_next(walker* w) {
	walk_result* next;
	walk_result* split;
	int limitbit, i;

	if (mbh_size(w->heap) == 0)
		return (walk_result*)NULL;

	split = w_pick_arena(w);
	/* must pick_arena first, since arena may move */
	next = wr_concrete(w, mbh_shift(w->heap));
	limitbit = next->nextbit;
	split->discard = next->next_discard;
	split->invsum = next->invsum;
	memcpy(&split->vec[0], &next->vec[0], w->vecsize * sizeof(int));

	if (next->nextbit == w->pp->valsize) {
		/* first */
		--next->nextbit;
		next->next_discard = next->discard + w->pp->value[next->nextbit].value;
		if (!qualify(w, next)) {
			w_free_arena(w, next);
		} else {
			mbh_insert(w->heap, wr_abstract(w, next));
		}
		w_free_arena(w, split);
		return next;
	}
		
	--next->nextbit;
	if (next->nextbit < 0) {
		w_free_arena(w, next);
	} else if (wr_testbit(next, next->nextbit)) {
		w_free_arena(w, next);
	} else {
		next->next_discard = next->discard + w->pp->value[next->nextbit].value;
		if (!qualify(w, next)) {
			w_free_arena(w, next);
		} else {
			mbh_insert(w->heap, wr_abstract(w, next));
		}
	}

	split->invsum = (split->invsum + w->pp->value[limitbit].inv) % w->pp->p;
	wr_setbit(split, limitbit);

	for (i = w->pp->valsize - 1; i > limitbit; --i) {
		if (!wr_testbit(split, i))
			break;
	}
	if (i <= limitbit) {
		w_free_arena(w, split);
	} else {
		split->nextbit = i;
		split->next_discard = split->discard + w->pp->value[i].value;
		if (!qualify(w, split)) {
			w_free_arena(w, split);
		} else {
			mbh_insert(w->heap, wr_abstract(w, split));
		}
	}
	return split;
}

walk_result* wr_clone(walker* w, walk_result* wr) {
	walk_result* wr2 = malloc(wr_size(w));
	memcpy(wr2, wr, wr_size(w));
	return wr2;
}

int mbh_compare(bhp h, bhsize_t left, bhsize_t right) {
	walker* w = (walker*)BHP(h)->context;
	walk_result* wl = wr_concrete(w, BHP(h)->heap[left]);
	walk_result* wr = wr_concrete(w, BHP(h)->heap[right]);
	return (wl->next_discard < wr->next_discard) ? -1
			: (wl->next_discard == wr->next_discard) ? 0 : +1;
}
