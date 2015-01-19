#include <stdlib.h>
#include <stdio.h>
#include "walker.h"
#include "pp.h"
#include "mbh.h"

static inline void* wr_abstract(walker* w, walk_result* wr) {
	return I2P((char*)wr - w->arena);
}

static inline walk_result* wr_concrete(walker* w, void* abstract) {
	return (walk_result*)(w->arena + P2I(abstract));
}

int mbh_compare_wr(void* context, void* left, void* right) {
	walker* w = (walker*)context;
	walk_result* wl = wr_concrete(w, left);
	walk_result* wr = wr_concrete(w, right);
	return w->cmper(NEXT_DISCARD(w, wl), NEXT_DISCARD(w, wr));
}

void setup_walker(void) {
	setup_mbh();
}

void teardown_walker(void) {
	teardown_mbh();
}

static inline int wr_size(walker* w) {
	return sizeof(walk_result) + w->vecsize * sizeof(int)
			+ w->numsize * sizeof(mp_limb_t) * 2;
}

static inline walk_result* wr_address(walker* w, int i) {
	return (walk_result*)&w->arena[i * wr_size(w)];
}

walk_result* w_pick_arena(walker* w) {
	walk_result* wr;
	int i;

	if (!w->arenanext) {
		int newsize = w->arenamax * 3 / 2;
		if (newsize < MINARENA)
			newsize = MINARENA;
		w->arena = realloc(w->arena, newsize * wr_size(w));
		for (i = w->arenamax; i < newsize; ++i) {
			wr = wr_address(w, i);
			wr->next = wr_address(w, i + 1);
		}
		wr->next = (walk_result*)NULL;
		w->arenanext = wr_address(w, w->arenamax);
		w->arenamax = newsize;
	}
	wr = w->arenanext;
	w->arenanext = wr->next;
	return wr;
}

walker* new_walker(pp_pp* pp, mpz_t limit, int invsum) {
	walker* w;
	walk_result* wr;
	int numsize = pp->valnumsize;

	w = malloc(sizeof(walker) + numsize * sizeof(mp_limb_t) * 2);
	w->heap = mbh_new((void*)w, &mbh_compare_wr);
	w->pp = pp;
	w->numsize = numsize;
	w->adder = pp->adder;
	w->cmper = pp->cmper;
	w->invsum = invsum;
	w->vecsize = (pp->valsize + 31) >> 5;
	w->arena = (char*)NULL;
	w->arenanext = (walk_result*)NULL;
	w->arenamax = 0;
	w->have_previous = 0;
	mpx_set_z(LIMIT(w), numsize, limit);

	wr = w_pick_arena(w);
	wr->invsum = 0;
	wr->nextbit = pp->valsize;
	mpx_set_ui(NEXT_DISCARD(w, wr), numsize, 0);
	mpx_set_ui(DISCARD(w, wr), numsize, 0);
	memset(VEC(w, wr), 0, w->vecsize * sizeof(int));
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

static inline void wr_setbit(walker* w, walk_result* wr, int bit) {
	int word = bit >> 5;
	int offset = 1 << (bit & 31);
	VEC(w, wr)[word] |= offset;
}

static inline int wr_testbit(walker* w, walk_result* wr, int bit) {
	int word = bit >> 5;
	int offset = 1 << (bit & 31);
	return (VEC(w, wr)[word] & offset) ? 1 : 0;
}

int qualify(walker* w, walk_result* wr) {
	if (w->cmper(LIMIT(w), NEXT_DISCARD(w, wr)) >= 0)
		return 1;
	return 0;
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
	mpx_set(DISCARD(w, split), w->numsize, NEXT_DISCARD(w, next), w->numsize);
	split->invsum = next->invsum;
	memcpy(VEC(w, split), VEC(w, next), w->vecsize * sizeof(int));

	if (next->nextbit == w->pp->valsize) {
		/* first */
		--next->nextbit;
		w->adder(NEXT_DISCARD(w, next), DISCARD(w, next),
				MPX(VALUE_I(w->pp, next->nextbit)));
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
	} else if (wr_testbit(w, next, next->nextbit)) {
		w_free_arena(w, next);
	} else {
		w->adder(NEXT_DISCARD(w, next), DISCARD(w, next),
				MPX(VALUE_I(w->pp, next->nextbit)));
		if (!qualify(w, next)) {
			w_free_arena(w, next);
		} else {
			mbh_insert(w->heap, wr_abstract(w, next));
		}
	}

	split->invsum = (split->invsum + VALUE_I(w->pp, limitbit)->inv) % w->pp->p;
	wr_setbit(w, split, limitbit);

	for (i = w->pp->valsize - 1; i > limitbit; --i) {
		if (!wr_testbit(w, split, i))
			break;
	}
	if (i <= limitbit) {
		w_free_arena(w, split);
	} else {
		split->nextbit = i;
		w->adder(NEXT_DISCARD(w, split), DISCARD(w, split),
				MPX(VALUE_I(w->pp, i)));
		if (!qualify(w, split)) {
			w_free_arena(w, split);
		} else {
			mbh_insert(w->heap, wr_abstract(w, split));
		}
	}
	return split;
}

walk_result* walker_findnext(walker* w) {
	while (1) {
		walk_result* next = walker_next(w);
		if (!next)
			return next;
		if (next->invsum != w->invsum)
			continue;
		if (w->have_previous && w->cmper(DISCARD(w, next), PREVIOUS(w)) == 0)
			continue;
		w->have_previous = 1;
		mpx_set(PREVIOUS(w), w->numsize, DISCARD(w, next), w->numsize);
		return next;
	}
}

walk_result* wr_clone(walker* w, walk_result* wr) {
	walk_result* wr2 = malloc(wr_size(w));
	memcpy(wr2, wr, wr_size(w));
	return wr2;
}

void wr_clone_free(walk_result* wr) {
	free(wr);
}
