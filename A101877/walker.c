#include <stdlib.h>
#include <stdio.h>
#include "walker.h"
#include "pp.h"
#include "mbh.h"

char* walk_arena;
int arena_size;
int arena_max;
#define MIN_WALK_ARENA 4096

inline walker* WP(whp wh) {
	return (walker*)&walk_arena[wh];
}
#define DWP(wh) ((walker*)&walk_arena[wh])

inline walk_result* WRP(whp wh, wrhp wrh) {
	return (walk_result*)&walk_arena[wrh];
}
#define DWRP(wh, wrh) ((walk_result*)&walk_arena[wrh])

static inline whp WHP(walker* w) {
	return (whp)((char*)w - walk_arena);
}

static inline wrhp WRHP(walker* w, walk_result* wr) {
	return (wrhp)((char*)wr - walk_arena);
}

static inline mpx_t w_limit(walker* w) {
	return (mpx_t)&w->tail[0];
}

static inline mpx_t w_previous(walker* w) {
	return (mpx_t)&w->tail[w->numsize];
}

static inline mpx_t wr_next_discard(walker* w, walk_result* wr) {
	return (mpx_t)&wr->tail[0];
}

inline mpx_t wr_discard(whp wh, wrhp wrh) {
	return (mpx_t)&WRP(wh, wrh)->tail[WP(wh)->numsize];
}

static inline mpx_t wr_discard_direct(walker* w, walk_result* wr) {
	return (mpx_t)&wr->tail[w->numsize];
}

inline int* wr_vec(whp wh, wrhp wrh) {
	return (int*)&WRP(wh, wrh)->tail[WP(wh)->numsize * 2];
}

static inline int* wr_vec_direct(walker* w, walk_result* wr) {
	return (int*)&wr->tail[w->numsize * 2];
}

int mbh_compare_wr(void* context, void* left, void* right) {
	whp wh = P2I(context);
	walker* w = DWP(wh);
	wrhp wlh = P2I(left);
	wrhp wrh = P2I(right);
	walk_result* wl = DWRP(wh, wlh);
	walk_result* wr = DWRP(wh, wrh);
	return w->cmper(wr_next_discard(w, wl), wr_next_discard(w, wr));
}

void setup_walker(void) {
	setup_mbh();
	arena_size = 4;	/* must not return 0 as a whp */
	arena_max = MIN_WALK_ARENA;
	walk_arena = malloc(arena_max);
}

void teardown_walker(void) {
	free(walk_arena);
	teardown_mbh();
}

static inline int wr_charsize(walker* w) {
	return (
		sizeof(walk_result)
		+ w->vecsize * sizeof(int)
		+ w->numsize * 2 * sizeof(mp_limb_t)
	);
}

#define grow_arena(w, size) \
	if ((size) > arena_max) { \
		int oldarena = ((char*)w - walk_arena); \
		arena_max = arena_max * 2; \
		walk_arena = realloc(walk_arena, arena_max); \
		w = (walker*)(walk_arena + oldarena); \
	}

#define w_pick_arena(w, wr) \
	{ \
		wrhp pick_wrh; \
		if (w->arenanext) { \
			pick_wrh = w->arenanext; \
			w->arenanext = DWRP(WHP(w), pick_wrh)->next; \
		} else { \
			pick_wrh = arena_size; \
			arena_size += wr_charsize(w); \
			grow_arena(w, arena_size); \
		} \
		wr = DWRP(WHP(w), pick_wrh); \
	}

static inline void push_heap(walker* w, walk_result* wr) {
	mbh_insert(w->heap, I2P(WRHP(w, wr)));
}

static inline walk_result* pop_heap(walker* w) {
	void* v = mbh_shift(w->heap);
	return DWRP(WHP(w), P2I(v));
}

static inline int walker_charsize(int numsize) {
	return sizeof(walker) + numsize * 2 * sizeof(mp_limb_t);
}

whp new_walker(pp_pp* pp, mpz_t limit, int invsum) {
	whp wh;
	walker* w;
	walk_result* wr;
	int numsize = pp->valnumsize;

	wh = arena_size;
	w = WP(wh);
	arena_size += walker_charsize(numsize);
	grow_arena(w, arena_size);
	w->heap = mbh_new(I2P(wh), &mbh_compare_wr);
	w->pp = pp;
	w->numsize = numsize;
	w->adder = pp->adder;
	w->cmper = pp->cmper;
	w->invsum = invsum;
	w->vecsize = (pp->valsize + 31) >> 5;
	w->arenanext = (wrhp)0;
	w->have_previous = 0;
	mpx_set_z(w_limit(w), numsize, limit);

	w_pick_arena(w, wr);
	wr->invsum = 0;
	wr->nextbit = pp->valsize;
	mpx_set_ui(wr_next_discard(w, wr), numsize, 0);
	mpx_set_ui(wr_discard_direct(w, wr), numsize, 0);
	memset(wr_vec_direct(w, wr), 0, w->vecsize * sizeof(int));
	push_heap(w, wr);

	return wh;
}

void delete_walker(whp wh) {
	walker* w = WP(wh);

	mbh_delete(w->heap);
	arena_size = wh;
}

static inline void w_free_arena(walker* w, walk_result* wr) {
	wr->next = w->arenanext;
	w->arenanext = WRHP(w, wr);
}

static inline void wr_setbit(walker* w, walk_result* wr, int bit) {
	int word = bit >> 5;
	int offset = 1 << (bit & 31);
	wr_vec_direct(w, wr)[word] |= offset;
}

static inline int wr_testbit(walker* w, walk_result* wr, int bit) {
	int word = bit >> 5;
	int offset = 1 << (bit & 31);
	return (wr_vec_direct(w, wr)[word] & offset) ? 1 : 0;
}

static inline int qualify(walker* w, walk_result* wr) {
	if (w->cmper(w_limit(w), wr_next_discard(w, wr)) >= 0)
		return 1;
	return 0;
}

wrhp walker_findnext(whp wh) {
	walker* w = DWP(wh);
	walk_result* next;
	walk_result* split;
	wrhp nexth;
	int limitbit, i;

	while (1) {
		if (mbh_size(w->heap) == 0)
			return (wrhp)0;

		/* pick_arena first, since arena may move */
		w_pick_arena(w, split);
		next = pop_heap(w);
		limitbit = next->nextbit;
		mpx_set(wr_discard_direct(w, split), w->numsize,
				wr_next_discard(w, next), w->numsize);
		split->invsum = next->invsum;
		memcpy(wr_vec_direct(w, split), wr_vec_direct(w, next),
				w->vecsize * sizeof(int));

		if (next->nextbit == w->pp->valsize) {
			/* first */
			--next->nextbit;
			w->adder(wr_next_discard(w, next), wr_discard_direct(w, next),
					ppv_mpx(pp_value_i(w->pp, next->nextbit)));
			if (!qualify(w, next)) {
				w_free_arena(w, next);
			} else {
				push_heap(w, next);
			}
			w_free_arena(w, split);
			goto found_line;
		}

		--next->nextbit;
		if (next->nextbit < 0) {
			w_free_arena(w, next);
		} else if (wr_testbit(w, next, next->nextbit)) {
			w_free_arena(w, next);
		} else {
			w->adder(wr_next_discard(w, next), wr_discard_direct(w, next),
					ppv_mpx(pp_value_i(w->pp, next->nextbit)));
			if (!qualify(w, next)) {
				w_free_arena(w, next);
			} else {
				push_heap(w, next);
			}
		}

		split->invsum = (split->invsum + pp_value_i(w->pp, limitbit)->inv)
				% w->pp->p;
		wr_setbit(w, split, limitbit);

		for (i = w->pp->valsize - 1; i > limitbit; --i) {
			if (!wr_testbit(w, split, i))
				break;
		}
		if (i <= limitbit) {
			w_free_arena(w, split);
		} else {
			split->nextbit = i;
			w->adder(wr_next_discard(w, split), wr_discard_direct(w, split),
					ppv_mpx(pp_value_i(w->pp, i)));
			if (!qualify(w, split)) {
				w_free_arena(w, split);
			} else {
				push_heap(w, split);
			}
		}
		next = split;
	  found_line:
		nexth = WRHP(w, next);
		if (next->invsum != w->invsum)
			continue;
		if (w->have_previous
				&& w->cmper(wr_discard_direct(w, next), w_previous(w)) == 0)
			continue;
		w->have_previous = 1;
		mpx_set(w_previous(w), w->numsize,
				wr_discard_direct(w, next), w->numsize);
		return nexth;
	}
}

wrhp wr_clone(whp wh, wrhp wrh) {
	walker* w = WP(wh);
	walk_result *wr, *wr2;

	w_pick_arena(w, wr2);
	/* we may attempt to clone something already released to arena, in which
	 * case pick_arena may return the same address. If so, there's nothing
	 * to do */
	wr = WRP(wh, wrh);
	if (wr2 != wr)
		memcpy(wr2, wr, wr_charsize(w));
	return WRHP(w, wr2);
}

void wr_clone_free(whp wh, wrhp wrh) {
	w_free_arena(WP(wh), WRP(wh, wrh));
}
