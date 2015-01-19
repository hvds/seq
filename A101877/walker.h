#ifndef WALKER_H
#define WALKER_H

typedef int wrhp;
typedef int whp;

#include "pp.h"
#include "mbh.h"
#include "mygmp.h"

typedef struct s_walk_result {
	wrhp next;
	int invsum;
	int nextbit;
	mp_limb_t tail[0];
	/* tail consists of:
		mp_limb_t next_discard_limbs[w->numsize];
		mp_limb_t discard_limbs[w->numsize];
		int vec[w->vecsize];
	*/
} walk_result;

typedef struct s_walker {
	bhp heap;
	struct pp_s_pp* pp;
	int numsize;
	mpx_add_func* adder;
	mpx_cmp_func* cmper;
	int invsum;
	int vecsize;
	wrhp arenanext;
	int have_previous;
	mp_limb_t tail[0];
	/* tail consists of:
		mp_limb_t limit_limbs[w->numsize];
		mp_limb_t previous[w->numsize];
	*/
} walker;

extern char* walk_arena;

#ifdef ALL_C
inline walker* WP(whp wh);
inline walk_result* WRP(whp wh, wrhp wrh);
inline mpx_t wr_discard(whp wh, wrhp wrh);
inline int* wr_vec(whp wh, wrhp wrh);
#else /* ALL_C */
extern inline walker* WP(whp wh) {
	return (walker*)&walk_arena[wh];
}
extern inline walk_result* WRP(whp wh, wrhp wrh) {
	return (walk_result*)&walk_arena[wrh];
}

extern inline mpx_t wr_discard(whp wh, wrhp wrh) {
	return (mpx_t)&WRP(wh, wrh)->tail[WP(wh)->numsize];
}
extern inline int* wr_vec(whp wh, wrhp wrh) {
	return (int*)&WRP(wh, wrh)->tail[WP(wh)->numsize * 2];
}
#endif /* ALL_C */

#define MINARENA 16

extern void setup_walker(void);
extern void teardown_walker(void);
extern whp new_walker(struct pp_s_pp* pp, mpz_t limit, int invsum);
extern wrhp walker_next(whp wh);
extern wrhp walker_findnext(whp wh);
extern void delete_walker(whp wh);
extern wrhp wr_clone(whp wh, wrhp wrh);
extern void wr_clone_free(whp wh, wrhp wrh);

#endif /* WALKER_H */
