#ifndef WALKER_H
#define WALKER_H

#include "pp.h"
#include "mbh.h"
#include "mygmp.h"

typedef struct s_walk_result {
	struct s_walk_result* next;
	int invsum;
	int nextbit;
	mp_limb_t tail[0];
	/* tail consists of:
		mp_limb_t discard_limbs[w->numsize];
		mp_limb_t next_discard_limbs[w->numsize];
		int vec[w->vecsize];
	*/
} walk_result;
#define DISCARD(w, wr) ((mpx_t)&(wr)->tail[0])
#define NEXT_DISCARD(w, wr) ((mpx_t)&(wr)->tail[(w)->numsize])
#define VEC(w, wr) ((int*)&(wr)->tail[(w)->numsize * 2])
#define VEC_N(wr, size) ((int*)&(wr)->tail[size * 2])

typedef struct s_walker {
	bhp heap;
	struct pp_s_pp* pp;
	int numsize;
	mpx_add_func* adder;
	mpx_cmp_func* cmper;
	int invsum;
	int vecsize;
	char* arena;
	walk_result* arenanext;
	int arenamax;
	int have_previous;
	mp_limb_t tail[0];
	/* tail consists of:
		mp_limb_t limit_limbs[w->numsize];
		mp_limb_t previous[w->numsize];
	*/
} walker;
#define LIMIT(w) ((mpx_t)&(w)->tail[0])
#define PREVIOUS(w) ((mpx_t)&(w)->tail[(w)->numsize])

#define MINARENA 16

extern void setup_walker(void);
extern void teardown_walker(void);
extern walker* new_walker(struct pp_s_pp* pp, mpz_t limit, int invsum);
extern void w_push_heap(walker* w, mpx_t sum, int nextbit);
extern walk_result* walker_next(walker* w);
extern walk_result* walker_findnext(walker* w);
extern void delete_walker(walker* w);
extern walk_result* wr_clone(walker* w, walk_result* wr);
extern void wr_clone_free(walk_result* wr);

#endif /* WALKER_H */
