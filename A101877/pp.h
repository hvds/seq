#ifndef PP_H
#define PP_H

#include <string.h>
#include "mygmp.h"
#include "walker.h"

typedef struct pp_s_value {
	int self;		/* boolean, TRUE if this is the entry for 1/parent, FALSE
					 * if it is the resolution of pp[parent] */
	int parent;		/* if not self, parent is the PP structure this resolves */
	int inv;		/* (MPX(vp) / pp->denominator)^-1 (mod pp->p) */
	/* the actual rational represented is MPX(vp) / pp->denominator */
	mp_limb_t limbs[0];
} pp_value;
#define MPX(vp) ((mpx_t)&((vp)->limbs[0]))
#define VALSIZE_N(n) (sizeof(pp_value) + (n) * sizeof(mp_limb_t))
#define VALSIZE(pp) VALSIZE_N((pp)->valnumsize)
#define VALUE_N_I(v, n, i) ((pp_value*)((char*)(v) + (i) * VALSIZE_N(n)))
#define VALUE_I(pp, i) VALUE_N_I((pp)->value, (pp)->valnumsize, i)

typedef struct pp_s_pp {
	int p;
	int pp;
	int depend;
	int valsize;
	int valmax;
	int valnumsize;
	mpx_add_func* adder;
	mpx_cmp_func* cmper;
	pp_value* value;
	mpz_t min_discard;
	mpz_t total;
	mpz_t denominator;
	int invtotal;
	int invdenom;
	struct s_walker* w;
	struct s_walk_result* wr;
	int wrnum;
	int wrcount;
	mpq_t spare;
} pp_pp;

extern pp_pp* pppp;
extern pp_pp** pplist;
extern int pplistsize;
extern int pplistmax;
#define MINPPSET 10

extern volatile char diag_signal_seen;

extern void setup_pp(int k);
extern void teardown_pp(void);
extern void pp_study(int target);
extern int pp_find(int target);
extern void pp_save_r(int n, int prime, int power);

#endif /* PP_H */
