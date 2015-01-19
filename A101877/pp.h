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

typedef struct pp_s_pp {
	int p;				/* prime p */
	int pp;				/* prime power p^k */
	int depend;			/* TRUE if dependent on higher pp */
	int valsize;		/* number of values stored */
	int valmax;			/* max number of values storable without realloc */
	int valnumsize;		/* mpx size of values */
	mpx_add_func* adder; /* mpx add function */
	mpx_cmp_func* cmper; /* mpx compare function */
	pp_value* value;	/* container for the values stored */
	mpz_t min_discard;	/* minimum discard is always zero if dependent */
	mpz_t total;		/* sum of the values */
	mpz_t denominator;	/* common denominator for values */
	int invtotal;		/* sum of inverse of values (mod p) */
	int invdenom;		/* inverse of denominator (mod p) */
	whp wh;				/* handle of current walker */
	wrhp wrh;			/* handle of latest walk_result */
	int wrnum;			/* index of latest walk_result */
	int wrcount;		/* count of results for previous walk */
	mpq_t spare;		/* total to discard at this level */
} pp_pp;

#ifdef ALL_C
inline mpx_t ppv_mpx(pp_value* ppv);
inline int pp_valsize_n(int numsize);
inline int pp_valsize(pp_pp* pp);
inline pp_value* ppv_value_n_i(pp_value* ppv, int numsize, int index);
inline pp_value* pp_value_i(pp_pp* pp, int index);
#else /* ALL_C */
extern inline mpx_t ppv_mpx(pp_value* ppv) {
	return (mpx_t)&ppv->limbs[0];
}
extern inline int pp_valsize_n(int numsize) {
	return sizeof(pp_value) + numsize * sizeof(mp_limb_t);
}
extern inline int pp_valsize(pp_pp* pp) {
	return pp_valsize_n(pp->valnumsize);
}
extern inline pp_value* ppv_value_n_i(pp_value* ppv, int numsize, int index) {
	return (pp_value*)((char*)ppv + index * pp_valsize_n(numsize));
}
extern inline pp_value* pp_value_i(pp_pp* pp, int index) {
	return ppv_value_n_i(pp->value, pp->valnumsize, index);
}
#endif /* ALL_C */

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
