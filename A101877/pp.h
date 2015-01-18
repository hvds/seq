#ifndef PP_H
#define PP_H

#include <string.h>
#include "mygmp.h"
#include "walker.h"

typedef struct pp_s_value {
	mpz_t value;
	int parent;
	int inv;
} pp_value;

typedef struct pp_s_pp {
	int p;
	int pp;
	int depend;
	int valsize;
	int valmax;
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
