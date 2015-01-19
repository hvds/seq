#include <stdlib.h>
#include <stdio.h>
#include "pp.h"
#include "prime.h"
#include "inverse.h"
#include "walker.h"

extern double timing(void);

#ifdef DEBUG
#include <stdarg.h>
static inline Dprintf(char* format, ...) {
	va_list ap;
	va_start(ap, format);
	gmp_vprintf(format, ap);
	va_end(ap);
}
#else
static inline Dprintf(char* format, ...) {
}
#endif

/* We rely on calculating ab (mod p) by simple multiplication */
#define MAX_P 65536

pp_pp* pppp = (pp_pp*)NULL;
pp_pp** pplist;
int pplistsize;
int pplistmax;
int k0;
mpz_t z_no_previous;

volatile char diag_signal_seen;

/* I can never remember the name of this GMP function */
static inline unsigned int mod_ui(mpz_t z, unsigned int u) {
	return mpz_fdiv_ui(z, u);
}

/* A surprising omission from the exported mpq_* functions */
static inline void mpq_mul_z(mpq_t dest, mpq_t src1, mpz_t src2) {
	mpz_mul(mpq_numref(dest), mpq_numref(src1), src2);
	if (&dest != &src1)
		mpz_set(mpq_denref(dest), mpq_denref(src1));
	mpq_canonicalize(dest);
}

void setup_pp(int k) {
	int i, d, q;

	k0 = k;
	setup_inverse();
	setup_walker();
	pppp = calloc(k + 1, sizeof(pp_pp));
	pplist = (pp_pp**)NULL;
	pplistsize = 0;
	pplistmax = 0;
	ZINIT(&z_no_previous, "z_no_previous");
	mpz_set_si(z_no_previous, -1);
	
	for (i = 1; i <= k; ++i) {
		int prime, power;
		power = greatest_prime_power(i, &prime);
		pp_save_r(i, prime, power);
	}
}

void pp_free(pp_pp* pp) {
	int i;
	free(pp->value);
	QCLEAR(&pp->spare, "pp_%d.spare", pp->pp);
	ZCLEAR(&pp->min_discard, "pp_%d.min_discard", pp->pp);
	ZCLEAR(&pp->denominator, "pp_%d.denominator", pp->pp);
	ZCLEAR(&pp->total, "pp_%d.total", pp->pp);
	if (pp->wr)
		wr_clone_free(pp->wr);
}

void teardown_pp(void) {
	int i;
	pp_pp* pp;

	ZCLEAR(&z_no_previous, "z_no_previous");
	/* non-NULL pp->wr is a clone (needing free()) iff pp is not in pplist */
	for (i = 0; i < pplistsize; ++i)
		pplist[i]->wr = (walk_result*)NULL;
	for (i = 0; i <= k0; ++i)
		if (pppp[i].pp)
			pp_free(&pppp[i]);
	free(pplist);
	free(pppp);
	teardown_walker();
	teardown_inverse();
}

/*
 * Qsort comparison function for pplist
 * Dependent pp sort in descending order; independent pp sort in ascending
 * order; dependent pp sort after only those independent pp on which they
 * depend.
 */
int pp_listcmp(const void* a, const void* b) {
	pp_pp* left = *(pp_pp**)a;
	pp_pp* right = *(pp_pp**)b;

	if (left->depend) {
		if (right->depend) {
			/* dependent PP sort in descending order */
			return right->pp - left->pp;
		} else {
			/* dependent PP sorts after any independent it depends on */
			return (left->pp * right->pp <= k0) ? +1 : -1;
		}
	} else { /* !left->depend */
		if (right->depend) {
			/* independent PP sorts before any PP dependent on it */
			return (left->pp * right->pp <= k0) ? -1 : +1;
		} else {
			/* independent PP sort in ascending order */
			return left->pp - right->pp;
		}
	}
}

void pp_grow(pp_pp* pp, int size) {
	int i;
	if (pp->valmax < size) {
		pp->value = realloc(pp->value, size * VALSIZE(pp));
		pp->valmax = size;
	}
}

void pp_grownum(pp_pp* pp, int newsize) {
	int i;
	pp_value *newpv, *newi, *oldi;
	mpx_support* xsup = mpx_support_n(newsize);

	newpv = malloc(pp->valmax * VALSIZE_N(newsize));
	for (i = 0; i < pp->valsize; ++i) {
		newi = VALUE_N_I(newpv, newsize, i);
		oldi = VALUE_I(pp, i);
		newi->self = oldi->self;
		newi->parent = oldi->parent;
		newi->inv = oldi->inv;
		mpx_set(MPX(newi), newsize, MPX(oldi), pp->valnumsize);
	}
	free(pp->value);
	pp->valnumsize = newsize;
	pp->value = newpv;
	pp->adder = xsup->adder;
	pp->cmper = xsup->cmper;
}

/*
 * Include a PP structure in the pplist.
 * Inputs:
 *   pp_pp* pp: a pointer to the PP structure to include
 * Returns:
 *   Nothing.
 * Notes:
 *   The list is kept in descending order of pp, and sorted later using
 *   pp_listcmp().
 */
void pp_pushlist(pp_pp* pp) {
	if (pplistsize + 1 >= pplistmax) {
		int old = pplistmax;
		pplistmax = pplistmax * 3 / 2;
		if (pplistmax < MINPPSET)
			pplistmax = MINPPSET;
		pplist = realloc(pplist, pplistmax * sizeof(pplist[0]));
	}
	if (pplistsize)
		memmove(&pplist[1], &pplist[0], pplistsize * sizeof(pplist[0]));
	++pplistsize;
	pplist[0] = pp;
}

/*
 * Remove a PP structure from the pplist.
 * Inputs:
 *   int i: the index in the list to be removed
 * Returns:
 *   Nothing.
 * Notes:
 *   The order of the remaining list items is preserved.
 */
void pp_listsplice(int i) {
	pp_pp* pp = pplist[i];
	if (i < --pplistsize)
		memmove(&pplist[i], &pplist[i + 1],
				(pplistsize - i) * sizeof(pplist[0]));
	/* Cannot free pp, the values may still be referenced */
}

/*
 * Sort the values in a PP structure, given that only the last-added value
 * may be out of order.
 * Input:
 *   pp_pp* pp: pointer to the PP structure to sort.
 * Returns:
 *   Nothing.
 * Caveats:
 *   No guarantees are made about a stable sort: pp_value structures with
 * equal values may be sorted in any order.
 */
void pp_insert_last(pp_pp* pp) {
	int size = pp->valsize - 1;
	pp_value* v = malloc(VALSIZE(pp));
	mpx_support* xsup = mpx_support_n(pp->valnumsize);
	int low, high, med, sign;

	memcpy(v, VALUE_I(pp, size), VALSIZE(pp));
	/* invariant: VALUE_I(pp, low) <= v <= VALUE_I(pp, high) */
	low = -1;
	high = size + 1;
	while (low + 1 < high) {
		med = (low + high) >> 1;
		sign = xsup->cmper(MPX(v), MPX(VALUE_I(pp, med)));
		if (sign < 0)
			low = med;
		else
			high = med;
	}
	if (high + 1 < size) {
		memmove(VALUE_I(pp, high + 1), VALUE_I(pp, high),
				(size - high) * VALSIZE(pp));
		memcpy(VALUE_I(pp, high), v, VALSIZE(pp));
	}
	free(v);
}

/*
 * Store a value in a PP structure.
 * Input:
 *   pp_pp* pp: pointer to the PP structure to store in
 *   mpq_t value: the value to store
 *   int parent: the parent to refer to
 *   int self: TRUE if value represents 1/parent, FALSE if it represents the
 *     resolution of pp[parent]
 * Returns:
 *   Nothing.
 * Effect:
 *   The specified value is stored in the specified PP structure, with a
 * reference to the specified parent. Internal variables are updated
 * accordingly.
 *
 * For a given value q, the caller is expected to have determined the
 * greatest prime power p^k dividing the denominator of q. In this case,
 * the pp structure should be that for p^k.
 *
 * The parent is used solely for determining a final solution vector.
 * If the value is a source value (i.e. 1/r representing the integer
 * r in the solution vector), the parent should be set to r, and self should
 * be TRUE. If not, the parent is assumed to be a fully resolved PP structure
 * storing a walk result that specifies the actual vector contribution
 * (possibly recursively), and self should be FALSE.
 */
void pp_save_any(pp_pp* pp, mpq_t value, int parent, int self) {
	int i;
	mpz_t zvalue;
	pp_value* v;
	ZINIT(&zvalue, "pp_save_any zvalue");

	/* new value = mpq_numref(value) * (denominator / mpq_denref(value))
	 * new_pp() setting of denominator guarantees exactness
	 */
	mpz_divexact(zvalue, pp->denominator, mpq_denref(value));
	mpz_mul(zvalue, zvalue, mpq_numref(value));

	/* adding to total may change required fixed size, at worst once per pp */
	mpz_add(pp->total, pp->total, zvalue);
	if (pp->total->_mp_size > pp->valnumsize) {
		Dprintf("pp_%d: growing to %d with total %Zd (denom %Zd)\n",
				pp->pp, pp->total->_mp_size, pp->total, pp->denominator);
		pp_grownum(pp, pp->total->_mp_size);
	}

	++pp->valsize;
	pp_grow(pp, pp->valsize);
	v = VALUE_I(pp, pp->valsize - 1);
	v->parent = parent;
	v->self = self;
	v->inv = (pp->invdenom * mod_ui(zvalue, pp->p)) % pp->p; /* ref: MAX_P */
	mpx_set_z(MPX(v), pp->valnumsize, zvalue);
	pp->invtotal = (pp->invtotal + v->inv) % pp->p;
	pp_insert_last(pp);
	ZCLEAR(&zvalue, "pp_save_any zvalue");
}

void new_pp(pp_pp* pp, int prime, int power) {
	int i, set;
	mpz_t reduced_denom;
	mpx_support* xsup;

	if (prime > MAX_P) {
		fprintf(stderr, "prime %d exceeds max %d\n", prime, MAX_P);
		exit(1);
	}
	pp->p = prime;
	pp->pp = power;
	if (pp->p == pp->pp)
		inverse_table(pp->p);
	pp_pushlist(pp);
	ZINIT(&pp->total, "pp_%d.total", pp->pp);
	ZINIT(&pp->denominator, "pp_%d.denominator", pp->pp);
	ZINIT(&pp->min_discard, "pp_%d.min_discard", pp->pp);
	QINIT(&pp->spare, "pp_%d.spare", pp->pp);

	/* required denominator is the lcm of (power, S) where S contains those
	 * pp such that pp * power <= k.
	 */
	ZINIT(&reduced_denom, "new_pp reduced_denom");
	set = 0;
	/* start from pplist[1], because we've already been inserted at [0] */
	for (i = 1; i < pplistsize; ++i) {
		if (pplist[i]->pp * power <= k0) {
			mpz_mul_ui(pp->denominator, pplist[i]->denominator,
					power / mpz_gcd_ui(NULL, pplist[i]->denominator, power));
			set = 1;
			break;
		}
	}
	if (!set) {
		/* power == 1 */
		if (power != 1) {
			fprintf(stderr, "new_pp: did not find denominator for %d\n", power);
			exit(1);
		}
		mpz_set_ui(pp->denominator, 1);
	}
	mpz_divexact_ui(reduced_denom, pp->denominator, power);
	pp->invdenom = invfast(mod_ui(reduced_denom, pp->p), pp->p);
	ZCLEAR(&reduced_denom, "new_pp reduced_denom");

	xsup = mpx_support_z(pp->denominator);
	pp->valnumsize = xsup->size;
	pp->adder = xsup->adder;
	pp->cmper = xsup->cmper;
	pp_grow(pp, MINPPSET);
}

void pp_save_r(int n, int prime, int power) {
	pp_pp* pp = &pppp[power];
	pp_value* v;
	int i, raise;
	mpq_t q;
	if (!pp->p)
		new_pp(pp, prime, power);
	QINIT(&q, "pp_save_r temp");
	mpq_set_ui(q, 1, n);
	pp_save_any(pp, q, n, 1);
	QCLEAR(&q, "pp_save_r temp");
}

void pp_save_w(walker* w, walk_result* wr, pp_pp* from) {
	mpq_t actual;
	mpz_t discard;
	int pp_effective = from->pp;
	int pp_applied = 1;
	int pp_p = from->p;
	int q, to, shared;
	pp_pp* pp;
	pp_value* v;

	ZINIT(&discard, "pp_save_w discard");
	QINIT(&actual, "pp_save_w temp");
	mpz_set_x(discard, DISCARD(w, wr), from->valnumsize);
	Dprintf("save wr: %Zd / %Zd [%d] from %d (%d)\n",
			discard, from->denominator, VEC(w, wr)[0], from->pp, from->p);
	from->wr = wr;
	mpz_sub(mpq_numref(actual), from->total, discard);
	mpz_set(mpq_denref(actual), from->denominator);
	mpq_canonicalize(actual);

	to = z_greatest_prime_power(mpq_denref(actual), (int*)NULL);
	pp_save_any(&pppp[to], actual, from->pp, 0);
	QCLEAR(&actual, "pp_save_w temp");
	ZCLEAR(&discard, "pp_save_w discard");
}

/*
 * Resolve simple cases from the pp list
 * For each PP structure, in descending order of pp value:
 * - if a walker finds no way to keep any of the values, discard the PP
 * - if a walker find precisely one way to keep some values:
 *   - save a clone of the walk_result representing that way in the PP structure
 *   - add the relevant sum to the appropriate PP structure
 *   - remove this PP from the pp list
 *   - note that a) this is safe to do only because we know this PP's
 *     total can no longer change, so valnumsize is fixed; b) we need to
 *     take some care to free the walk_result at the appropriate time
 * - for remaining PP structures, mark whether they are dependent on any
 *   higher pp in the pp list
 */
void pp_resolve_simple(void) {
	int i, dest;
	pp_pp *pp;
	walker* w;
	walk_result *wr1, *wr2;
	mpz_t ztotal;

	ZINIT(&ztotal, "pp_resolve_simple ztotal");
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (i == 0 || pp->pp * pplist[i - 1]->pp > k0) {
			/* we are only interested in subsets that keep at least one
			 * element, so we set the discard limit to (total - 1)
			 */
			mpz_sub_ui(ztotal, pp->total, 1);
			w = new_walker(pp, ztotal, pp->invtotal);
			wr1 = walker_findnext(w);
			if (!wr1) {
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			wr1 = wr_clone(w, wr1);
			wr2 = walker_findnext(w);
			if (!wr2) {
				pp_save_w(w, wr1, pp);
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			mpz_set_x(pp->min_discard, DISCARD(w, wr1), pp->valnumsize);
			wr_clone_free(wr1);
			delete_walker(w);
		}
		for (dest = pplistsize - 1; dest > i; --dest) {
			if (pp->pp * pplist[dest]->pp <= k0)
				pplist[dest]->depend = 1;
			else
				break;
		}
	}
	ZCLEAR(&ztotal, "pp_resolve_simple ztotal");
}

void pp_init_spare(int target) {
	int i;
	pp_pp *pp, *top;
	mpq_t spare;

	QINIT(&spare, "pp_study spare");
	top = pplist[0];
	mpq_set_si(top->spare, -target, 1);
	for (i = 0; i < pplistsize; ++i) {
		int g;

		pp = pplist[i];
		mpz_sub(mpq_numref(spare), pp->total, pp->min_discard);
		mpz_set(mpq_denref(spare), pp->denominator);
		mpq_canonicalize(spare);
		mpq_add(top->spare, top->spare, spare);
	}
	QCLEAR(&spare, "pp_study spare");
	Dprintf("study: spare = %Qd\n", top->spare);
}

void pp_study(int target) {
	int i;
	pp_pp *pp, *top;
	mpq_t spare;

	/* delete independent PP that provide no lines
	 * resolve independent PP that provide only one line
	 * mark dependent PP as dependent
	 */
	pp_resolve_simple();

	/* now sort pplist so as to satisfy dependencies as early as possible,
	 * and handle dependent PP as soon as dependencies are satisfied
	 */
	qsort(pplist, pplistsize, sizeof(pplist[0]), &pp_listcmp);

	/* set initial spare ready for pp_find() recursion
	 */
	pp_init_spare(target);
}

static inline void pp_setbit(int* vec, int bit) {
	int byte = bit >> 5;
	int offset = 1 << (bit & 31);
	vec[byte] |= offset;
}

static inline int pp_testbit(int* vec, int bit) {
	int byte = bit >> 5;
	int offset = 1 << (bit & 31);
	return (vec[byte] & offset) ? 1 : 0;
}

void pp_setvec(int* vec, pp_pp* pp, int* invec) {
	int i;
	pp_value* v;

	for (i = 0; i < pp->valsize; ++i) {
		if (invec && pp_testbit(invec, i)) 
			continue;
		v = VALUE_I(pp, i);
		if (v->self) {
			pp_setbit(vec, v->parent);
		} else {
			/* this entry comprises multiple values from the parent */
			pp_setvec(vec, &pppp[v->parent],
					VEC_N(pppp[v->parent].wr, pppp[v->parent].valnumsize));
		}
	}
}

int pp_solution(int index) {
	pp_pp *pp, *cur;
	int* v;
	int i, j, r;
	mpz_t discard;

	cur = pplist[index];
	if (mpz_cmp_ui(mpq_numref(cur->spare), 1) > 0)
		return 0;
	r = mpz_get_ui(mpq_numref(cur->spare));
	if (r == 1 && mpz_cmp_ui(mpq_denref(cur->spare), k0) > 0)
		return 0;
	r = (r == 0) ? 0 : mpz_get_ui(mpq_denref(cur->spare));

	/* probable solution: spare == 0 or spare == 1/r, 1 <= r <= k0 */
	ZINIT(&discard, "pp_solution discard");
	gmp_printf("probable solution, spare = %Qd\n", cur->spare);
	v = calloc((k0 + 32) >> 5, sizeof(int));
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (pp->wr) {
			pp_setvec(v, pp, VEC_N(pp->wr, pp->valnumsize));
		} else if (i >= index && mpz_cmp_ui(pp->min_discard, 0) != 0) {
			/* our solution assumes the unrecorded best line for this pp */
			walker* w = new_walker(pp, pp->total, pp->invtotal);
			walk_result* wr = walker_findnext(w);
			if (!wr) {
				gmp_printf("Solution cannot recreate min_discard for pp_%d\n",
						pp->pp);
				exit(1);
			}
			mpz_set_x(discard, DISCARD(w, wr), pp->valnumsize);
			if (mpz_cmp(pp->min_discard, discard) != 0) {
				gmp_printf("Solution cannot recreate min_discard %Zd for pp_%d "
						"(got %Zd)\n", pp->min_discard, pp->pp, discard);
				exit(1);
			}
			pp_setvec(v, pp, VEC(w, wr));
			delete_walker(w);
		} else {
			pp_setvec(v, pp, (int*)NULL);
		}
	}
	ZCLEAR(&discard, "pp_solution discard");

	if (pp_testbit(v, r)) {
		/* exact solution */
		v[r >> 5] &= ~(1 << (r & 31));
	} else {
		/* inexact, maybe fixable */
		printf("inexact(%d): ", r);
	}
	for (i = 0; i <= k0; ++i)
		if (pp_testbit(v, i))
			printf("%d ", i);
	printf("\nSuccess after %.2fs\n", timing());
	free(v);
	return 1;
}

void pp_diagnose(int level) {
	int j;
	int start = -1, end;
	pp_pp* pp;

	for (j = 0; j <= level; ++j) {
		pp = pplist[j];
		if (pp->wrnum != 1 || pp->wrcount != 0) {
			start = j;
			break;
		}
	}
	end = start + 10;
	if (end > level)
		end = level;
	printf("t=%.2f level %d..%d at %d of %d: ",
			timing(), start, end, level, pplistsize);
	for (j = start; j <= end; ++j) {
		pp = pplist[j];
		printf("%d:%d/%d ", pp->pp, pp->wrnum, pp->wrcount);
	}
	printf("\n");
}

int pp_find(int target) {
	pp_pp *pp, *nextpp;
	int i, g, invsum;
	int success = 0;
	int have_discard;
	mpz_t z, zr;
	mpq_t q, limit;

	pp_study(target);

	i = 0;
	pplist[i]->w = (walker*)NULL;
	if (mpq_sgn(pplist[i]->spare) < 0) {
		printf("n=%d, k=%d: optimizer finds no solution is possible. [%.2fs]\n",
				target, k0, timing());
		return 0;
	}
	ZINIT(&z, "pp_find z");
	ZINIT(&zr, "pp_find zr");
	QINIT(&q, "pp_find q");
	QINIT(&limit, "pp_find limit");
	while (i >= 0) {
		pp = pplist[i];
		if (!pp->w) {
			if (diag_signal_seen) {
				diag_signal_seen = 0;
				pp_diagnose(i);
			}

			mpq_set(limit, pp->spare);
			Dprintf("at pp_%d spare = %Qd\n", pp->pp, limit);

			have_discard = (!pp->depend && mpz_sgn(pp->min_discard) != 0);
			if (have_discard) {
				/* (limit_r) += min_discard / denominator */
				mpz_set(mpq_numref(q), pp->min_discard);
				mpz_set(mpq_denref(q), pp->denominator);
				mpq_canonicalize(q);
				mpq_add(limit, limit, q);
			}
			/* actual limit = floor(limit_num * denominator / limit_den) */
			mpq_mul_z(limit, limit, pp->denominator);
			mpz_fdiv_q(mpq_numref(limit), mpq_numref(limit), mpq_denref(limit));
			/* limit is now denormal, but we won't be using it again */

			if (pp->depend) {
				/* min_discard is 0, so nothing extra factored into spare:
				 * must discard spare (mod p) = spare_num * inv(spare_den) */
				mpz_fdiv_qr_ui(z, zr, mpq_denref(pp->spare), pp->pp);
				if (mpz_sgn(zr) == 0) {
					mpz_mul_ui(z, mpq_numref(pp->spare),
							invfast(mod_ui(z, pp->p), pp->p));
					invsum = mod_ui(z, pp->p);
				} else {
					invsum = 0;
				}
			} else {
				invsum = pp->invtotal;
			}
			Dprintf("effective limit is %Zd, invsum = %d\n",
					mpq_numref(limit), invsum);

			pp->w = new_walker(pp, mpq_numref(limit), invsum);
			pp->wrnum = 0;
		}

		pp->wr = walker_findnext(pp->w);
		if (!pp->wr) {
			delete_walker(pp->w);
			pp->w = (walker*)NULL;
			pp->wrcount = pp->wrnum;
			--i;
			Dprintf("no lines for %d, recurse back to %d\n",
					pp->pp, (i >= 0) ? pplist[i]->pp : -1);
			continue;
		}
		++pp->wrnum;

		++i;
		if (i >= pplistsize) {
			fprintf(stderr, "overflow\n");
			exit(1);
		}
		nextpp = pplist[i];

		/* new spare = old spare - (actual discard - min discard) / denominator
		 */
		mpz_set_x(mpq_numref(q), DISCARD(pp->w, pp->wr), pp->valnumsize);
		Dprintf("pp_%d walker found discard %Zd/%Zd\n",
				pp->pp, mpq_numref(q), pp->denominator);
		mpz_sub(mpq_numref(q), mpq_numref(q), pp->min_discard);
		mpz_set(mpq_denref(q), pp->denominator);
		mpq_canonicalize(q);
		mpq_sub(nextpp->spare, pp->spare, q);

		Dprintf("new spare %Qd\n", nextpp->spare);

		if (mpq_sgn(nextpp->spare) < 0) {
			/* we've overrun */
			--i;
			continue;
		}
		nextpp->w = (walker*)NULL;
		if (pp_solution(i)) {
			success = 1;
			break;
		}
	}
	QCLEAR(&limit, "pp_find limit");
	QCLEAR(&q, "pp_find q");
	ZCLEAR(&zr, "pp_find zr");
	ZCLEAR(&z, "pp_find z");
	for (i = 0; i < pplistsize; ++i) {
		if (pplist[i]->w)
			delete_walker(pplist[i]->w);
	}
	if (!success)
		printf("n=%d k=%d: no solution found. [%.2fs]\n",
                target, k0, timing());
	return success;
}
