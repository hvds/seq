#include <stdlib.h>
#include <stdio.h>
#include "pp.h"
#include "usable.h"
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

pp_n* ppn = (pp_n*)NULL;
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

static inline void mpq_mul_z(mpq_t dest, mpq_t src1, mpz_t src2) {
	mpz_mul(mpq_numref(dest), mpq_numref(src1), src2);
	if (&dest != &src1)
		mpz_set(mpq_denref(dest), mpq_denref(src1));
	mpq_canonicalize(dest);
}

void setup_pp(int k) {
	int i, d, q;
	pp_n* ppi;

	k0 = k;
	setup_usable(k + 1);
	setup_inverse();
	setup_walker();
	ppn = calloc(k + 1, sizeof(pp_n));
	pppp = calloc(k + 1, sizeof(pp_pp));
	pplist = (pp_pp**)NULL;
	pplistsize = 0;
	pplistmax = 0;
	ZINIT(&z_no_previous, "z_no_previous");
	mpz_set_si(z_no_previous, -1);
	
	for (i = 1; i <= k; ++i) {
		ppi = &ppn[i];
		ppi->n = i;
		ppi->pp = greatest_prime_power(i, &ppi->p);
		ppi->d = i / ppi->pp;
		ppi->usable = is_usable(ppi->p, ppi->pp, k);
		if (ppi->usable)
			pp_save_r(ppi);
	}
}

void pp_free(pp_pp* pp) {
	int i;
	for (i = 0; i < pp->valmax; ++i)
		ZCLEAR(&pp->value[i].value, "pp_%d.value[%d]", pp->pp, i);
	free(pp->value);
	QCLEAR(&pp->spare, "pp_%d.spare", pp->pp);
	QCLEAR(&pp->need, "pp_%d.need", pp->pp);
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
	free(ppn);
	for (i = 0; i < pplistsize; ++i)
		pplist[i]->wr = (walk_result*)NULL;
	for (i = 0; i <= k0; ++i)
		if (pppp[i].pp)
			pp_free(&pppp[i]);
	free(pplist);
	free(pppp);
	teardown_walker();
	teardown_inverse();
	teardown_usable();
}

void pp_grow(pp_pp* pp, int size) {
	int i;
	if (pp->valmax < size) {
		pp->value = realloc(pp->value, size * sizeof(pp_value));
		for (i = pp->valmax; i < size; ++i)
			ZINIT(&pp->value[i].value, "pp_%d.value[%d]", pp->pp, i);
		pp->valmax = size;
	}
}

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

void pp_listsplice(int i) {
	pp_pp* pp = pplist[i];
	if (i < --pplistsize)
		memmove(&pplist[i], &pplist[i + 1], (pplistsize - i) * sizeof(pplist[0]));
	/* Cannot free pp, the values are still referenced */
}

/*
 * Sort the values in a PP structure, given that only the last-added value
 * may be out of order.
 * Input:
 *   pp_pp* pp: pointer to the PP structure to sort.
 * Returns:
 *   Nothing.
 * Caveats:
 *   This copies a pp_value structure (which contains mpz_t) both by
 * structure assignment and by memmove. I don't know if this is safe to do.
 *   No guarantees are made about a stable sort: pp_value structures with
 * equal values may be sorted in any order.
 */
void pp_insert_last(pp_pp* pp) {
	int size = pp->valsize - 1;
	pp_value v = pp->value[size];
	int low, high, med, sign;

	/* invariant: pp->value[low].value <= v.value <= pp->value[high].value */
	low = -1;
	high = size + 1;
	while (low + 1 < high) {
		med = (low + high) >> 1;
		sign = mpz_cmp(v.value, pp->value[med].value);
		if (sign < 0)
			low = med;
		else
			high = med;
	}
	if (high + 1 < size) {
		memmove(&pp->value[high + 1], &pp->value[high], (size - high) * sizeof(pp_value));
		pp->value[high] = v;
	}
}

/*
 * Store a value in a PP structure.
 * Input:
 *   pp_pp* pp: pointer to the PP structure to store in
 *   mpq_t value: the value to store
 *   int parent: the parent to refer to
 * Returns:
 *   Nothing.
 * Effect:
 *   The specified value is stored in the specified PP structure, with a
 * reference to the specified parent. Internal variables are updated
 * accordingly.
 *
 * For a given value z, the caller is expected to have determined the
 * greatest prime power p^k dividing z. In this case, the pp structure
 * should be that for p^k, and the value passed should be the rational
 * z * p.
 *
 * The parent is used solely for determining a final solution vector.
 * If the value is a source value (i.e. 1/r representing the integer
 * r in the solution vector), the parent should be set to r. If not,
 * the parent is assumed to be a fully resolved PP structure storing
 * a walk result that specifies the actual vector contribution (possibly
 * recursively).
 */
void pp_save_any(pp_pp* pp, mpq_t value, int parent) {
	int i;
	mpz_t g, raise;
	pp_value* v;

	ZINIT(&g, "pp_save_any gcd");
	ZINIT(&raise, "pp_save_any raise");
	mpz_gcd(g, mpq_denref(value), pp->denominator);
	mpz_divexact(raise, mpq_denref(value), g);
	if (mpz_cmp_ui(raise, 1) > 0) {
		mpz_mul(pp->denominator, pp->denominator, raise);
		pp->invdenom = invfast(mod_ui(pp->denominator, pp->p), pp->p);
		mpz_mul(pp->total, pp->total, raise);
		for (i = 0; i < pp->valsize; ++i)
			mpz_mul(pp->value[i].value, pp->value[i].value, raise);
	}
	ZCLEAR(&raise, "pp_save_any raise");
	ZCLEAR(&g, "pp_save_any g");
	pp_grow(pp, pp->valsize + 1);
	v = &pp->value[pp->valsize++];
	mpz_mul(v->value, mpq_numref(value), pp->denominator);
	mpz_divexact(v->value, v->value, mpq_denref(value));
	v->parent = parent;
	v->inv = (pp->invdenom * mod_ui(v->value, pp->p)) % pp->p;
	mpz_add(pp->total, pp->total, v->value);
	pp->invtotal = (pp->invtotal + v->inv) % pp->p;

	pp_insert_last(pp);
}

void pp_save_r(pp_n* ppi) {
	pp_pp* pp = &pppp[ppi->pp];
	pp_value* v;
	int i, raise;
	mpq_t q;
	if (!pp->p) {
		pp->p = ppi->p;
		pp->pp = ppi->pp;
		pp_grow(pp, MINPPSET);
		if (pp->p == pp->pp)
			inverse_table(pp->p);
		pp_pushlist(pp);
		ZINIT(&pp->total, "pp_%d.total", pp->pp);
		ZINIT(&pp->denominator, "pp_%d.denominator", pp->pp);
		ZINIT(&pp->min_discard, "pp_%d.min_discard", pp->pp);
		QINIT(&pp->need, "pp_%d.need", pp->pp);
		QINIT(&pp->spare, "pp_%d.spare", pp->pp);
		pp->invdenom = 1;
		mpz_set_ui(pp->denominator, 1);
	}
	QINIT(&q, "pp_save_r temp");
	mpq_set_ui(q, 1, (unsigned int)ppi->d);
	pp_save_any(pp, q, ppi->n);
	QCLEAR(&q, "pp_save_r temp");
}

void pp_save_w(walk_result* wr, pp_pp* from) {
	mpq_t actual;
	int pp_effective = from->pp;
	int pp_applied = 1;
	int pp_p = from->p;
	int q, to, shared;
	pp_pp* pp;
	pp_value* v;

Dprintf("save wr: %Zd / %Zd [%d] from %d (%d)\n", wr->discard, from->denominator, wr->vec[0], from->pp, from->p);
	QINIT(&actual, "pp_save_w temp");
	from->wr = wr;
	mpz_sub(mpq_numref(actual), from->total, wr->discard);
	mpz_mul_ui(mpq_denref(actual), from->denominator, from->pp);
	mpq_canonicalize(actual);

	to = z_greatest_prime_power(mpq_denref(actual), (int*)NULL);
	mpz_divexact_ui(mpq_denref(actual), mpq_denref(actual), to);
	pp_save_any(&pppp[to], actual, from->pp);
	QCLEAR(&actual, "pp_save_w temp");
}

void pp_study(int target) {
	int i, dest;
	pp_pp *pp, *top;
	walker* w;
	walk_result *wr1, *wr2;
	mpz_t maxtotal;
	mpq_t spare;

	ZINIT(&maxtotal, "pp_study maxtotal");
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (i == 0 || pp->pp * pplist[i - 1]->pp > k0) {
			/* we are only interested in subsets that keep at least one
			 * element, so we set the limit to (total - 1)
			 */
			mpz_sub_ui(maxtotal, pp->total, 1);
			w = new_walker(pp, maxtotal, pp->invtotal);
			wr1 = walker_find(w, z_no_previous);
			if (!wr1) {
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			wr1 = wr_clone(w, wr1);
			wr2 = walker_find(w, wr1->discard);
			if (!wr2) {
				pp_save_w(wr1, pp);
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			mpz_set(pp->min_discard, wr1->discard);
			wr_clone_free(wr1);
			delete_walker(w);
		}
		dest = z_greatest_prime_power(pp->denominator, (int*)NULL);
		pppp[dest].depend = 1;
	}
	ZCLEAR(&maxtotal, "pp_study maxtotal");

	QINIT(&spare, "pp_study spare");
	top = pplist[0];
	mpq_set_ui(top->need, target, 1);
	mpq_set_si(top->spare, -target, 1);
	for (i = 0; i < pplistsize; ++i) {
		int g;

		pp = pplist[i];
		mpz_sub(mpq_numref(spare), pp->total, pp->min_discard);
		mpz_mul_ui(mpq_denref(spare), pp->denominator, pp->pp);
		mpq_canonicalize(spare);
		mpq_add(top->spare, top->spare, spare);
	}
	QCLEAR(&spare, "pp_study spare");
Dprintf("study: spare = %Qd\n", top->spare);
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
	int i, q, keepall;
	pp_value* v;
	mpz_t prod;

	ZINIT(&prod, "pp_setvec prod");
	for (i = 0; i < pp->valsize; ++i) {
		if (invec && pp_testbit(invec, i)) 
			continue;
		v = &pp->value[i];
		q = v->parent / pp->pp;
		keepall = 1;
		if (q * pp->pp != v->parent) {
			keepall = 0;
		} else {
			mpz_mul_ui(prod, v->value, q);
			if (mpz_cmp(prod, pp->denominator) != 0)
				keepall = 0;
		}
		if (keepall) {
			pp_setbit(vec, v->parent);
		} else {
			pp_setvec(vec, &pppp[v->parent], &pppp[v->parent].wr->vec[0]);
		}
	}
	ZCLEAR(&prod, "pp_setvec prod");
}

int pp_solution(int index) {
	pp_pp *pp, *cur;
	int* v;
	int i, j, r;

	cur = pplist[index];
	if (mpz_cmp_ui(mpq_numref(cur->spare), 1) > 0)
		return 0;
	r = mpz_get_ui(mpq_numref(cur->spare));
	if (r == 1 && mpz_cmp_ui(mpq_denref(cur->spare), k0) > 0)
		return 0;
	r = (r == 0) ? 0 : mpz_get_ui(mpq_denref(cur->spare));

	/* probable solution: spare == 0 or spare == 1/r, 1 <= r <= k0 */
	gmp_printf("probable solution, spare = %Qd (need = %Qd)\n", cur->spare, cur->need);
	v = calloc((k0 + 32) >> 5, sizeof(int));
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (i >= index && pp->depend) {
			pp_setvec(v, pp, (int*)NULL);
		} else {
			pp_setvec(v, pp, &pp->wr->vec[0]);
		}
	}

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
	printf("Show levels %d..%d at %d of %d\n", start, end, level, pplistsize);
	for (j = start; j <= end; ++j) {
		pp = pplist[j];
		printf("%d:%d/%d ", pp->pp, pp->wrnum, pp->wrcount);
	}
	printf("\n\n");
}

int pp_find(int target) {
	pp_pp *pp, *nextpp;
	int i, g, inv;
	int success = 0;
	int invsum;
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
			invsum = pp->invtotal;

Dprintf("at pp_%d spare = %Qd\n", pp->pp, limit);
			mpz_mul_ui(z, pp->denominator, pp->pp);
			if (!pp->depend && mpz_sgn(pp->min_discard) != 0) {
				/* (limit_r) += min_discard / denominator / pp */
				mpz_set(mpq_numref(q), pp->min_discard);
				mpz_set(mpq_denref(q), z);
				mpq_canonicalize(q);
				mpq_add(limit, limit, q);
			}
			/* actual limit = floor(limit_num * denominator * pp / limit_den) */
			mpq_mul_z(limit, limit, z);
			mpz_fdiv_q(mpq_numref(limit), mpq_numref(limit), mpq_denref(limit));
			/* limit is now denormal, but we won't be using it again */

			if (pp->depend) {
				/* must additionally discard -need_num * inv(need_den) */
				mpz_fdiv_qr_ui(z, zr, mpq_denref(pp->need), pp->pp);
				if (mpz_sgn(zr) == 0) {
					inv = invfast(mod_ui(z, pp->p), pp->p);
					mpz_mul_ui(z, mpq_numref(pp->need), inv);
					inv = pp->p - mod_ui(z, pp->p);
					invsum = (invsum + inv) % pp->p;
				}
			}
Dprintf("effective limit is %Zd, invsum = %d\n", mpq_numref(limit), invsum);

			pp->w = new_walker(pp, mpq_numref(limit), invsum);
			pp->wrnum = 0;
		}

		pp->wr = walker_find(pp->w, pp->wr ? pp->wr->discard : z_no_previous);
		if (!pp->wr) {
			delete_walker(pp->w);
			pp->w = (walker*)NULL;
			pp->wrcount = pp->wrnum;
			--i;
Dprintf("no lines for %d, recurse back to %d\n", pp->pp, (i >= 0) ? pplist[i]->pp : -1);
			continue;
		}
		++pp->wrnum;
Dprintf("pp_%d walker found discard %Zd/%Zd\n", pp->pp, pp->wr->discard, pp->denominator);

		++i;
		if (i >= pplistsize) {
			fprintf(stderr, "overflow\n");
			exit(1);
		}
		nextpp = pplist[i];

		mpz_sub(mpq_numref(q), pp->wr->discard, pp->min_discard);
		mpz_mul_ui(mpq_denref(q), pp->denominator, pp->pp);
		mpq_canonicalize(q);
		mpq_sub(nextpp->spare, pp->spare, q);

		mpz_sub(mpq_numref(q), pp->total, pp->wr->discard);
		mpz_mul_ui(mpq_denref(q), pp->denominator, pp->pp);
		mpq_canonicalize(q);
		mpq_sub(nextpp->need, pp->need, q);

Dprintf("new spare %Qd; new need %Qd\n", nextpp->spare, nextpp->need);

		if (mpq_sgn(nextpp->spare) < 0) {
			/* we've overrun */
			--i;
			continue;
		}
		nextpp->w = (walker*)NULL;
		if (nextpp->depend && pp_solution(i)) {
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
