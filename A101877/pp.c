#include "pp.h"
#include "usable.h"
#include "prime.h"
#include "inverse.h"
#include "walker.h"
#include <stdlib.h>
#include <stdio.h>

pp_n* ppn = (pp_n*)NULL;
pp_pp* pppp = (pp_pp*)NULL;
pp_list* pplist = (pp_list*)NULL;
int pplistsize = 0;
int pplistmax = 0;
int k0;

void init_pp(int k) {
	int i, d, q;
	pp_n* ppi;

	k0 = k;
	init_usable(k + 1);
	setup_inverse();
	walker_init();
	ppn = calloc(k + 1, sizeof(pp_n));
	pppp = calloc(k + 1, sizeof(pp_pp));
	
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

void pp_grow(pp_pp* pp, int size) {
	if (pp->valmax < size) {
		pp->valmax = size;
		pp->value = realloc(pp->value, size * sizeof(pp_value));
	}
}

void pp_pushlist(pp_pp* pp) {
	if (pplistsize + 1 >= pplistmax) {
		pplistmax = pplistmax * 3 / 2;
		if (pplistmax < MINPPSET)
			pplistmax = MINPPSET;
		pplist = realloc(pplist, pplistmax * sizeof(pp_list));
	}
	pplist[pplistsize++].pp = pp;
}

void pp_free(pp_pp* pp) {
	free(pp->value);
	memset(pp, 0, sizeof(pp_pp));
}

void pp_listsplice(int i) {
	pp_pp* pp = pplist[i].pp;
	/* FIXME: do we need to keep order, or can we just move down the last? */
	if (i < --pplistsize)
		memmove(&pplist[i], &pplist[i + 1], (pplistsize - i) * sizeof(pp_list));
	/* Cannot free pp, since parentvec elsewhere may hold refs */
	/* pp_free(pp); */
}

void pp_insert_last(pp_pp* pp) {
	int size = pp->valsize - 1;
	pp_value v = pp->value[size];
	int low, high, med;

	/* invariant: pp->value[low].value <= v.value <= pp->value[high].value */
	low = -1;
	high = size + 1;
	while (low + 1 < high) {
		med = (low + high) >> 1;
		if (v.value < pp->value[med].value)
			low = med;
		else
			high = med;
	}
	if (high + 1 < size) {
		memmove(&pp->value[high + 1], &pp->value[high], (size - high) * sizeof(pp_value));
		pp->value[high] = v;
	}
}

void pp_save_any(pp_pp* pp, int vnum, int vden, int parent) {
	int i, raise;
	pp_value* v;

	raise = vden / gcd(vden, pp->denominator);
	if (raise > 1) {
		pp->denominator *= raise;
		pp->invdenom = invfast(pp->denominator % pp->p, pp->p);
		pp->total *= raise;
		for (i = 0; i < pp->valsize; ++i)
			pp->value[i].value *= raise;
	}
	pp_grow(pp, pp->valsize + 1);
	v = &pp->value[pp->valsize++];
	v->value = vnum * (pp->denominator / vden);
	v->parent = parent;
	v->inv = (pp->invdenom * (v->value % pp->p)) % pp->p;
	pp->total += v->value;
	pp->invtotal = (pp->invtotal + v->inv) % pp->p;

	pp_insert_last(pp);
}

void pp_save_r(pp_n* ppi) {
	pp_pp* pp = &pppp[ppi->pp];
	pp_value* v;
	int i, raise;
	if (!pp->p) {
		pp->p = ppi->p;
		pp->pp = ppi->pp;
		pp_grow(pp, MINPPSET);
		if (pp->p == pp->pp)
			inverse_table(pp->p);
		pp_pushlist(pp);
		pp->total = 0;
		pp->denominator = 1;
		pp->invdenom = 1;
	}
	pp_save_any(pp, 1, ppi->d, ppi->n);
}

void pp_save_w(walk_result* wr, pp_pp* from) {
	int actual_num;
	int actual_den;
	int pp_effective = from->pp;
	int pp_applied = 1;
	int pp_p = from->p;
	int q, to, shared;
	pp_pp* pp;
	pp_value* v;

printf("save wr: %d / %d [%d] from %d (%d)\n", wr->discard, from->denominator, wr->vec[0], from->pp, from->p);
	from->wr = wr;
	actual_num = from->total - wr->discard;
	actual_den = from->denominator;
	shared = gcd(actual_num, actual_den);
	if (shared > 1) {
		actual_num /= shared;
		actual_den /= shared;
	}
printf("actual %d / %d (shared %d)\n", actual_num, actual_den, shared);
	while (pp_applied != pp_effective) {
		q = actual_num / pp_p;
		if (q * pp_p != actual_num)
			break;
		actual_num = q;
		pp_applied *= pp_p;
	}
	if (pp_applied != pp_effective)
		actual_den *= pp_effective / pp_applied;
printf("actual %d / %d\n", actual_num, actual_den);

	to = greatest_prime_power(actual_den, (int*)NULL);
	actual_den /= to;
	pp_save_any(&pppp[to], actual_num, actual_den, from->pp);
}

void pp_study(void) {
	int i, target;
	pp_pp* pp;
	walker* w;
	walk_result *wr1, *wr2;
	for (i = pplistsize - 1; i >= 0; --i) {
		pp = pplist[i].pp;
		if (!pp->depend) {
			w = new_walker(pp, pp->total - 1, pp->invtotal);
			wr1 = walker_find(w);
			if (!wr1) {
printf("discard %d, no solutions\n", pp->pp);
				delete_walker(w);
				pp_listsplice(i);
				continue;
			}
			wr1 = wr_clone(w, wr1);
			wr2 = walker_find(w);
			if (!wr2) {
				pp_save_w(wr1, pp);
				delete_walker(w);
				pp_listsplice(i);
				continue;
			}
			pp->min_discard = wr1->discard;
			free(wr1);
			delete_walker(w);
		}
		target = greatest_prime_power(pp->denominator, (int*)NULL);
		pppp[target].depend = 1;
printf("for %d set first dependent %d (denominator %d)\n", pp->pp, target, pp->denominator);
	}
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
	int i, q;
	pp_value* v;

	for (i = 0; i < pp->valsize; ++i) {
		if (invec && pp_testbit(invec, i)) 
			continue;
		v = &pp->value[i];
		q = v->parent / pp->pp;
		if (q * pp->pp == v->parent && v->value * q == pp->denominator) {
			pp_setbit(vec, v->parent);
		} else {
			pp_setvec(vec, &pppp[v->parent], &pppp[v->parent].wr->vec[0]);
		}
	}
}

int pp_solution(int index) {
	pp_list* ppl = &pplist[index];
	pp_pp* pp;
	int* v;
	int i, j;

	if (ppl->spare_num > 1)
		return 0;
	if (ppl->spare_num == 1 && ppl->spare_den > k0)
		return 0;
	/* probable solution: spare == 0 or spare == 1/r, 1 <= r <= k0 */
	printf("probable solution, spare = %d/%d (need = %d/%d)\n", ppl->spare_num, ppl->spare_den, ppl->need_num, ppl->need_den);
	v = calloc((k0 + 31) >> 5, sizeof(int));
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i].pp;
		if (i <= index && pp->depend) {
			pp_setvec(v, pp, (int*)NULL);
		} else {
			pp_setvec(v, pp, &pp->wr->vec[0]);
		}
	}

	if (pp_testbit(v, ppl->spare_den)) {
		/* exact solution */
		v[ppl->spare_den >> 5] &= ~(1 << (ppl->spare_den & 31));
	} else {
		/* inexact, maybe fixable */
		printf("inexact(%d): ", ppl->spare_den);
	}
	for (i = 0; i <= k0; ++i)
		if (pp_testbit(v, i))
			printf("%d ", i);
	printf("\n");
	return 1;
}

void pp_find(int target) {
	pp_list* ppl = &pplist[pplistsize - 1];
	pp_pp* pp;
	int i, num, den, g;

	ppl->need_num = target;
	ppl->need_den = 1;
	ppl->spare_num = -target;
	ppl->spare_den = 1;
	for (i = pplistsize - 1; i >= 0; --i) {
		pp = pplist[i].pp;
		num = pp->total - (pp->depend ? 0 : pp->min_discard);
		den = pp->denominator * pp->pp;
		/* rat_add(&ppl->spare_num, &ppl->spare_den, max, den); */
		g = gcd(num, den);
		if (g > 1) {
			num /= g;
			den /= g;
		}
		g = gcd(ppl->spare_den, den);
		den /= g;
		ppl->spare_num = ppl->spare_num * den + num * (ppl->spare_den / g);
		ppl->spare_den *= den;
		g = gcd(abs(ppl->spare_num), ppl->spare_den);
		if (g > 1) {
			ppl->spare_num /= g;
			ppl->spare_den /= g;
		}
printf("spare at pp_%d = %d/%d\n", pp->pp, ppl->spare_num, ppl->spare_den);
	}

	i = pplistsize - 1;
	pplist[i].w = (walker*)NULL;
	while (i < pplistsize) {
		ppl = &pplist[i];
		pp = ppl->pp;
		if (!ppl->w) {
			int limit_num = ppl->spare_num;
			int limit_den = ppl->spare_den;
			int invsum = pp->invtotal;

printf("at pp_%d spare = %d / %d\n", pp->pp, limit_num, limit_den);
			if (!pp->depend && pp->min_discard) {
				/* (limit_r) -= min_discard / denominator / pp */
				int num = pp->min_discard;
				int den = pp->denominator * pp->pp;
				g = gcd(num, den);
				if (g > 1) {
					num /= g;
					den /= g;
				}
				g = gcd(limit_den, den);
				limit_num = limit_num * (den / g) + num * (limit_den / g);
				limit_den *= den / g;
printf("with min_discard = %d / %d / %d, actual spare = %d / %d\n", pp->min_discard, pp->denominator, pp->pp, limit_num, limit_den);
			}
			/* effective limit = limit * denominator * pp */
			g = gcd(limit_den, pp->denominator);
			limit_den /= g;
			limit_num *= pp->denominator / g;
			limit_num = limit_num * pp->pp / limit_den;
printf("effective limit is %d / %d\n", limit_num, pp->denominator);

printf("at pp_%d invsum = %d\n", pp->pp, invsum);
			if (pp->depend) {
				/* must additionally discard -need_num * inv(need_den) */
				int need_num = ppl->need_num;
				int need_den = ppl->need_den / pp->pp;
				if (need_den * pp->pp == ppl->need_den) {
					int inv = invfast(need_den % pp->p, pp->p);
					inv = pp->p - ((need_num * inv) % pp->p);
					invsum = (invsum + inv) % pp->p;
				}
printf("with need = %d / %d, actual invsum = %d\n", ppl->need_num, ppl->need_den, invsum);
			}

			ppl->w = new_walker(pp, limit_num, invsum);
		}

		pp->wr = walker_find(ppl->w);
		if (!pp->wr) {
			delete_walker(ppl->w);
			++i;
printf("no lines for %d, recurse back to %d\n", pp->pp, pplist[i].pp->pp);
			continue;
		}
		--i;
if (i < 0) {
	fprintf(stderr, "underflow\n");
	exit(1);
}
		{
			int num = pp->wr->discard - (pp->depend ? 0 : pp->min_discard);
			int den = pp->denominator * pp->pp;
			g = gcd(ppl->spare_den, den);
			pplist[i].spare_num = ppl->spare_num * (den / g) - num * (ppl->spare_den / g);
			pplist[i].spare_den = ppl->spare_den * (den / g);
			g = gcd(pplist[i].spare_num, pplist[i].spare_den);
			pplist[i].spare_num /= g;
			pplist[i].spare_den /= g;
printf("new spare %d / %d = old spare %d / %d - (new discard %d - min discard %d) / %d / %d\n", pplist[i].spare_num, pplist[i].spare_den, ppl->spare_num, ppl->spare_den, pp->wr->discard, (pp->depend ? 0 : pp->min_discard), pp->denominator, pp->pp);
		}

		{
			int num = pp->total - pp->wr->discard;
			int den = pp->denominator * pp->pp;
			g = gcd(ppl->need_den, den);
			pplist[i].need_num = ppl->need_num * (den / g) - num * (ppl->need_den / g);
			pplist[i].need_den = ppl->need_den * (den / g);
			g = gcd(pplist[i].need_num, pplist[i].need_den);
			pplist[i].need_num /= g;
			pplist[i].need_den /= g;
printf("new need %d / %d = old need %d / %d - (total %d - discard %d) / %d / %d\n", pplist[i].need_num, pplist[i].need_den, ppl->need_num, ppl->need_den, pp->total, pp->wr->discard, pp->denominator, pp->pp);
		}
		pplist[i].w = (walker*)NULL;

		if (pplist[i].need_num < 0) {
			/* we've overrun */
			++i;
			continue;
		}
		if (pplist[i].pp->depend && pp_solution(i))
			return;
	}
	printf("no solution found\n");
}
