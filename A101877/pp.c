#include "pp.h"
#include "usable.h"
#include "prime.h"
#include "inverse.h"
#include "walker.h"
#include <stdlib.h>
#include <stdio.h>

pp_n* ppn = (pp_n*)NULL;
pp_pp* pppp = (pp_pp*)NULL;
pp_pp** pplist;
int pplistsize;
int pplistmax;
int k0;

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

void teardown_pp(void) {
	int i;
	pp_pp* pp;

	free(ppn);
	for (i = 0; i < pplistsize; ++i) {
		pplist[i]->wr = (walk_result*)NULL;
	}
	for (i = 0; i <= k0; ++i) {
		pp = &pppp[i];
		if (pp->value)
			free(pp->value);
		if (pp->wr)
			free(pp->wr);
	}
	free(pplist);
	free(pppp);
	teardown_walker();
	teardown_inverse();
	teardown_usable();
}

void pp_grow(pp_pp* pp, int size) {
	if (pp->valmax < size) {
		pp->valmax = size;
		pp->value = realloc(pp->value, size * sizeof(pp_value));
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

void pp_free(pp_pp* pp) {
	free(pp->value);
	memset(pp, 0, sizeof(pp_pp));
}

void pp_listsplice(int i) {
	pp_pp* pp = pplist[i];
	if (i < --pplistsize)
		memmove(&pplist[i], &pplist[i + 1], (pplistsize - i) * sizeof(pplist[0]));
	/* Cannot free pp, the values are still referenced */
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
	while (pp_applied != pp_effective) {
		q = actual_num / pp_p;
		if (q * pp_p != actual_num)
			break;
		actual_num = q;
		pp_applied *= pp_p;
	}
	if (pp_applied != pp_effective)
		actual_den *= pp_effective / pp_applied;

	to = greatest_prime_power(actual_den, (int*)NULL);
	actual_den /= to;
	pp_save_any(&pppp[to], actual_num, actual_den, from->pp);
}

void pp_study(int target) {
	int i, dest;
	pp_pp *pp, *top;
	walker* w;
	walk_result *wr1, *wr2;

	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (i == 0 || pp->pp * pplist[i - 1]->pp > k0) {
			w = new_walker(pp, pp->total - 1, pp->invtotal);
			wr1 = walker_find(w);
			if (!wr1) {
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			wr1 = wr_clone(w, wr1);
			wr2 = walker_find(w);
			if (!wr2) {
				pp_save_w(wr1, pp);
				delete_walker(w);
				pp_listsplice(i);
				--i;
				continue;
			}
			pp->min_discard = wr1->discard;
			free(wr1);
			delete_walker(w);
		}
		dest = greatest_prime_power(pp->denominator, (int*)NULL);
		pppp[dest].depend = 1;
	}

	top = pplist[0];
	top->need_num = target;
	top->need_den = 1;
	top->spare_num = -target;
	top->spare_den = 1;
	for (i = 0; i < pplistsize; ++i) {
		int num, den, g;

		pp = pplist[i];
		num = pp->total - pp->min_discard;
		den = pp->denominator * pp->pp;
		/* rat_add(&top->spare_num, &top->spare_den, max, den); */
		g = gcd(num, den);
		if (g > 1) {
			num /= g;
			den /= g;
		}
		g = gcd(top->spare_den, den);
		den /= g;
		top->spare_num = top->spare_num * den + num * (top->spare_den / g);
		top->spare_den *= den;
		g = gcd(abs(top->spare_num), top->spare_den);
		if (g > 1) {
			top->spare_num /= g;
			top->spare_den /= g;
		}
	}
printf("study: spare = %d/%d\n", top->spare_num, top->spare_den);
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
	pp_pp *pp, *cur;
	int* v;
	int i, j;

	cur = pplist[index];
	if (cur->spare_num > 1)
		return 0;
	if (cur->spare_num == 1 && cur->spare_den > k0)
		return 0;
	/* probable solution: spare == 0 or spare == 1/r, 1 <= r <= k0 */
	printf("probable solution, spare = %d/%d (need = %d/%d)\n", cur->spare_num, cur->spare_den, cur->need_num, cur->need_den);
	v = calloc((k0 + 31) >> 5, sizeof(int));
	for (i = 0; i < pplistsize; ++i) {
		pp = pplist[i];
		if (i >= index && pp->depend) {
			pp_setvec(v, pp, (int*)NULL);
		} else {
			pp_setvec(v, pp, &pp->wr->vec[0]);
		}
	}

	if (pp_testbit(v, cur->spare_den)) {
		/* exact solution */
		v[cur->spare_den >> 5] &= ~(1 << (cur->spare_den & 31));
	} else {
		/* inexact, maybe fixable */
		printf("inexact(%d): ", cur->spare_den);
	}
	for (i = 0; i <= k0; ++i)
		if (pp_testbit(v, i))
			printf("%d ", i);
	printf("\n");
	free(v);
	return 1;
}

int pp_find(int target) {
	pp_pp *pp, *nextpp;
	int i, num, den, g, inv;
	int success = 0;

	pp_study(target);

	i = 0;
	pplist[i]->w = (walker*)NULL;
	if (pplist[i]->spare_num < 0) {
		printf("Optimizer finds no solution is possible.\n");
		return 0;
	}
	while (i >= 0) {
		pp = pplist[i];
		if (!pp->w) {
			int limit_num = pp->spare_num;
			int limit_den = pp->spare_den;
			int invsum = pp->invtotal;

printf("at pp_%d spare = %d / %d\n", pp->pp, limit_num, limit_den);
			if (!pp->depend && pp->min_discard) {
				/* (limit_r) -= min_discard / denominator / pp */
				num = pp->min_discard;
				den = pp->denominator * pp->pp;
				g = gcd(num, den);
				if (g > 1) {
					num /= g;
					den /= g;
				}
				g = gcd(limit_den, den);
				limit_num = limit_num * (den / g) + num * (limit_den / g);
				limit_den *= den / g;
			}
			/* effective limit = limit * denominator * pp */
			g = gcd(limit_den, pp->denominator);
			limit_den /= g;
			limit_num *= pp->denominator / g;
			limit_num = limit_num * pp->pp / limit_den;

			if (pp->depend) {
				/* must additionally discard -need_num * inv(need_den) */
				num = pp->need_num;
				den = pp->need_den / pp->pp;
				if (den * pp->pp == pp->need_den) {
					inv = invfast(den % pp->p, pp->p);
					inv = pp->p - ((num * inv) % pp->p);
					invsum = (invsum + inv) % pp->p;
				}
			}
printf("effective limit is %d / %d, invsum = %d\n", limit_num, pp->denominator, invsum);

			pp->w = new_walker(pp, limit_num, invsum);
		}

		pp->wr = walker_find(pp->w);
		if (!pp->wr) {
			delete_walker(pp->w);
			pp->w = (walker*)NULL;
			--i;
printf("no lines for %d, recurse back to %d\n", pp->pp, (i >= 0) ? pplist[i]->pp : -1);
			continue;
		}

		++i;
if (i >= pplistsize) {
	fprintf(stderr, "overflow\n");
	exit(1);
}
		nextpp = pplist[i];

		num = pp->wr->discard - (pp->depend ? 0 : pp->min_discard);
		den = pp->denominator * pp->pp;
		g = gcd(pp->spare_den, den);
		nextpp->spare_num = pp->spare_num * (den / g) - num * (pp->spare_den / g);
		nextpp->spare_den = pp->spare_den * (den / g);
		g = gcd(nextpp->spare_num, nextpp->spare_den);
		nextpp->spare_num /= g;
		nextpp->spare_den /= g;

		num = pp->total - pp->wr->discard;
		den = pp->denominator * pp->pp;
		g = gcd(pp->need_den, den);
		nextpp->need_num = pp->need_num * (den / g) - num * (pp->need_den / g);
		nextpp->need_den = pp->need_den * (den / g);
		g = gcd(nextpp->need_num, nextpp->need_den);
		nextpp->need_num /= g;
		nextpp->need_den /= g;
printf("new spare %d / %d; new need %d / %d\n", nextpp->spare_num, nextpp->spare_den, nextpp->need_num, nextpp->need_den);

		nextpp->w = (walker*)NULL;
		if (nextpp->need_num < 0) {
			/* we've overrun */
			--i;
			continue;
		}
		if (nextpp->depend && pp_solution(i)) {
			success = 1;
			break;
		}
	}
	for (i = 0; i < pplistsize; ++i) {
		if (pplist[i]->w)
			delete_walker(pplist[i]->w);
	}
	if (!success)
		printf("no solution found\n");
	return success;
}
