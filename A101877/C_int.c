#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <string.h>

#include "vector.h"
#include "match.h"
#include "total.h"

/* in-place add_ui preserves canonicalisation */
#define mpqi_add_ui(q, u) mpz_addmul_ui(mpq_numref(q), mpq_denref(q), u)

#ifdef DEBUG_GMP_LEAK
#define QINIT(q, format, args...) ({ \
	mpq_t *_q = q; \
	mpq_init(*_q); \
	fprintf(stderr, "%p init: " format "\n", _q, args); \
})
#define QCLEAR(q, format, args...) ({ \
	mpq_t *_q = q; \
	mpq_clear(*_q); \
	fprintf(stderr, "%p clear: " format "\n", _q, args); \
})
#else
#define QINIT(q, format, args...) mpq_init(*q);
#define QCLEAR(q, format, args...) mpq_clear(*q);
#endif

/* resizing arrays */
typedef struct array_s {
	void* array;		/* pointer to currently malloced area */
	ulong count;		/* number of used entries */
	ulong space;		/* number of malloced entries */
} array_t;

/* pp cache */
typedef struct cache_s {
	match_t *lfm;	/* growing array of match_t */
	array_t nfm;	/* list of numbers used */
	int havelfm;	/* boolean */
	int complete;	/* boolean */
} cache_t;

/* data per prime power */
typedef struct pp_s {
	ulong pp;		/* the prime power represented here */
	ulong p;		/* the prime */
	int count;		/* number of values accounted */
	array_t val;	/* list of (ulong) values accounted */
	array_t rec;	/* list of (ulong) modular inverses */
	array_t span;	/* list of (vector_t *) bitvectors of reachable mod values */
	total_t *totals;	/* triangle of running totals */
	int index;		/* index of this struct in known[] */
	mpq_t sum;		/* sum of reciprocals */
	ulong modsum;	/* sum of self.rec[] */
	array_t match;	/* array of p dynamic arrays of (match_t *) */
	array_t cache;	/* (i mod p) array of growing array of (cache_t) */
	vector_t* complete;	/* p-bit cache of known completions */
	int cansolve;	/* boolean: true if solvable for mod = 0 */
} pp_t;

/* array of pp_t structures */
array_t known_pp;

/* basic globals */
ulong n;			/* A101877() index, and target sum */
ulong k;			/* possible value of A101877(n) being tested */
mpq_t H_k;			/* sum_1^k{1/i} */
vector_t *solution;	/* bit vector constructed by pp_span */

/* for timings */
long clock_tick;

/* scratch variables */
mpq_t q0;
mpz_t z0, z1, z2;
mpz_t mi_r, mi_p, mi_k;	/* mod_invert */
mpq_t rlimit;			/* pp_consider: 1/k */
array_t consider_sum;	/* pp_consider: mpz_t of sums by pp */
array_t consider_kept;	/* pp_consider: mpz_t of total kept by pp */
array_t consider_tried;	/* pp_consider: match objects of current trial by pp */
int diagnose_retry;		/* pp_diagnose: pp index above which to show */
int diagnose_len1;		/* pp_diagnose: lengths written */
int diagnose_len2;

/* factorisation */
array_t prime;			/* list of (ulong) primes */
array_t factor;			/* list giving prime p for each prime power p^k */

void array_resize(array_t* a, size_t size, size_t count) {
	if (a->space >= count) return;
	if (a->space == 0) {
		if (count < 10) count = 10;
		a->array = malloc(size * count);
		a->space = count;
	} else {
		if (count < a->space * 2) count = a->space * 2;
		a->array = realloc(a->array, size * count);
		a->space = count;
	}
}

/* return time so far in seconds */
double timing(void) {
	struct tms ttd;
	times(&ttd);
	return ((double)ttd.tms_utime) / clock_tick;
}

/* find and return next prime not yet in the prime array */
ulong nextprime(void) {
	ulong* p0 = (ulong *)prime.array;
	ulong p = p0[prime.count - 1];
	ulong d;
	int pc;

	while (1) {
		p += 2;
		pc = 0;
		while (1) {
			d = p0[pc++];
			if (d * d > p) {
				array_resize(&prime, sizeof(ulong), prime.count + 1);
				p0 = (ulong *)prime.array;
				p0[prime.count++] = p;
				array_resize(&factor, sizeof(ulong), p + 1);
				((ulong *)factor.array)[p] = p;
				return p;
			}
			if ((p % d) == 0) break;
		}
	}
}

/*
	return the greatest prime power dividing n0
	ensures along the way that the prime and factor arrays are extended
	far enough
*/
ulong gpp(ulong n0) {
	ulong best = 1;
	ulong* p0 = (ulong *)prime.array;
	int pc = 0, pz = prime.count;
	ulong n = n0, d, s;
	while (n > 1) {
		d = p0[pc++];
		if (d * d > n) break;
		if ((n % d) != 0) continue;
		s = 1;
		while ((n % d) == 0) {
			s *= d;
			n /= d;
		}
		if (best < s) best = s;
		if (best == n0) {
			array_resize(&factor, sizeof(ulong), best + 1);
			((ulong *)factor.array)[best] = d;
		}
	}
	if (n > 1) {
		/* remnant is prime */
		if (best < n) best = n;
		if (n > p0[pz - 1]) {
			while (nextprime() < n)
				;
		}
	}
	return best;
}

/*
	return inverse of k (mod p)
*/
ulong mod_invert(ulong k, ulong p) {
	mpz_set_ui(mi_p, p);
	mpz_set_ui(mi_k, k);
	mpz_invert(mi_r, mi_k, mi_p);
	return mpz_get_ui(mi_r);
}

void free_cache_array(array_t *ap) {
	int i;
	cache_t *cp = (cache_t *)ap->array;
	for (i = 0; i < ap->count; ++i) {
		if (cp[i].havelfm) free_match(cp[i].lfm);
		cp[i].nfm.count = cp[i].havelfm = cp[i].complete = 0;
		if (cp[i].nfm.space) free(cp[i].nfm.array);
		cp[i].nfm.space = 0;
	}
	ap->count = 0;
}

/* locate the existing pp structure for prime power gpp */
pp_t* pp_locate(ulong gpp) {
	ulong low = 0, high = known_pp.count, med;
	pp_t* pp0 = (pp_t *)known_pp.array;
	while (low < high) {
		med = (low + high) >> 1;
		if (pp0[med].pp < gpp)
			low = med + 1;
		else
			high = med;
	}
	return &pp0[low];
}

void pp_free(pp_t *self) {
	int i, j;
	vector_t **vp;
	array_t *map, *cap;
	match_t **mp;

	free(self->val.array);
	free(self->rec.array);
	vp = (vector_t **)self->span.array;
	for (i = 0; i < self->span.count; ++i) free_vector(vp[i]);
	free(self->span.array);
	free_total(self->totals);
	QCLEAR(&self->sum, "pp_free pp=%lu sum", self->pp);
	map = (array_t *)self->match.array;
	for (i = 0; i < self->p; ++i) {
		mp = (match_t **)map[i].array;
		for (j = 0; j < map[i].count; ++j) free_match(mp[j]);
		free(map[i].array);
	}
	free(self->match.array);
	cap = (array_t *)self->cache.array;
	for (i = 0; i < self->p; ++i) {
		free_cache_array(cap + i);
		free(cap[i].array);
	}
	free(self->cache.array);
	if (self->complete) free_vector(self->complete);
}

/* create a new pp structure for prime power gpp and return it */
pp_t* pp_new(ulong gpp) {
	pp_t* self;
	ulong p;
	array_t *match0;
	match_t **mp;
	vector_t *v;
	array_resize(&known_pp, sizeof(pp_t), known_pp.count + 1);

	self = ((pp_t *)known_pp.array) + known_pp.count;
	memset(self, 0, sizeof(pp_t));
	self->index = known_pp.count++;
	self->pp = gpp;
	p = self->p = ((ulong *)factor.array)[gpp];
	QINIT(&self->sum, "new gpp=%lu sum", gpp);

	array_resize(&self->span, sizeof(vector_t *), 1);
	v = ((vector_t **)self->span.array)[0] = new_vector(p);
	fill_vector(v, 0);
	set_vector(v, 0, 1);
	self->span.count = 1;

	self->totals = new_total();

	v = self->complete = new_vector(p);
	fill_vector(v, 0);
	/* mod 0 is initially completable */
	set_vector(v, 0, 1);

/* match is an array of p arrays of (match_t *) */
	array_resize(&self->match, sizeof(array_t), p);
	memset(self->match.array, 0, sizeof(array_t) * p);
	self->match.count = p;
/* match[0] is an array of (match_t *) */
	match0 = (array_t *)self->match.array;
	array_resize(match0, sizeof(match_t *), 1);
/* match[0][0] = { sum => 0, vec => NULL } */
	mp = (match_t **)match0->array;
	mp[0] = new_match();
	match0->count = 1;

	array_resize(&self->cache, sizeof(array_t), p);
	memset(self->cache.array, 0, sizeof(array_t) * p);

	return self;
}

void pp_undiagnose(int length) {
#ifndef DEBUG_GMP_LEAK
	while (length--) fprintf(stderr, "\x08 \x08");
#endif
}

void pp_diagnose(pp_t *self, int match_index) {
#ifndef DEBUG_GMP_LEAK
	int index = self->index;
	if (index - 1 > diagnose_retry && match_index != 0)
		diagnose_retry = index - 1;
	if (index > diagnose_retry && match_index == 0)
		return;
	if (index > diagnose_retry) {
		pp_undiagnose(diagnose_len1 + diagnose_len2);
		diagnose_len1 = fprintf(stderr, "%lu %d", self->pp, match_index);
		diagnose_len2 = 0;
	} else {
		pp_undiagnose(diagnose_len2);
		diagnose_len2 = fprintf(stderr, " %lu %d", self->pp, match_index);
	}
#endif
}

/* construct a (possible) solution set as a bit vector */
void pp_span(match_t **tried, int offset) {
	vector_t *v;
	int i, j;
	pp_t *pp = (pp_t *)known_pp.array;

	if (solution) free_vector(solution);
	solution = new_vector(k);
	fill_vector(solution, 1);
	for (i = offset; i < known_pp.count; ++i) {
		v = tried[i]->v;	/* v is NULL if keeping none */
		for (j = 0; j < pp[i].count; ++j)
			if (!v || !test_vector(v, j))
				set_vector(solution, ((ulong *)pp[i].val.array)[j], 0);
	}
}

/* print a found solution */
void pp_sol(void) {
	int i;

	printf("a(%lu) = %lu : { ", n, k);
	for (i = 1; i <= k; ++i)
		if (test_vector(solution, i))
			printf("%d ", i);
	printf("}\n");
}

/* we have a solution: find the set, and print it */
int pp_solution(match_t **tried, int offset) {
	pp_span(tried, offset);
	pp_sol();
	return 1;
}

/* we have a solution iff the set includes fly; if so, find and print the set */
int pp_maybe_solution(match_t **tried, int offset, ulong fly) {
	pp_span(tried, offset);
	if (!test_vector(solution, fly)) return 0;
	set_vector(solution, fly, 0);
	pp_sol();
	return 1;
}

/*
	given an array of (match_t *) sorted by sum, find the index of the first
	match that has sum >= boundary
*/
int pp_find_chop(array_t *array, mpq_t boundary) {
	int min, max, med;
	match_t **match = (match_t **)array->array;
	min = -1;
	max = array->count;
	while (max - min > 1) {
		med = (min + max) >> 1;
		if (mpq_cmp(match[med]->sum, boundary) < 0)
			max = med;
		else
			min = med;
	}
	return max;
}

/*
	given an array of (match_t *) sorted by sum, find the index of the
	match with sum == boundary. Returns -1 if no match found.
*/
int pp_exact_chop(array_t *array, mpq_t boundary) {
	int min, max, med;
	match_t **match = (match_t **)array->array;
	min = -1;
	max = array->count;
	while (max - min > 1) {
		med = (min + max) >> 1;
		if (mpq_cmp(match[med]->sum, boundary) < 0)
			max = med;
		else
			min = med; 
	}
	return mpq_equal(match[max]->sum, boundary) ? max : -1;
}

match_t *pp_find(pp_t *self, ulong mod, mpq_t max) {
	ulong p, dismod, *cmodp, *recp, curv, *curp;
	int count, number, i, j, setcount, besti, havelfm, subseq, cmp;
	mpq_t min, *csump, psum;
	array_t *match, cur, csum, cmod, *cachearray;
	cache_t *cache;
	match_t *best, *lfm, *result;

#ifdef DEBUG
	gmp_fprintf(stderr, "pp_find(pp[%lu], mod=%lu, max=%Qd)\n", self->pp, mod, max);
#endif
	match = ((array_t *)self->match.array) + mod;
	if (match->count) {
		match_t **mp = (match_t **)match->array;
		if (mpq_cmp(mp[match->count - 1]->sum, max) < 0) {
			int index = pp_find_chop(match, max);
			if (self->index >= diagnose_retry) pp_diagnose(self, index);
			return mp[index];
		}
	}
	if (test_vector(self->complete, mod)) return (match_t *)NULL;

	QINIT(&min, "find pp=%lu k=%lu min", self->pp, k);
	mpq_sub(min, self->sum, max);
	QINIT(&psum, "find pp=%lu k=%lu psum", self->pp, k);

	cur.count = cur.space = 0;
	csum.count = csum.space = 0;
	cmod.count = cmod.space = 0;

	p = self->p;
	recp = (ulong *)self->rec.array;
	dismod = (self->modsum + p - mod) % p;
	number = self->count;

	best = new_match();
	mpq_set(best->sum, self->sum);
	mpqi_add_ui(best->sum, 1);
	besti = -1;

	cachearray = ((array_t *)self->cache.array) + mod;
	cache = (cache_t *)cachearray->array;
	for (count = 0; count < cachearray->count; ++count) {
		if (!cache[count].havelfm) continue;
		if (mpq_cmp(cache[count].lfm->sum, best->sum) >= 0) continue;
		free_match(best);
		best = copy_match(cache[count].lfm);
		besti = count;
	}
	array_resize(cachearray, sizeof(cache_t), number + 1);
	cache = (cache_t *)cachearray->array;
	for (i = cachearray->count; i <= number; ++i) {
		cache[i].nfm.count = cache[i].nfm.space = 0;
		cache[i].complete = cache[i].havelfm = 0;
	}
	cachearray->count = number + 1;
	for (count = 0; count <= number; ++count) {
/* RIGHTFIX? array underflow copied from perl code: */
/*   last if $self->{total}[$number - 1][$count - 1] > $best->[0]; */
		if (count && mpq_cmp(*get_total(self->totals, number - 1, count - 1), best->sum) > 0)
			break;
		if (cache[count].havelfm || cache[count].complete) continue;
		array_resize(&cur, sizeof(ulong), cache[count].nfm.count);
		memcpy(cur.array, cache[count].nfm.array, sizeof(ulong) * cache[count].nfm.count);
		cur.count = cache[count].nfm.count;
		cache[count].nfm.count = 0;
		if (csum.count < cur.count + 1) {
			array_resize(&csum, sizeof(mpq_t), cur.count + 1);
			csump = (mpq_t *)csum.array;
			for (i = csum.count; i < cur.count + 1; ++i)
				QINIT(&csump[i], "find pp=%lu k=%lu csum[%d]", self->pp, k, i);
			csum.count = cur.count + 1;
			array_resize(&cmod, sizeof(ulong), cur.count + 1);
			cmodp = (ulong *)cmod.array;
		}

		{
			mpq_t rtcsum;
			ulong rtcmod = 0;
			QINIT(&rtcsum, "find pp=%lu k=%lu rtcsum", self->pp, k);
			mpq_set(csump[0], rtcsum);
			cmodp[0] = rtcmod;
			for (i = 0; i < cur.count; ++i) {
				curv = curp[i];
				mpq_add(rtcsum, rtcsum, *get_total(self->totals, curv, 0));
				mpq_set(csump[i + 1], rtcsum);
				rtcmod = (rtcmod + recp[curv]) % p;
				cmodp[i + 1] = rtcmod;
			}
			QCLEAR(&rtcsum, "find pp=%lu k=%lu rtcsum", self->pp, k);
		}
		havelfm = 0;
		setcount = cur.count;
		subseq = 1;
		while (subseq) {
			if (setcount == count && cmodp[setcount] == dismod && mpq_cmp(csump[setcount], min) > 0) {
				if (cache[count].nfm.count == 0) {
					array_resize(&(cache[count].nfm), sizeof(ulong), setcount);
					memcpy(cache[count].nfm.array, cur.array, sizeof(ulong) * setcount);
					cache[count].nfm.count = cur.count;
				}
				if (!havelfm || mpq_cmp(lfm->sum, csump[setcount]) > 0) {
					if (havelfm) free_match(lfm);
					lfm = new_match();
					lfm->v = new_vector(self->count);
					fill_vector(lfm->v, 1);
					for (i = 0; i < setcount; ++i)
						set_vector(lfm->v, ((ulong *)cur.array)[i], 0);
					mpq_set(lfm->sum, csump[setcount]);
					havelfm = 1;
				}
			}
			if (setcount < count) {
				array_resize(&cur, sizeof(ulong), setcount + 1);
				((ulong *)cur.array)[setcount] = setcount ? ((ulong *)cur.array)[setcount - 1] : number;
				++setcount;
			}
			while (1) {
				if (setcount == 0) {
					if (cache[count].nfm.count == 0) cache[count].complete = 1;
					subseq = 0;
					break;	/* last SUBSEQ */
				}
				i = --(((ulong *)cur.array)[setcount - 1]);
				if (i < count - setcount) {
					--setcount;
					continue;	/* next INDEX */
				}
				mpq_set(psum, csump[setcount - 1]);
				if (havelfm) {
					mpq_add(q0, psum, *get_total(self->totals, i, count - setcount));
					if (mpq_cmp(q0, lfm->sum) > 0) {
						--setcount;
						continue;	/* next INDEX */
					}
				}
				if (csum.count < setcount + 1) {
					array_resize(&csum, sizeof(mpq_t), setcount + 1);
					csump = (mpq_t *)csum.array;
					array_resize(&cmod, sizeof(ulong), setcount + 1);
					cmodp = (ulong *)cmod.array;
					for (j = csum.count; j < setcount + 1; ++j)
						QINIT(&csump[j], "find pp=%lu k=%lu csum[%d++]", self->pp, k, j);
					csum.count = setcount + 1;
				}
				mpq_add(csump[setcount], psum, *get_total(self->totals, i, 0));
				cmodp[setcount] = (cmodp[setcount - 1] + recp[i]) % p;
				break;	/* next SUBSEQ */
			}
		}
		if (!havelfm) continue;
		if (cache[count].havelfm) free_match(cache[count].lfm);
		cache[count].havelfm = 1;
		cache[count].lfm = lfm;
		cmp = mpq_cmp(lfm->sum, best->sum);
		if (cmp > 0 || (cmp == 0 && count > besti)) continue;
		free_match(best);
		best = copy_match(lfm);
		besti = count;
	}
	if (besti >= 0) {
		array_resize(match, sizeof(match_t *), match->count + 1);
		result = ((match_t **)match->array)[match->count++] = new_match();
		if (cache[besti].havelfm) {
			free_match(cache[besti].lfm);
			cache[besti].havelfm = 0;
		}
		for (i = besti + 1; i < cachearray->count; ++i) {
			if (cache[i].havelfm && mpq_equal(cache[i].lfm->sum, best->sum)) {
				free_match(cache[i].lfm);
				cache[i].havelfm = 0;
			}
		}
		mpq_sub(result->sum, self->sum, best->sum);
		if (best->v) result->v = copy_vector(best->v);
		if (self->index >= diagnose_retry) pp_diagnose(self, match->count - 1);
	} else {
		set_vector(self->complete, mod, 1);
		result = (match_t *)NULL;
	}

	QCLEAR(&min, "find pp=%lu k=%lu min", self->pp, k);
	QCLEAR(&psum, "find pp=%lu k=%lu psum", self->pp, k);
	csump = (mpq_t *)csum.array;
	for (i = 0; i < csum.count; ++i)
		QCLEAR(&csump[i], "find pp=%lu k=%lu csum[%d]", self->pp, k, i);
	if (csum.space) free(csum.array);
	if (cmod.space) free(cmod.array);
	if (cur.space) free(cur.array);
	free_match(best);

	return result;
}

match_t *pp_find_exact(pp_t *self, ulong mod, mpq_t discard) {
	match_t *result;
	int count, number, i, setcount, subseq;
	ulong p, dismod, curv, *curp, *cmodp, *recp;
	mpq_t sum, keep, *csump;
	array_t *match, *cachearray, cur, csum, cmod;
	cache_t *cache;

#ifdef DEBUG
	gmp_fprintf(stderr, "pp_find_exact(pp[%lu], mod=%lu, discard=%Qd)\n", self->pp, mod, discard);
#endif
	QINIT(&sum, "find_exact pp=%lu k=%lu sum", self->pp, k);
	mpq_set(sum, self->sum);
	QINIT(&keep, "find exact pp=%lu k=%lu kept", self->pp, k);
	mpq_sub(keep, sum, discard);

	recp = (ulong *)self->rec.array;
	cur.count = cur.space = 0;
	csum.count = csum.space = 0;
	cmod.count = cmod.space = 0;

	match = ((array_t *)self->match.array) + mod;
	if (match->count) {
		match_t **mp = (match_t **)match->array;
		if (mpq_cmp(mp[match->count - 1]->sum, keep) <= 0) {
			int index = pp_exact_chop(match, keep);
			/* FIXME: should this copy_match()? */
			result = (index < 0) ? (match_t *)NULL : mp[index];
			goto EXACT_DONE;
		}
	}
	if (test_vector(self->complete, mod)) {
		result = (match_t *)NULL;
		goto EXACT_DONE;
	}
	p = self->p;
	dismod = (self->modsum + p - mod) % p;
	number = self->count;
	cachearray = ((array_t *)self->cache.array) + mod;
	array_resize(cachearray, sizeof(cache_t), number + 1);
	cache = (cache_t *)cachearray->array;
	for (i = cachearray->count; i <= number; ++i) {
		cache[i].nfm.count = cache[i].nfm.space = 0;
		cache[i].complete = cache[i].havelfm = 0;
	}
	cachearray->count = number + 1;
	for (count = 0; count <= number; ++count) {
		/* RIGHTFIX? as above */
		if (count && mpq_cmp(*get_total(self->totals, number - 1, count - 1), discard) > 0)
			break;
		if (cache[count].havelfm && mpq_equal(cache[count].lfm->sum, discard)) {
			result = new_match();
			mpq_sub(result->sum, self->sum, cache[count].lfm->sum);
			result->v = copy_vector(cache[count].lfm->v);
			goto EXACT_DONE;
		}
		if (cache[count].complete) continue;

		array_resize(&cur, sizeof(ulong), cache[count].nfm.count);
		memcpy(cur.array, cache[count].nfm.array, sizeof(ulong) * cache[count].nfm.count);
		cur.count = cache[count].nfm.count;
		cache[count].nfm.count = 0;

		array_resize(&csum, sizeof(mpq_t), cur.count + 1);
		csump = (mpq_t *)csum.array;
		for (i = 0; i < cur.count + 1; ++i)
			QINIT(&csump[i], "find_exact pp=%lu k=%lu csum[%d]", self->pp, k, i);
		csum.count = cur.count + 1;
		array_resize(&cmod, sizeof(ulong), cur.count + 1);
		cmodp = (ulong *)cmod.array;

		{
			mpq_t rtcsum;
			ulong rtcmod = 0;
			QINIT(&rtcsum, "find_exact pp=%lu k=%lu rtcsum", self->pp, k);
			mpq_set(csump[0], rtcsum);
			cmodp[0] = rtcmod;
			for (i = 0; i < cur.count; ++i) {
				curv = curp[i];
				mpq_add(rtcsum, rtcsum, *get_total(self->totals, curv, 0));
				mpq_set(csump[i + 1], rtcsum);
				rtcmod = (rtcmod + recp[curv]) % p;
				cmodp[i + 1] = rtcmod;
			}
			QCLEAR(&rtcsum, "find_exact pp=%lu k=%lu rtcsum", self->pp, k);
		}
		setcount = cur.count;
		subseq = 1;
		while (subseq) {
			if (setcount == count && cmodp[setcount] == dismod && mpq_equal(csump[setcount], discard)) {
				result = new_match();
				mpq_sub(result->sum, self->sum, csump[setcount]);
				result->v = new_vector(number);
				fill_vector(result->v, 1);
				for (i = 0; i < cur.count; ++i)
					set_vector(result->v, curp[i], 0);
				goto EXACT_DONE;
			}
			if (setcount < count) {
				array_resize(&cur, sizeof(ulong), setcount + 1);
				curp[setcount] = cur.count ? curp[setcount - 1] : number;
				++setcount;
			}
			while (1) {
				if (setcount == 0) {
					subseq = 0;
					break;
				}
				i = --curp[setcount - 1];
				if (i < count - setcount) {
					--setcount;
					continue;
				}
				mpq_add(q0, csump[setcount - 1], *get_total(self->totals, i, count - setcount));
				if (mpq_cmp(q0, discard) > 0) {
					--setcount;
					continue;
				}
				array_resize(&csum, sizeof(mpq_t), setcount + 1);
				csump = (mpq_t *)csum.array;
				while (setcount + 1 < csum.count)
					QINIT(&csump[csum.count++], "find_exact pp=%lu k=%lu csum[%lu++]", self->pp, k, csum.count);
				array_resize(&cmod, sizeof(ulong), setcount + 1);
				cmodp = (ulong *)cmod.array;
				mpq_add(csump[setcount], csump[setcount - 1], *get_total(self->totals, i, 0));
				cmodp[setcount] = (cmodp[setcount - 1] + recp[i]) % p;
				break;
			}
		}
	}
	result = (match_t *)NULL;
  EXACT_DONE:
	QCLEAR(&sum, "find_exact pp=%lu k=%lu sum", self->pp, k);
	QCLEAR(&keep, "find exact pp=%lu k=%lu kept", self->pp, k);
	csump = (mpq_t *)csum.array;
	for (i = 0; i < csum.count; ++i)
		QCLEAR(&csump[i], "find exact pp=%lu k=%lu csum[%d]", self->pp, k, i);
	if (csum.space) free(csum.array);
	if (cmod.space) free(cmod.array);
	if (cur.space) free(cur.array);
	return result;
}

pp_t* pp_account(ulong k) {
	pp_t* self;
	ulong pp, p, inv;
	vector_t *pspan, **spanp;
	int count, i, j, cansolve, complete_0;
	array_t *match, *cap;
	match_t **mp;

#ifdef DEBUG
	fprintf(stderr, "pp_account(k=%lu)\n", k);
#endif
	/* H_k += 1/k */
	mpq_set_ui(q0, 1, k);
	mpq_add(H_k, H_k, q0);

	pp = gpp(k);
	self = (pp == k) ? pp_new(k) : pp_locate(pp);
	count = ++self->count;
	p = self->p;

	/* push k on val[] */
	array_resize(&self->val, sizeof(ulong), count);
	((ulong *)self->val.array)[count - 1] = k;

	inv = mod_invert(k / pp, p);
	self->modsum = (self->modsum + inv) % p;

	/* push inv on rec */
	array_resize(&self->rec, sizeof(ulong), count);
	((ulong *)self->rec.array)[count - 1] = inv;

	/* self->sum += 1/k */
	mpq_add(self->sum, self->sum, q0);
	expand_total(self->totals, q0);

	/* push merged span of reachable moduli */
	array_resize(&self->span, sizeof(vector_t *), count + 1);
	spanp = (vector_t **)self->span.array;
	pspan = spanp[count - 1];
	spanp[count] = merge_vector(pspan, inv);
	++self->span.count;

	cansolve = self->cansolve = test_vector(pspan, (p - inv) % p);
	complete_0 = test_vector(self->complete, 0);
	fill_vector(self->complete, 0);
	if (!cansolve && complete_0) set_vector(self->complete, 0, 0);

	cap = (array_t *)self->cache.array;
	for (i = 0; i < p; ++i) free_cache_array(cap + i);

	match = (array_t *)self->match.array;
	for (i = cansolve ? 0 : 1; i < p; ++i) {
		mp = (match_t **)match[i].array;
		for (j = 0; j < match[i].count; ++j) free_match(mp[j]);
		match[i].count = 0;
	}

	if (cansolve) return self;
	for (i = 1; i + self->index < known_pp.count; ++i) {
		/* FIXME: ?? not needed if higher pp cannot affect me, so:
			if gpp > sqrt(k) and p.gpp >= k return null
			or if !grep {
				(pp[i].p == p && pp[i].cansolve)
				|| (pp[i].pp * gpp < k && pp[i].cansolve)
			} all i > index
		*/
		if (self[i].cansolve) return self;
	}
	return (pp_t *)NULL;
}

int pp_consider(pp_t *self) {
	/* ulong reqpp = self->pp; */
	/* int reqind = self->count; */
	int known = known_pp.count, i, mod, result;
	int try = known - 1;
	ulong botseen = k;
	mpq_t *sum, *kept, spare, tkept, max, discard;
	match_t **tried;
	pp_t *pp;

#ifdef DEBUG
	fprintf(stderr, "pp_consider(pp[%lu]) at k=%lu\n", self->pp, k);
#endif
	QINIT(&spare, "consider k=%lu spare", k);
	QINIT(&tkept, "consider k=%lu tkept", k);
	QINIT(&max, "consider k=%lu max", k);
	QINIT(&discard, "consider k=%lu discard", k);
	mpq_set_ui(rlimit, 1, k);
	array_resize(&consider_sum, sizeof(mpq_t), known + 1);
	array_resize(&consider_kept, sizeof(mpq_t), known + 1);
	sum = (mpq_t *)consider_sum.array;
	kept = (mpq_t *)consider_kept.array;
	for (i = consider_sum.count; i < known + 1; ++i) {
		QINIT(&sum[i], "consider k=%lu sum[%d]", k, i);
		QINIT(&kept[i], "consider k=%lu kept[%d]", k, i);
	}
	consider_sum.count = consider_kept.count = known + 1;
	array_resize(&consider_tried, sizeof(match_t *), known);
	tried = (match_t **)consider_tried.array;
	tried[try] = (match_t *)NULL;
	mpq_set_ui(q0, n, 1);
	mpq_sub(sum[try + 1], H_k, q0);
	mpq_set_ui(kept[try + 1], 0, 1);

	while (try < known) {
#ifdef DEBUG
		fprintf(stderr, "PP_consider(pp[%lu]): while (try=%u < known=%u)\n", self->pp, try, known);
#endif
		mpq_set(spare, sum[try + 1]);
		mpq_set(tkept, kept[try + 1]);
		pp = &(((pp_t *)known_pp.array)[try]);
		if (tried[try]) {
			mpq_set(max, tried[try]->sum);
		} else {
			mpq_set(max, pp->sum);
			mpqi_add_ui(max, 1);
		}
		if (pp->pp * pp->p <= k && mpz_fdiv_ui(mpq_denref(tkept), pp->pp) == 0) {
			mpz_set_ui(z0, pp->p);
			mpz_divexact_ui(z1, mpq_denref(tkept), pp->pp);
			mpz_invert(z1, z1, z0);
			mpz_mod_ui(z2, mpq_numref(tkept), pp->p);
			mpz_neg(z2, z2);
			mpz_mul(z1, z1, z2);
			mod = mpz_fdiv_ui(z1, pp->p);
		} else {
			mod = 0;
		}
		if (!(tried[try] = pp_find(pp, mod, max))) {
			++try;
			continue;
		}
		mpq_sub(discard, pp->sum, tried[try]->sum);
		mpq_add(kept[try], tkept, tried[try]->sum);
		mpq_sub(spare, spare, discard);
		mpq_set(sum[try], spare);
		if (mpq_sgn(spare) < 0) {
			++try;
			continue;
		}
		/* if (pp->pp == reqpp && !test_vector(tried[try]->v, reqind)) continue; */
		if (mpq_sgn(spare) > 0 && mpq_cmp(spare, rlimit) < 0) {
			if (tried[try] == pp_find_exact(pp, mod, spare)) {
				mpq_set_ui(spare, 0, 1);
			} else {
				++try;
				continue;
			}
		}
		if (mpq_sgn(spare) == 0) {
			result = pp_solution(tried, try);
			goto DONE;
		}
		if (mpz_cmp_ui(mpq_numref(spare), 1) == 0
			&& mpz_cmp_ui(mpq_denref(spare), k) <= 0
		) {
			gmp_printf("  maybe solution: %Qd\n", spare);
			if (pp_maybe_solution(tried, try, mpz_get_ui(mpq_denref(spare)))) {
				result = 1;
				goto DONE;
			}
		}
		if (try) tried[--try] = (match_t *)NULL;
		if (botseen >= pp->pp) botseen = ((pp_t *)known_pp.array)[try].pp;
	}
	pp_undiagnose(diagnose_len1 + diagnose_len2);
	printf("... reached %lu\n", botseen);
	result = 0;
  DONE:
	QCLEAR(&spare, "consider k=%lu spare", k);
	QCLEAR(&tkept, "consider k=%lu tkept", k);
	QCLEAR(&max, "consider k=%lu max", k);
	QCLEAR(&discard, "consider k=%lu discard", k);
	return result;
}

/*
	Try to find a sequence K of integers 1 <= i <= k such that sum 1/i = n
	and with firstk <= k <= lastk (or continuing until a sequence is found
	if lastk = 0).
	If no sequence is found for a given k, this guarantees A101877(n) > k.
	If a sequence is found for a given k, this guarantees A101877(n) <= k.
	On finding a sequence, it is printed and the program exits.
*/
void A101877(ulong firstk, ulong lastk) {
	pp_t* pp;

#ifdef DEBUG
	fprintf(stderr, "A101877(firstk=%lu, lastk=%lu) for n=%lu\n", firstk, lastk, n);
#endif
	for (k = 1; !lastk || k <= lastk; ++k) {
		pp = pp_account(k);
		if (k < firstk) continue;
		if (mpq_cmp_ui(H_k, n, 1) < 0) continue;
		if (pp) {
			printf("%lu: S = %.6f (t=%.2f)\n", k, mpq_get_d(H_k), timing());
			diagnose_retry = -1;
			diagnose_len1 = diagnose_len2 = 0;
			if (pp_consider(pp)) return;
		}
	}
}

void init(void) {
	/* init timing */
	clock_tick = sysconf(_SC_CLK_TCK);

	/* init primes */
	prime.count = 2;
	prime.space = 0;
	factor.count = 0;
	factor.space = 0;
	array_resize(&prime, sizeof(ulong), 2);
	((ulong *)prime.array)[0] = 2;
	((ulong *)prime.array)[1] = 3;
	prime.count = 2;
	array_resize(&factor, sizeof(ulong), 4);
	((ulong *)factor.array)[1] = 1;
	((ulong *)factor.array)[2] = 2;
	((ulong *)factor.array)[3] = 3;

	/* init fixed mp*_t */
	QINIT(&H_k, "init%d H_k", 0);
	QINIT(&q0, "init%d q0", 0);
	QINIT(&rlimit, "init%d rlimit", 0);

	mpz_init(z0);
	mpz_init(z1);
	mpz_init(z2);
	mpz_init(mi_r);
	mpz_init(mi_p);
	mpz_init(mi_k);

	/* init global arrays */
	known_pp.count = 0;
	known_pp.space = 0;

	solution = (vector_t *)NULL;
}

void clear(void) {
	int i;
	pp_t *ppp;
	mpq_t *qp;
	match_t **mp;

	/* clear fixed mp*_t */
	QCLEAR(&H_k, "clear%d H_k", 0);
	QCLEAR(&q0, "clear%d q0", 0);
	QCLEAR(&rlimit, "clear%d rlimit", 0);

	mpz_clear(z0);
	mpz_clear(z1);
	mpz_clear(z2);
	mpz_clear(mi_r);
	mpz_clear(mi_p);
	mpz_clear(mi_k);
	/* free lots of arrays and other malloced areas */
	if (solution) free_vector(solution);
	free(factor.array);
	free(prime.array);

	ppp = (pp_t *)known_pp.array;
	for (i = 0; i < known_pp.count; ++i) pp_free(ppp + i);
	free(known_pp.array);

	qp = (mpq_t *)consider_sum.array;
	for (i = 0; i < consider_sum.count; ++i)
		QCLEAR(&qp[i], "clear csum[%d]", i);
	free(consider_sum.array);
	qp = (mpq_t *)consider_kept.array;
	for (i = 0; i < consider_kept.count; ++i)
		QCLEAR(&qp[i], "clear kept[%d]", i);
	free(consider_kept.array);
	mp = (match_t **)consider_tried.array;
	for (i = 0; i < consider_tried.count; ++i) free_match(mp[i]);
	free(consider_tried.array);
}

int main(int argc, char** argv) {
	ulong firstk, lastk;
	if (argc < 2 || argc > 4) {
		fprintf(stderr, "Usage: %s <n> [ <first> [ <last> ]]\n", argv[0]);
		exit(-1);
	}
	n = atoi(argv[1]);
	if (n < 1) n = 4;
	firstk = 1;
	if (argc > 2) firstk = atoi(argv[2]);
	if (firstk < 1) firstk = 1;
	lastk = 0;
	if (argc > 3) lastk = atoi(argv[3]);
	if (lastk < 0) lastk = 0;
	if (setvbuf(stdout, (char*)NULL, _IONBF, (size_t)0)) {
		fprintf(stderr, "setvbuf failed\n");
		exit(-1);
	}
	init();
	A101877(firstk, lastk);
	clear();
	return 0;
}

/*
a(n): The least integer such that there is a sum of distinct unit fractions
equal to _n_, the greatest denominator being a(n).

If a(n) = k then there exist S = [ s_1, s_2, ... s_m ] such that s_m = k,
sum_1^m{1/s_i} = n, and 1 <= i < j <= m => s_i < s_j.

eg f(1) = 1: [ 1 ]
   f(2) = 6: [ 1, 2, 3, 6 ]
   f(3) = 24: [ 1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 18, 20, 24 ]

Proof of f(3): for k = 24, we have a candidate set of 1 .. 24, of which the
prime powers greater than 12 can immediately be discarded as unusable;
the multiples of 11 are unavailable since no partial sum of [ 2/1, 2/2 ]
is divisible by 11; the multiples of 7 are unavailable since no partial
sum of [ 6/1, 6/2, 6/3 ] is divisible by 7.

That leaves the candidate set as S U { 12 }, with a sum of 3 1/12. It
immediately follows that S is a candidate set with the right sum; further,
since the greatest denominator is 24, no two fractions from this candidate
set can be <= 1/12, so no candidate set with a lower maximal element can
sum to 3.

*********************************************************
***                                                   ***
*** see H6 for correct/up to date timings and results ***
***                                                   ***
*********************************************************

Multiple ways to get the same total for a given prime power are common.
I'm not sure if/when the code currently checks and skips for those cases,
or when it is safe to - I suspect this is the reqpp problem mentioned
below.
FIXME: we may skip a reqpp when a different set with the same sum would
pass the test - we'll never see the other set, so this may miss a solution.

Hmm, if when searching for LFM for some count we see an LFM already
recorded for another count, skip it and keep looking.

Hmm, when looping in find, we could notice *each* mod we pass, and
update other caches as we go. But that loses a lot of memory to store
matches relevant only 1/pp of the time.

Hmm, when processing $known[$cur], consider the minimum discard from
$known[$cur - 1] - we must either leave at least this much spare, or
choose our set to keep to be the value (mod $known[$cur - 1]->{f}) so
as to change the requirements for $cur - 1.

Hmm, since no pp can affect another when both are above sqrt(limit)
(at least if this is the highest represented power of p), we can run
through all pp > sqrt(limit) to find the least discard from each, and
reduce spare accordingly. This will be an effective prune when
backtracking reaches up into this range.

The approach used to test optimality of a(N) is:
- set Limit successively to 1, 2, ...
- calculate Q_{p^k} as the subsequence of 1..Limit for which gpp(n) = p^k
- set H(Limit) = sum_{i=1..Limit}{ 1/i }
- set Kept = 0/1, Discarded = 0/1
- recursively for each p^k such that Q_{p^k} is non-empty, in decreasing order
  - calculate Req = the required sum modulo p :
    - if denominator(Kept) is not divisible by p^k then Req = 0
      else Req = (- numerator(Kept) * ( denominator(Kept) / p^k ) ^ -1) (mod p)
  - for Q_{p^k}, calculate R, the sequence (q_i/p^k) ^ -1 (mod p)
  - save the values of Kept and Discarded
  - for each subsequence S of R such that sum(S) == Req (mod p)
    - for each q_i in Q, add 1/q_i to Kept or Discarded depending
      on whether the corresponding r_i was used ('kept') in S
    - if H(Limit) - Discarded = N then terminate successfully
    - if H(Limit) - Discarded > N then recurse to the next lower prime power
    - restore Kept and Discarded to the saved values
  - derecurse
- and try the next value for Limit

(All calculations are done using exact rationals; various optimisations
are elided.)

Note that this approach guarantees that denominator(H(Limit) - Discarded)
is free of any prime power >= p^k before we recurse to the next lower
prime power.

*/

