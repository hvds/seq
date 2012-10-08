#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>

typedef unsigned int uint;
typedef unsigned char uchar;
#define UINT_MAX 0xffffffff
#define VALIGN 4

typedef struct pp_s {
	/* fixed for life */
	uint p;			/* this prime */
	uint vecbase;	/* location of bit vector in vector stack */

	/* varying when k changes */
	uint max;		/* ceil(k / p) */
	uint stolen;	/* floor(k / p / p_0) for primes other than first */

	/* varying when offset changes */
	int offset;		/* try run starting at -offset (mod p) */
	uint flags;
/* set if all offsets up to this point are symmetric for the run */
#define F_SYMMETRIC 1
/* set if we've skipped one or more offsets setting a single value */
#define F_SINGLETON 2
/* set if we've set a single value this time through */
#define F_IS_SINGLE 4

	/* varying when k or offset changes */
	uint needed;	/* number of bits clear in the vector */
	int excess;		/* degrees of freedom */
} pp_t;

pp_t* pp;
uint pn;
uint k = 0;
uchar* v;
uint vsize;
double t0;
long CLK_TCK;

double timer(void) {
	struct tms tbuf;
	times(&tbuf);
	return (double)tbuf.tms_utime / CLK_TCK;
}

int JdSetPP(int index) {
	int excess = -k;
	int i, offset, vbyte, vbit;
	uint needed = k, symmetric = 1, found, match;
	uchar w[vsize];
	pp_t* pi;

	for (i = 0; i < pn; ++i) {
		pi = &pp[i];
		pi->max = (int)((k + pi->p - 1) / pi->p);
		pi->stolen = i > 0 ? (int)(k / pi->p / pp[0].p) : 0;
		excess += pi->max - pi->stolen;
	}

	memset(w, 0, vsize);
	for (i = 0; i < index; ++i) {
		pi = &pp[i];
		memcpy(v + pi->vecbase, w, vsize);
		pi->needed = needed;
		pi->excess = excess;
		pi->flags = (symmetric ? F_SYMMETRIC : 0)
				| (pi->flags & F_SINGLETON);	/* F_IS_SINGLE cleared */

		found = 0;
		if (pi->offset + pi->p >= k) {
			/* if max = 2, but at this offset we only hit 1, we can treat
			   remaining offsets the same as max=1 - act as if we hit a
			   single value, but don't set any value on the vector */
			found = 1;
			symmetric = 0;
			pi->offset = pi->p;
			pi->flags |= F_IS_SINGLE;
		} else {
			for (offset = pi->offset; offset < k; offset += pi->p) {
				vbyte = offset / 8;
				vbit = 1 << (offset & 7);
				if ((w[vbyte] & vbit) == 0) {
					++found;
					w[vbyte] |= vbit;
					match = offset;
				}
			}
			if (found == 1) {
				w[match >> 3] &= ~(1 << (match & 7));
				symmetric = 0;
				pi->flags |= F_IS_SINGLE;
			}
		}
		excess += pi->stolen - pi->max + found;
		needed -= found;

		if (symmetric) {
			int reverse_offset = (k - 1 - pi->offset) % pi->p;
			if (pi->offset != reverse_offset) {
				symmetric = 0;
				if (pi->offset > reverse_offset) {
					index = i;
				}
			}
		}
	}
	pi = &pp[index];
	memcpy(v + pi->vecbase, w, vsize);
	pi->needed = needed;
	pi->excess = excess;
	pi->flags = (symmetric ? F_SYMMETRIC : 0)
			| (pi->flags & F_SINGLETON);	/* F_IS_SINGLE cleared */
	return index;
}

__attribute__((noinline)) int JdRecordSolution(int index) {
	int i;
	double t1;

	/* we found a solution, record it */
	t1 = timer();
	printf("[%.2f] Jd = %u [", t1 - t0, k);
	for (i = 0; i < pn; ++i) {
		if (i) printf(".");
		if (i >= index || pp[i].flags & F_IS_SINGLE) {
			printf("*");
		} else {
			printf("%i", pp[i].offset);
		}
	}
	printf("]\n");

	/* now try the next length */
	++k;
	return JdSetPP(index);
}

uint Jd(void) {
	int index;

	if (k <= pn) {
		k = pn + 1;
	}
	pp[0].offset = -1;
	index = JdSetPP(0);

	while (index >= 0) {
		pp_t* pi = &pp[index];
		pp_t* pj = &pp[index + 1];
		uint found;
		int offset;
		uchar* w;
		int vbyte, vbit;

		/* break out to record solution for this k if we have at least
		   as many unused primes as unfilled slots */
		if (pn - index >= pi->needed) {
			index = JdRecordSolution(index);
			/* retry at the next length */
			continue;
		}

		/* if this (and hence all remaining primes) can only appear once,
		   we have no further chance to find a solution here */
		if (pi->max == 1) {
			--index;
			continue;
		}

		/* for this prime, we want to try starting points from 0 to p-1 */
		if (++pi->offset >= pi->p) {
			--index;
			continue;
		}

		pj->flags = 0;	/* reset F_SYMMETRIC, F_SINGLETON and F_IS_SINGLE */
		if (pi->flags & F_SYMMETRIC) {
			int reverse_offset = (k - 1 - pi->offset) % pi->p;
			if (pi->offset > reverse_offset) {
				continue;
			} else if (pi->offset == reverse_offset) {
				pj->flags = F_SYMMETRIC;
			}
		}

		w = v + pj->vecbase;

		found = 0;
		if (pi->offset + pi->p >= k) {
			/* if max = 2, but at this offset we only hit 1, we can treat
			   remaining offsets the same as max=1 - act as if we hit a
			   single value, but don't set any value on the vector */
			pi->offset = pi->p;	/* bump offset to abort on backtracking */
			found = 1;
			/* fall through to found == 1 case below */
		} else {
			memcpy(w, v + pi->vecbase, vsize);
			for (offset = pi->offset; offset < k; offset += pi->p) {
				vbyte = offset / 8;
				vbit = 1 << (offset & 7);
				if ((w[vbyte] & vbit) == 0) {
					++found;
					w[vbyte] |= vbit;
				}
			}
			if (found == 0) {
				/* that was useless, we can do better */
				continue;
			}
		}

		if (found == 1) {
			if (pi->flags & F_SINGLETON) {
				continue;
			}
			pi->flags |= F_SINGLETON | F_IS_SINGLE;
			memcpy(w, v + pi->vecbase, vsize);
		} else {
			pi->flags &= ~F_IS_SINGLE;
		}

		/* if we now have negative excess, no solution is possible */
		pj->excess = pi->excess + pi->stolen - pi->max + found;
		if (pj->excess < 0) {
			continue;
		}
		pj->needed = pi->needed - found;

		/* we've managed to assign this prime without contradiction,
		   so let's try the next */
		if (pj->p == 0) {
			/* there are no more primes to try - do we have a solution? */
			if (pj->needed == 0) {
				int newindex;
				/* we do: advance the index so we record it correctly */
				newindex = JdRecordSolution(index + 1);
				if (newindex > index) {
					/* retry the same offset at the new length */
					--pi->offset;
				} else {
					/* we're continuing at some lower index */
					index = newindex;
				}
			}
			continue;
		}

		/* try the next prime */
		pj->offset = -1;
		++index;
	}
	/* we failed at length k */
	return k - 1;
}

uint Jd_raw_sorted(uint n, uint* rawp) {
	uint has_2 = 0, p = 0, i;
	uint result;

	pn = 0;
	pp = (pp_t*)malloc((n + 1) * sizeof(pp_t));

	for (i = 0; i < n; ++i) {
		if (rawp[i] != p) {
			p = rawp[i];
			if (p == 2) {
				has_2 = 1;
			} else {
				pp[pn++].p = p;
			}
		}
	}
	pp[pn].p = 0;

	/* working on the assumption that Jd(p_1, p_2, ..., p_n) cannot
	   exceed n^2 */
	vsize = (pn * pn + 7) / 8;
	if (vsize & (VALIGN - 1)) {
		vsize &= ~(VALIGN - 1);
		vsize += VALIGN;
	}
	v = (uchar*)malloc((pn + 1) * vsize);

	for (i = 0; i <= pn; ++i) {
		pp[i].vecbase = i * vsize;
		pp[i].flags = 0;
	}

	result = Jd();

	free(v);
	free(pp);
	return has_2 ? result * 2 + 1 : result;
}

int cmpuint(const void *p1, const void *p2) {
	return (*(uint*)p1) - (*(uint*)p2);
}

int main(int argc, char** argv) {
	uint* rawp;
	uint i, result;
	double t1;

	if (argc > 1 && argv[1][0] == '-' && argv[1][1] == 'k') {
		k = atoi(argv[1] + 2);
		--argc;
		++argv;
	}
	if (argc < 2) {
		fprintf(stderr, "Usage: cj <p_1> <p_2> ...\n");
		return 0;
	}
	CLK_TCK = sysconf(_SC_CLK_TCK);
	rawp = (uint*)malloc(argc * sizeof(uint));
	for (i = 1; i < argc; ++i) {
		rawp[i - 1] = (uint)atoi(argv[i]);
	}
	qsort(rawp, argc - 1, sizeof(uint), cmpuint);

	t0 = timer();
	result = Jd_raw_sorted(argc - 1, rawp);
	t1 = timer();
	printf("[%.2f] Jd = %u\n", t1 - t0, result);

	free(rawp);
	return 0;
}

/*
3#: [0.00] Jd = 4 [0.*.*]
4#: [0.00] Jd = 6 [1.0.*.*]
5#: [0.00] Jd = 10 [0.2.1.*.*]
6#: [0.00] Jd = 12 [1.3.2.0.*.*]
7#: [0.00] Jd = 16 [0.0.4.2.1.*.*]
8#: [0.00] Jd = 19 [0.1.0.2.4.*.*.*]
9#: [0.00] Jd = 22 [0.3.0.5.4.2.1.*.*]
10#: [0.00] Jd = 28 [0.1.3.8.7.5.4.2.*.*]
11#: [0.00] Jd = 32 [0.1.0.8.4.5.10.2.*.*.*]
12#: [0.00] Jd = 36 [1.0.0.1.11.9.8.6.3.2.*.*]
13#: [0.04] Jd = 44 [0.0.1.4.2.14.13.11.*.7.*.*.*]
14#: [0.34] Jd = 49 [0.0.1.4.2.14.0.11.17.16.7.*.*.*]
15#: [0.57] Jd = 52 [0.0.1.4.2.14.0.11.17.16.7.*.*.*.*]
16#: [2.24] Jd = 58 [0.0.0.0.0.0.0.0.8.1.16.2.4.*.*.*]
17#: [38.40] Jd = 65 [0.2.1.2.1.4.6.3.5.10.19.20.16.11.*.*.*]
18#: [157.18] Jd = 75 [1.1.4.2.7.12.8.0.9.14.17.3.5.15.*.*.*.*]
19#: [480.49] Jd = 86 [0.2.1.2.5.4.2.3.16.25.28.20.10.11.23.14.19.*.*]
20#: [2150.31] Jd = 94 [0.1.1.3.10.2.17.13.9.3.7.32.40.5.*.20.28.*.*.4]
21#: [10922.37] Jd = 99 [1.0.5.7.1.6.2.17.9.10.11.1.1.24.3.39.32.*.*.8.*]
22#: [29334.62] Jd = 107 [0.1.1.3.10.2.2.9.7.5.15.38.34.35.20.44.13.28.*.*.4.17]
23#: [57453.89] Jd = 116 [0.2.1.2.5.6.0.11.1.24.28.20.10.4.41.14.49.*.*.16.25.*.26]
24#: [105569.37] Jd = 128 [0.1.3.0.1.13.5.12.8.20.23.29.25.18.19.50.28.7.32.34.4.2.*.*]
*/
