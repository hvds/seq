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
	uint p;	/* this prime */
	uint remain;	/* number of primes from here to end of list */
	uint vecbase;	/* location of bit vector in vector stack */

	/* varying when k changes */
	uint max;	/* ceil(k / p) */
	uint stolen;	/* floor(k / p / p_0) for primes other than first */

	/* varying when offset changes */
	int offset;		/* try run starting at -offset (mod p) */

	/* varying when k or offset changes */
	uint needed;	/* number of bits clear in the vector */
	int excess;		/* degrees of freedom */
} pp_t;

pp_t* pp;
uint pn;
uint k;
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
	uint needed = k;
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

		excess += pi->stolen - pi->max;
		for (offset = pi->offset; offset < k; offset += pi->p) {
			vbyte = offset / 8;
			vbit = 1 << (offset & 7);
			if ((w[vbyte] & vbit) == 0) {
				++excess;
				--needed;
				w[vbyte] |= vbit;
			}
		}
	}
	pi = &pp[index];
	memcpy(v + pi->vecbase, w, vsize);
	pi->needed = needed;
	pi->excess = excess;
	return index;
}

__attribute__((noinline)) int JdRecordSolution(int index) {
	int i;
	double t1;

	/* we found a solution, record it */
	t1 = timer();
	printf("[%.2f] Jd = %u [", t1 - t0, k);
	for (i = 0; i < index; ++i) {
		printf(i ? ", %u=%i" : "%u=%i", pp[i].p, pp[i].offset);
	}
	printf("] floating %u\n", pn - index);

	/* now try the next length */
	++k;
	return JdSetPP(index);
}

uint Jd(void) {
	int index;

	k = pn;
	pp[0].offset = -1;
	pp[0].needed = k;
	memset(v + pp[0].vecbase, 0, vsize);

	++k;
	index = JdSetPP(0);

	while (index >= 0) {
		pp_t* pi = &pp[index];
		pp_t* pj = &pp[index + 1];
		int offset;
		uchar* w;
		int vbyte, vbit;

		/* break out to record solution for this k if we have at least
		   as many unused primes as unfilled slots */
		if (pi->remain >= pi->needed) {
			index = JdRecordSolution(index);
			/* retry at the next length */
			continue;
		}

		/* for this prime, we want to try starting points from 0 to p-1 */
		if (++pi->offset >= pi->p) {
			--index;
			continue;
		}

		/* if this (and hence all remaining primes) can only appear once,
		   we have no further chance to find a solution here */
		if (pi->max == 1) {
			--index;
			continue;
		}

		w = v + pj->vecbase;
		memcpy(w, v + pi->vecbase, vsize);
		pj->excess = pi->excess + pi->stolen - pi->max;
		pj->needed = pi->needed;
		for (offset = pi->offset; offset < k; offset += pi->p) {
			vbyte = offset / 8;
			vbit = 1 << (offset & 7);
			if ((w[vbyte] & vbit) == 0) {
				++pj->excess;
				--pj->needed;
				w[vbyte] |= vbit;
			}
		}

		/* if we now have negative excess, no solution is possible */
		if (pj->excess < 0) {
			continue;
		}

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
		pp[i].remain = pn - i;
		pp[i].vecbase = i * vsize;
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
3#: [0.00] Jd = 4 [3=0] floating 2
4#: [0.00] Jd = 6 [3=1, 5=0] floating 2
5#: [0.00] Jd = 10 [3=0, 5=2, 7=1] floating 2
6#: [0.00] Jd = 12 [3=1, 5=3, 7=2, 11=0] floating 2
7#: [0.00] Jd = 16 [3=0, 5=0, 7=4, 11=2, 13=1] floating 2
8#: [0.00] Jd = 19 [3=0, 5=1, 7=0, 11=2, 13=4] floating 3
9#: [0.00] Jd = 22 [3=0, 5=3, 7=0, 11=5, 13=4, 17=2, 19=1] floating 2
10#: [0.00] Jd = 28 [3=0, 5=1, 7=3, 11=8, 13=7, 17=5, 19=4, 23=2] floating 2
11#: [0.02] Jd = 32 [3=0, 5=1, 7=0, 11=8, 13=4, 17=5, 19=10, 23=2] floating 3
12#: [0.02] Jd = 36 [3=1, 5=0, 7=0, 11=1, 13=11, 17=9, 19=8, 23=6, 29=3, 31=2] floating 2
13#: [1.04] Jd = 44 [3=0, 5=0, 7=1, 11=4, 13=2, 17=14, 19=13, 23=11, 29=16, 31=7] floating 3
14#: [6.70] Jd = 49 [3=0, 5=0, 7=1, 11=4, 13=2, 17=14, 19=0, 23=11, 29=17, 31=16, 37=7] floating 3
15#: [15.00] Jd = 52 [3=0, 5=0, 7=1, 11=4, 13=2, 17=14, 19=0, 23=11, 29=17, 31=16, 37=7] floating 4
16#: [54.10] Jd = 58 [3=0, 5=0, 7=0, 11=0, 13=0, 17=0, 19=0, 23=0, 29=8, 31=1, 37=16, 41=2, 43=4] floating 3
17#: [3827.42] Jd = 65 [3=0, 5=2, 7=1, 11=2, 13=1, 17=4, 19=6, 23=3, 29=5, 31=10, 37=19, 41=20, 43=16, 47=11] floating 3
18#: [4255.09] Jd = 75 [3=1, 5=1, 7=4, 11=2, 13=7, 17=12, 19=8, 23=0, 29=9, 31=14, 37=17, 41=3, 43=5, 47=15] floating 4
19#: [14351.82] Jd = 86 [3=0, 5=2, 7=1, 11=2, 13=5, 17=4, 19=2, 23=3, 29=16, 31=25, 37=28, 41=20, 43=10, 47=11, 53=23, 59=14, 61=19] floating 2
20#: [59021.680] Jd = 94 [3=0, 5=1, 7=1, 11=3, 13=10, 17=2, 19=17, 23=13, 29=9, 31=3, 37=7, 41=32, 43=40, 47=5, 53=15, 59=20, 61=28, 67=35, 71=37, 73=4] floating 0
*/
