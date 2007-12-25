#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <string.h>

#define DEBUG_CMAX 0
#define MAXN 16
#define MAX_ULONG 0xfffffffful

/* global */
mpz_t total, pmax;

/* storage for count_* */
mpz_t tz;

/* variables for calc_max() */
mpz_t chigh, clow, t1, t2;

/* non-recursed variables for A118085_r() */
mpz_t d, tempz, cur, cmax, g, p2, q2, minf;

/* stack of recursed variables for A118085_r() */
typedef struct A118085_r_varstack {
	mpz_t cur, cmax, p2, q2;
} A118085_r_varstack_t;
A118085_r_varstack_t va[MAXN];

#define PLIM 100000
unsigned long small_factor[PLIM];
unsigned long prime[PLIM];	/* overkill */
int pz = 0;

unsigned long *divisors = (unsigned long *)NULL;
int zdivisors = 0;
int ndivisors;

#define DEBUG_LOOP 1000000
int debug_counter = DEBUG_LOOP;
int debug_size = 0;

void debug_clear(void) {
	while (debug_size) {
		printf("\x08 \x08");
		--debug_size;
	}
}

void debug_print(int n, mpz_t p, mpz_t q, mpz_t cmax, mpz_t cur) {
	debug_clear();
	debug_size = gmp_printf("%d %Zd %Zd %Zd %Zd ", n, p, q, cmax, cur);
}

void initp(void) {
	int i, lastp = 2;
	memset((void *)small_factor, 0, sizeof(small_factor));
	small_factor[1] = 1;
	while (lastp * lastp < PLIM) {
		small_factor[lastp] = lastp;
		prime[pz++] = lastp;
		for (i = lastp * lastp; i < PLIM; i += lastp) {
			small_factor[i] || (small_factor[i] = lastp);
		}
		for (++lastp; small_factor[lastp]; ++lastp)
			;
	}
	while (lastp < PLIM) {
		small_factor[lastp] = prime[pz++] = lastp;
		for (++lastp; small_factor[lastp]; ++lastp)
			;
	}
	mpz_init(pmax);
	mpz_set_ui(pmax, PLIM + 1);
	mpz_mul(pmax, pmax, pmax);
}
	
void initz(void) {
	int i;
	mpz_init(total); mpz_init(tz);
	mpz_init(chigh); mpz_init(clow); mpz_init(t1); mpz_init(t2);
	mpz_init(d); mpz_init(tempz); mpz_init(cur); mpz_init(cmax);
	mpz_init(g); mpz_init(p2); mpz_init(q2); mpz_init(minf);
	for (i = 0; i < MAXN; ++i) {
		mpz_init(va[i].cur);
		mpz_init(va[i].cmax);
		mpz_init(va[i].p2);
		mpz_init(va[i].q2);
	}
}

void cleanz(void) {
	int i;
	mpz_clear(total); mpz_clear(tz);
	mpz_clear(chigh); mpz_clear(clow); mpz_clear(t1); mpz_clear(t2);
	mpz_clear(d); mpz_clear(tempz); mpz_clear(cur); mpz_clear(cmax);
	mpz_clear(g); mpz_clear(p2); mpz_clear(q2); mpz_clear(minf);
	for (i = 0; i < MAXN; ++i) {
		mpz_clear(va[i].cur);
		mpz_clear(va[i].cmax);
		mpz_clear(va[i].p2);
		mpz_clear(va[i].q2);
	}
}

/*
  find the next prime power prime[pi]^k dividing *n0, pi >= *pi0
  returns k, and sets *pi0 = pi, *n0 = *n0 / prime[pi]^k
  if n is a prime > PLIM, returns 1 and sets *pi0 = -1, n = prime
*/
int nextpp_i(unsigned long *n0, int *pi0) {
	int pi = *pi0, k;
	unsigned long n = *n0;
	while (pi < pz) {
		if ((n % prime[pi]) == 0) {
			k = 1;
			n /= prime[pi];
			while ((n % prime[pi]) == 0) {
				++k;
				n /= prime[pi];
			}
			*pi0 = pi;
			*n0 = n;
			return k;
		}
		++pi;
	}
	/* assert PLIM * PLIM > MAX_ULONG */
	*pi0 = -1;
	return 1;
}

/*
  find the next prime power prime[pi]^k dividing n, pi >= *pi0
  returns k, and sets *pi0 = pi, n = n / prime[pi]^k
  if n is a prime > PLIM, returns 1 and sets *pi0 = -1, n = prime
*/
int nextpp_z(mpz_t n, int *pi0) {
	int pi = *pi0, k;
	while (pi < pz) {
		if (mpz_divisible_ui_p(n, prime[pi])) {
			k = 1;
			mpz_divexact_ui(n, n, prime[pi]);
			while (mpz_divisible_ui_p(n, prime[pi])) {
				++k;
				mpz_divexact_ui(n, n, prime[pi]);
			}
			*pi0 = pi;
			return k;
		}
		++pi;
	}
	if (mpz_cmp(n, pmax) < 0) {
		*pi0 = -1;
		return 1;
	}
	gmp_fprintf(stderr, "nextpp_z(%Zd) out of range\n", n);
}

void resize_divisors(int size) {
	divisors = (unsigned long *)realloc((void *)divisors, size * sizeof(unsigned long));
	zdivisors = size;
}

void init_divisors(void) {
	if (!divisors) {
		resize_divisors(100);
	}
	divisors[0] = 1;
	ndivisors = 1;
}

/*
  combine current list of divisors with prime powers (p^0 .. p^k)
  ignores multiples that exceed maxf
*/
void append_divisors_pp(int k, unsigned long p, unsigned long maxf) {
	int i = 0, x, y;
	int maxd = ndivisors * (k + 1);
	unsigned long maxn = maxf / p;
	int divc = ndivisors;
	if (zdivisors < maxd) resize_divisors(maxd * 2);
	for (x = 0; x < k; ++x) {
		for (y = divc; y > 0; --y) {
			if (divisors[i] <= maxn) {
				divisors[ndivisors++] = divisors[i] * p;
			} else {
				--divc;
			}
			++i;
		}
	}
}

/*
  combine current list of divisors with all divisors of n
  ignores multiples that exceed maxf
*/
void append_divisors(mpz_t n, unsigned long maxf) {
	int k, pi = 0, lastp = 0;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		k = nextpp_z(n, &pi);
		if (pi < 0) return;	/* the remaining prime is already > MAX_ULONG */
		append_divisors_pp(k, prime[pi], maxf);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		k = nextpp_i(&un, &pi);
		if (pi < 0) {
			append_divisors_pp(k, un, maxf);
			return;
		}
		append_divisors_pp(k, prime[pi], maxf);
		++pi;
	}
	k = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++k;
		} else {
			if (k) append_divisors_pp(k, lastp, maxf);
			k = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (k) append_divisors_pp(k, lastp, maxf);
}

/*
  count the number of divisors of n
  modifies n
*/
int tau(mpz_t n) {
	int k = 1, pi = 0, lastp = 0, pow;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		k *= 1 + nextpp_z(n, &pi);
		if (pi < 0) return k;
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		k *= 1 + nextpp_i(&un, &pi);
		if (pi < 0) return k;
	}
	pow = 1;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++pow;
		} else {
			k *= pow;
			pow = 2;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	return k * pow;
}

/*
	add to total the number of divisors <= sqrt(pq)
*/
void count_all(mpz_t p, mpz_t q) {
	mpz_set_ui(tz, tau(p));
	mpz_mul_ui(tz, tz, tau(q));
	mpz_add_ui(tz, tz, 1);
	mpz_fdiv_q_2exp(tz, tz, 1);
	mpz_add(total, total, tz);
}

/*
	add to total the number of divisors <= sqrt(pq) that are == -q (mod d)
*/
void count_mod(mpz_t p, mpz_t q, mpz_t d) {
	int c = 0, i;
	unsigned long maxf;
	mpz_mul(tz, p, q);
	mpz_sqrt(tz, tz);
	maxf = mpz_fits_ulong_p(d) ? mpz_get_ui(d) : MAX_ULONG;
	init_divisors();
	append_divisors(p, maxf);
	mpz_set(tz, q);
	append_divisors(tz, maxf);
	for (i = 0; i < ndivisors; ++i) {
		mpz_add_ui(tz, q, divisors[i]);
		mpz_fdiv_r(tz, tz, d);
		if (mpz_cmp_ui(tz, 0) == 0) ++c;
	}
	if (c) mpz_add_ui(total, total, c);
}

/*
	add to total the number of divisors iminf <= d <= sqrt(pq) that are == -q (mod d)
*/
void count_min_mod(mpz_t p, mpz_t q, mpz_t d, unsigned long minf) {
	int c = 0, i;
	unsigned long maxf;
	mpz_mul(tz, p, q);
	mpz_sqrt(tz, tz);
	maxf = mpz_fits_ulong_p(tz) ? mpz_get_ui(tz) : MAX_ULONG;
	init_divisors();
	append_divisors(p, maxf);
	mpz_set(tz, q);
	append_divisors(tz, maxf);
	for (i = 0; i < ndivisors; ++i) {
		if (divisors[i] < minf) continue;
		mpz_add_ui(tz, q, divisors[i]);
		mpz_fdiv_r(tz, tz, d);
		if (mpz_cmp_ui(tz, 0) == 0) ++c;
	}
	if (c) mpz_add_ui(total, total, c);
}

/*
	add to total the number of divisors iminf <= d <= sqrt(pq)
*/
void count_min(mpz_t p, mpz_t q, mpz_t d, unsigned long minf) {
	int c = 0, i;
	unsigned long maxf;
	mpz_mul(tz, p, q);
	mpz_sqrt(tz, tz);
	maxf = mpz_fits_ulong_p(tz) ? mpz_get_ui(tz) : MAX_ULONG;
	init_divisors();
	append_divisors(p, maxf);
	append_divisors(q, maxf);
	for (i = 0; i < ndivisors; ++i) {
		if (divisors[i] >= minf) ++c;
	}
	if (c) mpz_add_ui(total, total, c);
}

void calc_max(mpz_t cmid, int n, mpz_t p, mpz_t q, mpz_t cmin) {
	mpz_set(cmid, cmin);
	int state = 0, cmp;
	while (1) {
		  /* t1 = p . c^n */
		mpz_pow_ui(t1, cmid, n);
		mpz_mul(t1, t1, p);
		  /* t2 = q . (c+1)^n */
		mpz_add_ui(t2, cmid, 1);
		mpz_pow_ui(t2, t2, n);
		mpz_mul(t2, t2, q);
		cmp = mpz_cmp(t1, t2);
#if DEBUG_CMAX
		gmp_fprintf(stderr, "cmax(n=%d, p=%Zd, q=%Zd, min=%Zd): try c=%Zd (state=%d): pc^n=%Zd, q(c+1)^n=%Zd, cmp=%d\n", n, p, q, cmin, cmid, state, t1, t2, cmp);
#endif
		switch (state) {
		  case 0:
			if (cmp > 0) {
				/* cmin is already too high, so return result < cmin */
				mpz_set_ui(cmid, 0);
				return;
			} else if (cmp == 0) {
				/* cmin is exactly right, return cmid == cmin */
				return;
			} else {
				/* cmin is low, mark it and enter doubling phase */
				mpz_set(clow, cmid);
			}
			state = 1;
			mpz_mul_2exp(cmid, cmid, 1);
			break;
		  case 1:
			if (cmp > 0) {
				/* finally found an upper bound, enter binary chop phase */
				mpz_set(chigh, cmid);
				mpz_add(cmid, clow, chigh);
				mpz_fdiv_q_2exp(cmid, cmid, 1);
				state = 2;
			} else if (cmp == 0) {
				/* cmid is exactly right, return it */
				return;
			} else {
				/* still too low, mark and continue */
				mpz_set(clow, cmid);
				mpz_mul_2exp(cmid, cmid, 1);
			}
			break;
		  case 2:
			if (cmp > 0) {
				mpz_set(chigh, cmid);
			} else if (cmp == 0) {
				return;
			} else {
				mpz_set(clow, cmid);
			}
			/* binary chop */
			mpz_sub(cmid, chigh, clow);
			if (mpz_cmp_ui(cmid, 1) <= 0) {
				/* we are converged, clow is the floor */
				mpz_set(cmid, clow);
				return;
			}
			mpz_fdiv_q_2exp(cmid, cmid, 1);
			mpz_add(cmid, cmid, clow);
			break;
		}
	}
}

void A118085_r(int n, mpz_t p, mpz_t q, mpz_t cmin) {
	int i;
	unsigned long iminf, imaxf;
	  /* p, q /= gcd(p, q) */
	mpz_gcd(g, p, q);
	if (mpz_cmp_ui(g, 1) > 0) {
		mpz_div(p, p, g);
		mpz_div(q, q, g);
	}
	  /* d = p - q */
	mpz_sub(d, p, q);
	  /* cur = max(cmin, q \ d + 1) */
	mpz_fdiv_q(tempz, q, d);
	mpz_add_ui(tempz, tempz, 1);
	mpz_set(cur, (mpz_cmp(cmin, tempz) < 0) ? tempz : cmin);

	  /* cmax = floor(1 / (sqrtn(p / q, n) - 1)) */
	  /* i.e. find max(c): p.c^n <= q.(c+1)^n    */
	calc_max(cmax, n, p, q, cur);
#if DEBUG_CMAX
	gmp_fprintf(stderr, "For n=%d, p=%Zd, q=%Zd, cmin=%Zd, got cmax=%Zd\n", n, p, q, cur, cmax);
#endif

	if (n > 3) {
		A118085_r_varstack_t *v = &va[n];
		mpz_set(v->cur, cur);
		mpz_set(v->cmax, cmax);
		for (; mpz_cmp(v->cur, v->cmax) <= 0; mpz_add_ui(v->cur, v->cur, 1)) {
			if (!--debug_counter) {
				debug_counter = DEBUG_LOOP;
				debug_print(n, p, q, cmax, cur);
			}
			mpz_mul(v->p2, p, v->cur);
			mpz_add_ui(v->q2, v->cur, 1);
			mpz_mul(v->q2, v->q2, q);
			A118085_r(n - 1, v->p2, v->q2, v->cur);
		}
		return;
	}
	for (; mpz_cmp(cur, cmax) <= 0; mpz_add_ui(cur, cur, 1)) {
		if (!--debug_counter) {
			debug_counter = DEBUG_LOOP;
			debug_print(n, p, q, cmax, cur);
		}
		mpz_mul(p2, p, cur);
		mpz_add_ui(q2, cur, 1);
		mpz_mul(q2, q2, q);
		mpz_gcd(g, p2, q2);
		if (mpz_cmp_ui(g, 1) > 0) {
			mpz_div(p2, p2, g);
			mpz_div(q2, q2, g);
		}
		mpz_sub(d, p2, q2);
		mpz_mul(minf, cur, d);
		mpz_sub(minf, minf, q2);
/* gmp_printf("%Zd/%Zd (%Zd) f>%Zd\n", p2, q2, d, minf); continue; */
		if (mpz_cmp_ui(minf, 1) <= 0) {
			if (mpz_cmp_ui(d, 1) > 0) {
				count_mod(p2, q2, d);
			} else {
				count_all(p2, q2);
			}
		} else {
			if (!mpz_fits_ulong_p(minf)) continue;
			iminf = mpz_get_ui(minf);
			if (mpz_cmp_ui(d, 1) > 0) {
				count_min_mod(p2, q2, d, iminf);
			} else {
				count_min(p2, q2, d, iminf);
			}
		}
	}
}

int main(int argc, char** argv) {
	int i, n;
	long clock_tick;
	mpz_t p, q, cur, tempz;
	clock_t t0;

	if (argc < 2) {
		fprintf(stderr, "Usage: gmseq <n> [ <a_1> <a_2> ... ]\n");
		exit(-1);
	}
	n = atoi(argv[1]);
	if (n < 1 || n > MAXN) {
		fprintf(stderr, "Error: 1 <= n <= %d required (got n = %d)\n", MAXN, n);
		exit(-1);
	}
	if (setvbuf(stdout, (char*)NULL, _IONBF, (size_t)0)) {
		fprintf(stderr, "setvbuf failed\n");
		exit(-1);
	}
	initp();
	initz();
	mpz_init_set_ui(p, 2);
	mpz_init_set_ui(q, 1);
	mpz_init_set_ui(cur, 0);
	mpz_init(tempz);
	for (i = 2; i < argc; ++i) {
		if (mpz_set_str(tempz, argv[i], 10) != 0) {
			fprintf(stderr, "Error parsing argument '%s'\n", argv[i]);
			exit(-1);
		}
		if (mpz_cmp(cur, tempz) > 0) {
			fprintf(stderr, "Argument '%s' is smaller than a previous value\n", argv[i]);
			exit(-1);
		}
		mpz_set(cur, tempz);
		mpz_mul(p, p, cur);
		mpz_add_ui(tempz, tempz, 1);
		mpz_mul(q, q, tempz);
		--n;
	}
	gmp_fprintf(stderr, "count A118085(n=%d, p=%Zd, q=%Zd, cmin=%Zd)\n",
			n, p, q, cur);
	if (n < 3) {
		fprintf(stderr, "This code requires at least 3 elements unspecified\n");
		exit(-1);
	}

	t0 = times((struct tms *)NULL);
	clock_tick = sysconf(_SC_CLK_TCK);
	mpz_set_ui(total, 0);
	A118085_r(n, p, q, cur);
	debug_clear();
	gmp_printf("result: %Zd (%.2f)\n", total, ((double)(times((struct tms *)NULL) - t0)) / clock_tick);
	cleanz();
	return 0;
}

/*
  A118085(n) counts the number of distinct multisets { a_1, a_2, ... a_n }
  such that prod{1+1/a_i} = 2.

  This program takes mandatory first parameter n, and optional additional
  parameters a_1, a_2, ... and counts the number of distinct n-element
  multisets satisfying the above property, with any specified a_i fixed,
  and with all new a_i greater than or equal to the largest specified a_i.

count A118085(n=7, p=2, q=1, cmin=0)
result: 6910688 (2.38)  

*/
