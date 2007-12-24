#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <string.h>

#define DEBUG_CMAX 0

#define MAXN 16

/* global */
mpz_t total, pmax;

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
	mpz_init(total);
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
	mpz_clear(total);
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

int tau(mpz_t n) {
	int k = 1, pi = 0, pow;
	unsigned long un;
	if (!mpz_fits_ulong_p(n)) {
		while (1) {
			if (pi == pz) {
				if (mpz_cmp(n, pmax) < 0) {
					pow = 2;
					goto tau_final;
				}
				gmp_fprintf(stderr, "tau(%Zdz) out of range\n", n);
				exit(1);
			}
			if (mpz_divisible_ui_p(n, prime[pi])) {
				pow = 2;
				mpz_divexact_ui(n, n, prime[pi]);
				while (mpz_divisible_ui_p(n, prime[pi])) {
					++pow;
					mpz_divexact_ui(n, n, prime[pi]);
				}
				k *= pow;
				if (mpz_fits_ulong_p(n)) {
					++pi;
					break;
				}
			}
			++pi;
		}
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		if (pi == pz) {
			/* PLIM >= 2^16 guarantees un is prime */
			pow = 2;
			goto tau_final;
		}
		if ((un % prime[pi]) == 0) {
			pow = 2;
			un /= prime[pi];
			while ((un % prime[pi]) == 0) {
				++pow;
				un /= prime[pi];
			}
			k *= pow;
		}
		if (prime[pi] * prime[pi] > un) {
			if (un > 1) k *= 2;
			return k;
		}
		++pi;
	}
	pi = 0;
	pow = 1;
	while (un > 1) {
		if (pi == small_factor[un]) {
			++pow;
		} else {
			k *= pow;
			pow = 2;
			pi = small_factor[un];
		}
		un /= pi;
	}
  tau_final:
	return k * pow;
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
/* fmax = (1 + #p2 = divisors(p2)) \ 2;
/* q2 = (-q2) % d;
/* total += sum(i = 1, fmax, q2 == p2[i] % d & 1);
 */
			} else {
				  /* total += ceil( || f: f | p.q || / 2 ) */
				mpz_set_ui(d, tau(p2));
				mpz_mul_ui(d, d, tau(q2));
				mpz_add_ui(d, d, 1);
				mpz_fdiv_q_2exp(d, d, 1);
				mpz_add(total, total, d);
			}
		} else {
/* fmax = (1 + #p2 = divisors(p2)) \ 2;
/* vmin = 1;
/* vmax = min(minf, fmax);
/* while (vmin < vmax,
/*   vmed = (vmin + vmax) \ 2;
/*   if (p2[vmed] < minf, vmin = vmed + 1, vmax = vmed);
/* );
/* if (d > 1, /* else d==1 since we know d>0 */
/*   q2 = (-q2) % d;
/*   sum(i = vmin, fmax, q2 == p2[i] % d & 1);
/* , 
/*   fmax + 1 - vmin;
/* );
 */
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

*/
