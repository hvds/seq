#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <string.h>

#define DEBUG_CMAX 0
#define MAXN 16
#define MAX_ULONG 0xfffffffful

int opt_v = 0;	/* print info on matching multisets if set */

/* global */
mpz_t total, pmax;

/* storage for append_* */
mpz_t maxn;

/* storage for count_* */
mpz_t tz, cm_maxf;

/* variables for calc_max() */
mpz_t chigh, clow, t1, t2;

/* non-recursed variables for A118085_r() */
mpz_t d, tempz, cur, cmax, g, p2, q2, minf;

/* stack of recursed variables for A118085_r() */
typedef struct A118085_r_varstack {
	mpz_t cur, cmax, p2, q2;
} A118085_r_varstack_t;
A118085_r_varstack_t va[MAXN];
int van;

/* for n=8, largest prime factor < 3.2^16 */
#define PLIM 200000
unsigned long small_factor[PLIM];
unsigned long prime[PLIM];	/* overkill */
int pz = 0;

mpz_t *zdivisors = (mpz_t *)NULL;
int zzdivisors = 0;
int ndivisors;
unsigned long *udivisors = (unsigned long *)NULL;
int zudivisors = 0;

#define SMALLMOD 32
unsigned long modlist_0[SMALLMOD];
unsigned long modlist_1[SMALLMOD];
unsigned long *modlist = modlist_0;
unsigned long *newmodlist = modlist_1;
unsigned long modulus;

#define DEBUG_LOOP 100000
int debug_counter = DEBUG_LOOP;
int debug_size = 0;

void debug_clear(void) {
	while (debug_size) {
		printf("\x08 \x08");
		--debug_size;
	}
}

int print_set(int n) {
	int i = 0, size = 0;
	for (i = van; i > n; --i) {
		if (i < van) size += printf(" ");
		size += gmp_printf("%Zd", va[i - 1].cur);
	}
	return size;
}

void debug_print(int n, mpz_t p, mpz_t q, mpz_t cmax, mpz_t cur) {
	debug_clear();
	debug_size += print_set(n);
	debug_size += gmp_printf(" %Zd: %d %Zd %Zd %Zd ", cur, n, p, q, cmax);
}

void initp(void) {
	int i, lastp = 2;
	memset((void *)small_factor, 0, sizeof(small_factor));
	small_factor[1] = 1;
	while (lastp * lastp < PLIM) {
		small_factor[lastp] = lastp;
		prime[pz++] = lastp;
		for (i = lastp * lastp; i < PLIM; i += lastp) {
			if (!small_factor[i]) small_factor[i] = lastp;
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
	mpz_init(total); mpz_init(maxn); mpz_init(tz); mpz_init(cm_maxf);
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
	mpz_clear(total); mpz_clear(maxn); mpz_clear(tz); mpz_clear(cm_maxf);
	mpz_clear(chigh); mpz_clear(clow); mpz_clear(t1); mpz_clear(t2);
	mpz_clear(d); mpz_clear(tempz); mpz_clear(cur); mpz_clear(cmax);
	mpz_clear(g); mpz_clear(p2); mpz_clear(q2); mpz_clear(minf);
	for (i = 0; i < MAXN; ++i) {
		mpz_clear(va[i].cur);
		mpz_clear(va[i].cmax);
		mpz_clear(va[i].p2);
		mpz_clear(va[i].q2);
	}
	for (i = 0; i < zzdivisors; ++i) mpz_clear(zdivisors[i]);
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
	exit(-1);
}

void resize_zdivisors(int size) {
	zdivisors = (mpz_t *)realloc((void *)zdivisors, size * sizeof(mpz_t));
	while (zzdivisors < size) mpz_init(zdivisors[zzdivisors++]);
}

void resize_udivisors(int size) {
	udivisors = (unsigned long *)realloc((void *)udivisors, size * sizeof(unsigned long));
	zudivisors = size;
}

void init_zdivisors(void) {
	if (!zdivisors) resize_zdivisors(100);
	mpz_set_ui(zdivisors[0], 1);
	ndivisors = 1;
}

void init_udivisors(void) {
	if (!udivisors) resize_udivisors(100);
	udivisors[0] = 1;
	ndivisors = 1;
}

/*
  combine current list of divisors with prime powers (p^0 .. p^k)
  ignores multiples that exceed maxf
*/
void append_zdivisors_pp(int k, unsigned long p, mpz_t maxf) {
	int i = 0, x, y;
	int maxd = ndivisors * (k + 1);
	int divc = ndivisors;
	mpz_fdiv_q_ui(maxn, maxf, p);
	if (zzdivisors < maxd) resize_zdivisors(maxd * 2);
	for (x = 0; x < k; ++x) {
		for (y = divc; y > 0; --y) {
			if (mpz_cmp(maxn, zdivisors[i]) >= 0) {
				mpz_mul_ui(zdivisors[ndivisors++], zdivisors[i], p);
			} else {
				--divc;
			}
			++i;
		}
	}
}

void append_udivisors_pp(int k, unsigned long p, unsigned long maxf) {
	int x, y;
	int maxd = ndivisors * (k + 1);
	int divc = ndivisors;
	unsigned long umaxn = maxf / p;
	unsigned long *dps, *dpd, f;
	if (zudivisors < maxd) resize_udivisors(maxd * 2);
	dps = udivisors;
	dpd = &udivisors[ndivisors];
	for (x = 0; x < k; ++x) {
		for (y = divc; y > 0; --y) {
			if (umaxn >= (f = *dps++)) {
				*dpd++ = f * p;
			} else {
				--divc;
			}
		}
	}
	ndivisors = dpd - udivisors;
}

/*
  combine current list of divisors with all divisors of n
  ignores multiples that exceed maxf
*/
void append_zdivisors(mpz_t n, mpz_t maxf) {
	int k, pi = 0, lastp = 0;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		k = nextpp_z(n, &pi);
		if (pi < 0) {
			gmp_fprintf(stderr, "append: prime %Zd missed from divisors\n", n);
			exit(-1);
		}
		append_zdivisors_pp(k, prime[pi], maxf);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		k = nextpp_i(&un, &pi);
		if (pi < 0) {
			append_zdivisors_pp(k, un, maxf);
			return;
		}
		append_zdivisors_pp(k, prime[pi], maxf);
		++pi;
	}
	k = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++k;
		} else {
			if (k) append_zdivisors_pp(k, lastp, maxf);
			k = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (k) append_zdivisors_pp(k, lastp, maxf);
}

void append_udivisors(mpz_t n, unsigned long maxf) {
	int k, pi = 0, lastp = 0;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		k = nextpp_z(n, &pi);
		if (pi < 0) return;	/* the remaining prime is already > MAX_ULONG */
		append_udivisors_pp(k, prime[pi], maxf);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		k = nextpp_i(&un, &pi);
		if (pi < 0) {
			append_udivisors_pp(k, un, maxf);
			return;
		}
		append_udivisors_pp(k, prime[pi], maxf);
		++pi;
	}
	k = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++k;
		} else {
			if (k) append_udivisors_pp(k, lastp, maxf);
			k = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (k) append_udivisors_pp(k, lastp, maxf);
}

void init_udivisors_mod(unsigned long mod) {
	modulus = mod;
	memset((void*)modlist, 0, mod * sizeof(unsigned long));
	modlist[1] = 1;
}

void append_udivisors_pp_mod(int k, unsigned long p) {
	int i, j, count;
	unsigned long vmod, *temp;
/*
	// This check is worth doing if gcd() is fast enough, and given
	// modulus <= SMALLMOD, it should be (we could even use a lookup table)
	if (p <= modulus && gcd(p, modulus) > 1) return;
*/
	p = p % modulus;
	memcpy((void*)newmodlist, (const void *)modlist, modulus * sizeof(unsigned long));
	for (i = 0; i < modulus; ++i) {
		if (!(count = modlist[i])) continue;
		vmod = i;
		for (j = 1; j <= k; ++j) {
			vmod = (vmod * p) % modulus;
			newmodlist[vmod] += count;
		}
	}
	temp = modlist; modlist = newmodlist; newmodlist = temp;
}

void append_udivisors_mod(mpz_t n) {
	int k, pi = 0, lastp = 0;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		k = nextpp_z(n, &pi);
		if (pi < 0) {
			un = mpz_fdiv_ui(n, modulus);
			append_udivisors_pp_mod(k, un);
			return;
		}
		append_udivisors_pp_mod(k, prime[pi]);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		k = nextpp_i(&un, &pi);
		if (pi < 0) {
			append_udivisors_pp_mod(k, un);
			return;
		}
		append_udivisors_pp_mod(k, prime[pi]);
		++pi;
	}
	k = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++k;
		} else {
			if (k) append_udivisors_pp_mod(k, lastp);
			k = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (k) append_udivisors_pp_mod(k, lastp);
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
  count the number of divisors of n, and accumulate actual divisors <= maxf
  modifies n
*/
int tau_zmax(mpz_t n, mpz_t maxf) {
	int k = 1, pi = 0, lastp = 0, pow;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		pow = nextpp_z(n, &pi);
		k *= pow + 1;
		if (pi < 0) {
			if (mpz_cmp(n, maxf) > 0) return k;
			gmp_fprintf(stderr, "tau_zmax: lost prime factor %Zd\n", n);
			exit(-1);
		}
		append_zdivisors_pp(pow, prime[pi], maxf);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		pow = nextpp_i(&un, &pi);
		k *= pow + 1;
		append_zdivisors_pp(pow, pi < 0 ? un : prime[pi], maxf);
		if (pi < 0) return k;
	}
	pow = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++pow;
		} else {
			if (pow) {
				k *= pow + 1;
				append_zdivisors_pp(pow, lastp, maxf);
			}
			pow = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (pow) append_zdivisors_pp(pow, lastp, maxf);
	return k * (pow + 1);
}

int tau_umax(mpz_t n, unsigned long maxf) {
	int k = 1, pi = 0, lastp = 0, pow;
	unsigned long un;
	while (!mpz_fits_ulong_p(n)) {
		pow = nextpp_z(n, &pi);
		k *= pow + 1;
		if (pi < 0) return k;	/* prime > MAX_ULONG */
		append_udivisors_pp(pow, prime[pi], maxf);
		++pi;
	}
	un = mpz_get_ui(n);
	while (un > PLIM) {
		pow = nextpp_i(&un, &pi);
		k *= pow + 1;
		append_udivisors_pp(pow, pi < 0 ? un : prime[pi], maxf);
		if (pi < 0) return k;
	}
	pow = 0;
	while (un > 1) {
		if (lastp == small_factor[un]) {
			++pow;
		} else {
			if (pow) {
				k *= pow + 1;
				append_udivisors_pp(pow, lastp, maxf);
			}
			pow = 1;
			lastp = small_factor[un];
		}
		un /= lastp;
	}
	if (pow) append_udivisors_pp(pow, lastp, maxf);
	return k * (pow + 1);
}

/*
	add to total the number of divisors <= sqrt(pq)
	== floor((tau(p) * tau(q) + 1) / 2)
*/
void count_all(mpz_t p, mpz_t q) {
	mpz_set_ui(tz, tau(p));
	mpz_mul_ui(tz, tz, tau(q));
	mpz_add_ui(tz, tz, 1);
	mpz_fdiv_q_2exp(tz, tz, 1);
	mpz_add(total, total, tz);
	if (opt_v) {
		debug_clear();
		print_set(2);
		gmp_printf(": %Zd (A)\n", tz);
	}
}

/*
	add to total the number of divisors minf <= d <= sqrt(pq)
	== count_all(p, q) - || d: d | pq, d < minf ||
*/
void count_min(mpz_t p, mpz_t q, mpz_t d, mpz_t minf) {
	if (mpz_fits_ulong_p(minf)) {
		unsigned long umaxf = mpz_get_ui(minf) - 1;
		init_udivisors();
		mpz_set_ui(tz, tau_umax(p, umaxf));
		mpz_mul_ui(tz, tz, tau_umax(q, umaxf));
	} else {
		init_zdivisors();
		mpz_sub_ui(minf, minf, 1);
		mpz_set_ui(tz, tau_zmax(p, minf));
		mpz_mul_ui(tz, tz, tau_zmax(q, minf));
	}
	mpz_add_ui(tz, tz, 1);
	mpz_fdiv_q_2exp(tz, tz, 1);
	mpz_sub_ui(tz, tz, ndivisors);
	mpz_add(total, total, tz);
	if (opt_v && mpz_cmp_ui(tz, 0) > 0) {
		debug_clear();
		print_set(2);
		gmp_printf(": %Zd (B >= %Zd+1)\n", tz, minf);
	}
}

/*
	add to total the number of divisors <= sqrt(pq) that are == -q (mod d)
*/
void count_mod(mpz_t p, mpz_t q, mpz_t d) {
	int c = 0, i;
	unsigned long umaxf, ud, qmod;
	mpz_mul(cm_maxf, p, q);
	mpz_sqrt(cm_maxf, cm_maxf);
	if (mpz_fits_ulong_p(cm_maxf)) {
		umaxf = mpz_get_ui(cm_maxf);
		if (mpz_fits_ulong_p(d)) {
			ud = mpz_get_ui(d);
			mpz_neg(tz, q);
			qmod = mpz_fdiv_ui(tz, ud);
			if (ud <= SMALLMOD) {
				init_udivisors_mod(ud);
				append_udivisors_mod(p);
				mpz_set(tz, q);
				append_udivisors_mod(tz);
				c = (modlist[qmod] + 1) >> 1;
			} else {
				init_udivisors();
				append_udivisors(p, umaxf);
				mpz_set(tz, q);
				append_udivisors(tz, umaxf);
				for (i = 0; i < ndivisors; ++i) {
					if (qmod == (udivisors[i] % ud)) ++c;
				}
			}
		} else {
			init_udivisors();
			append_udivisors(p, umaxf);
			mpz_set(tz, q);
			append_udivisors(tz, umaxf);
			for (i = 0; i < ndivisors; ++i) {
				mpz_add_ui(tz, q, udivisors[i]);
				if (mpz_divisible_p(tz, d)) ++c;
			}
		}
	} else {
		init_zdivisors();
		append_zdivisors(p, cm_maxf);
		mpz_set(tz, q);
		append_zdivisors(tz, cm_maxf);
		for (i = 0; i < ndivisors; ++i) {
			mpz_add(tz, q, zdivisors[i]);
			if (mpz_divisible_p(tz, d)) ++c;
		}
	}
	if (c) {
		if (opt_v) {
			debug_clear();
			print_set(2);
			gmp_printf(": %d (C %Zd%%%Zd)\n", c, q, d);
		}
		mpz_add_ui(total, total, c);
	}
}

/*
	add to total the number of divisors minf <= d <= sqrt(pq) that are == -q (mod d)
*/
void count_min_mod(mpz_t p, mpz_t q, mpz_t d, mpz_t minf) {
	int c = 0, i;
	unsigned long uminf, umaxf;
	mpz_mul(cm_maxf, p, q);
	mpz_sqrt(cm_maxf, cm_maxf);
	mpz_set(tz, q);
	if (mpz_fits_ulong_p(cm_maxf)) {
		if (!mpz_fits_ulong_p(minf)) return;	/* minf > maxf */
		umaxf = mpz_get_ui(cm_maxf);
		uminf = mpz_get_ui(minf);
		init_udivisors();
		append_udivisors(p, umaxf);
		append_udivisors(tz, umaxf);
		for (i = 0; i < ndivisors; ++i) {
			if (uminf > udivisors[i]) continue;
			mpz_add_ui(tz, q, udivisors[i]);
			if (mpz_divisible_p(tz, d)) ++c;
		}
	} else {
		init_zdivisors();
		append_zdivisors(p, cm_maxf);
		append_zdivisors(tz, cm_maxf);
		if (mpz_fits_ulong_p(minf)) {
			uminf = mpz_get_ui(minf);
			for (i = 0; i < ndivisors; ++i) {
				if (mpz_cmp_ui(zdivisors[i], uminf) < 0) continue;
				mpz_add(tz, q, zdivisors[i]);
				if (mpz_divisible_p(tz, d)) ++c;
			}
		} else {
			for (i = 0; i < ndivisors; ++i) {
				if (mpz_cmp(minf, zdivisors[i]) > 0) continue;
				mpz_add(tz, q, zdivisors[i]);
				if (mpz_divisible_p(tz, d)) ++c;
			}
		}
	}
	if (c) {
		if (opt_v) {
			debug_clear();
			print_set(2);
			gmp_printf(": %d (D %Zd%%%Zd > %Zd)\n", c, q, d, minf);
		}
		mpz_add_ui(total, total, c);
	}
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
				/* cmin is low, mark it and enter multiplying phase */
				mpz_set(clow, cmid);
			}
			state = 1;
			mpz_mul_2exp(cmid, cmid, 2);
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
				mpz_mul_2exp(cmid, cmid, 2);
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
		mpz_divexact(p, p, g);
		mpz_divexact(q, q, g);
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
		A118085_r_varstack_t *v = &va[n - 1];
		mpz_set(v->cur, cur);
		mpz_set(v->cmax, cmax);
		for (; mpz_cmp(v->cur, v->cmax) <= 0; mpz_add_ui(v->cur, v->cur, 1)) {
			if (!--debug_counter && !opt_v) {
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
		if (opt_v) mpz_set(va[2].cur, cur);
		if (!--debug_counter && !opt_v) {
			debug_counter = DEBUG_LOOP;
			debug_print(n, p, q, cmax, cur);
		}
		mpz_mul(p2, p, cur);
		mpz_add_ui(q2, cur, 1);
		mpz_mul(q2, q2, q);
		mpz_gcd(g, p2, q2);
		if (mpz_cmp_ui(g, 1) > 0) {
			mpz_divexact(p2, p2, g);
			mpz_divexact(q2, q2, g);
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
			if (mpz_cmp_ui(d, 1) > 0) {
				count_min_mod(p2, q2, d, minf);
			} else {
				count_min(p2, q2, d, minf);
			}
		}
	}
}

int main(int argc, char** argv) {
	int n, arg = 1;
	long clock_tick;
	mpz_t p, q, cur, tempz;
	clock_t t0;

	while (arg < argc && argv[arg][0] == '-') {
		if (strcmp(argv[arg], "-v") == 0) {
			opt_v = 1 - opt_v;
		} else if (strcmp(argv[arg], "--") == 0) {
			++arg;
			break;
		} else if (strcmp(argv[arg], "-?") == 0 || strcmp(argv[arg], "-h") == 0) {
			printf("Usage: %s [ -v ] <n> [ <a_1> <a_2> ... ]\n", argv[0]);
			printf("Calculates A118085(n), the number of n-element multisets with prod{1+1/a_i}=2\nAny specified a_i are fixed in the ordered list of elements\nIf -v is specified, information on matching multisets is also printed\n");
			exit(0);
		} else {
			fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
			exit(-1);
		}
		++arg;
	}
	if (!(arg < argc)) {
		fprintf(stderr, "Usage: %s [ -v ] <n> [ <a_1> <a_2> ... ]\n", argv[0]);
		exit(-1);
	}
	n = atoi(argv[arg++]);
	if (n < 1 || n > MAXN) {
		fprintf(stderr, "Error: 1 <= n <= %d required (got n = %d)\n", MAXN, n);
		exit(-1);
	}
	if (!opt_v) {
		if (setvbuf(stdout, (char*)NULL, _IONBF, (size_t)0)) {
			fprintf(stderr, "setvbuf failed\n");
			exit(-1);
		}
	}
	initp();
	initz();
	mpz_init_set_ui(p, 2);
	mpz_init_set_ui(q, 1);
	mpz_init_set_ui(cur, 0);
	mpz_init(tempz);
	van = n;
	while (arg < argc) {
		if (mpz_set_str(tempz, argv[arg], 10) != 0) {
			fprintf(stderr, "Error parsing argument '%s'\n", argv[arg]);
			exit(-1);
		}
		if (mpz_cmp(cur, tempz) > 0) {
			fprintf(stderr, "Argument '%s' is smaller than a previous value\n", argv[arg]);
			exit(-1);
		}
		mpz_set(cur, tempz);
		mpz_mul(p, p, cur);
		mpz_add_ui(tempz, tempz, 1);
		mpz_mul(q, q, tempz);
		mpz_set(va[n - 1].cur, cur);
		--n;
		++arg;
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

gmseq 3 - result: 5 (0.00)
gmseq 4 - result: 43 (0.00)
gmseq 5 - result: 876 (0.00)
gmseq 6 - result: 49513 (0.06)
gmseq 7 - result: 13005235 (171.21)
gmseq 8 - (est. 12000000)

gmseq 8 11 - result: 0 (0.00)
gmseq 8 10 - result: 0 (0.00)
gmseq 8 9 - result: 9 (0.00)
gmseq 8 8 - result: 130 (0.01)
gmseq 8 7 - result: 1683 (0.03)
gmseq 8 6 - result: 59606 (0.34)
gmseq 8 5 - result: 1470772 (14.11)
gmseq 8 4 - (est. 360000)
gmseq 8 3 - (est. 1250000)
gmseq 8 2 - (est. 10500000)

gmseq 8 4 14 - result: 0 (0.00)
gmseq 8 4 13 - result: 0 (0.00)
gmseq 8 4 12 - result: 5 (0.00)
gmseq 8 4 11 - result: 17 (0.00)
gmseq 8 4 10 - result: 345 (0.00)
gmseq 8 4 9 - result: 1159 (0.03)
gmseq 8 4 8 - result: 7169 (0.08)
gmseq 8 4 7 - result: 105009 (0.43)
gmseq 8 4 6 - result: 3881067 (560.24)
gmseq 8 4 5 - result: 3250390 (15.14)
gmseq 8 4 4 - (est. 360000)

gmseq 8 3 16 - result: 0 (0.00)
gmseq 8 3 15 - result: 2 (0.00)
gmseq 8 3 14 - result: 12 (0.00)
gmseq 8 3 13 - result: 28 (0.01)
gmseq 8 3 12 - result: 263 (0.00)
gmseq 8 3 11 - result: 708 (0.03)
gmseq 8 3 10 - result: 5616 (0.06)
gmseq 8 3 9 - result: 29605 (0.23)
gmseq 8 3 8 - result: 251077 (0.80)
gmseq 8 3 7 - result: 3443209 (401.17)
gmseq 8 3 6 - result: 18632445 (1504.26)
gmseq 8 3 5 - result: 55759629 (3316.43)
gmseq 8 3 4 - result: 159867712 (30100.49)
gmseq 8 3 3 - (est. 1200000)

gmseq 8 2 23 - result: 0 (0.00)
gmseq 8 2 22 - result: 0 (0.00)
gmseq 8 2 21 - result: 3 (0.00)
gmseq 8 2 20 - result: 13 (0.00)
gmseq 8 2 19 - result: 8 (0.00)
gmseq 8 2 18 - result: 232 (0.01)
gmseq 8 2 17 - result: 236 (0.05)
gmseq 8 2 16 - result: 1201 (0.18)
gmseq 8 2 15 - result: 7452 (0.21)
gmseq 8 2 14 - result: 19715 (0.69)
gmseq 8 2 13 - result: 71988 (1.47)
gmseq 8 2 12 - result: 234759 (3.53)
gmseq 8 2 11 - result: 907158 (7.83)
gmseq 8 2 10 - result: 10998677 (3562.45)
gmseq 8 2 9 - result: 18352365 (569.73)
gmseq 8 2 8 - (est. 8400)
gmseq 8 2 7 - (est. 34000)
gmseq 8 2 6 - (est. 170000)
gmseq 8 2 5 - (est. 1700000)
gmseq 8 2 4 - (est. 8400000)

*/
