#include <stdlib.h>
#include "coulfact.h"
#include "gmp_main.h"   /* prime_count */

/* for sorting */
int _mpz_comparator(const void *va, const void *vb) {
    return mpz_cmp(*(mpz_t *)va, *(mpz_t *)vb);
}

void init_fact(t_fact *f) {
    f->count = 0;
    f->size = 16;
    f->ppow = (t_ppow *)malloc(f->size * sizeof(t_ppow));
}
void free_fact(t_fact *f) {
    free(f->ppow);
}
void add_fact(t_fact *f, t_ppow pp) {
    uint count = f->count++;
	if (f->count > f->size) {
		uint size = f->size * 2;
		f->ppow = (t_ppow *)realloc(f->ppow, size * sizeof(t_ppow));
		f->size = size;
	}
	f->ppow[count] = pp;
}
void reverse_fact(t_fact *f) {
    t_ppow pp;
    uint c = f->count;
    for (uint i = 0; i + i < c - 1; ++i) {
        uint j = c - i - 1;
        pp = f->ppow[i];
        f->ppow[i] = f->ppow[j];
        f->ppow[j] = pp;
    }
}

void init_zfact(t_zfact *f) {
    f->count = 0;
    f->size = 16;
    f->ppow = (t_zpow *)malloc(f->size * sizeof(t_zpow));
	for (int i = 0; i < f->size; ++i)
		mpz_init(f->ppow[i].p);
}
void free_zfact(t_zfact *f) {
	for (int i = 0; i < f->size; ++i)
		mpz_clear(f->ppow[i].p);
	free(f->ppow);
}
void add_zfact(t_zfact *f, t_zpow pp) {
    uint count = f->count++;
	if (f->count > f->size) {
        uint size = f->size * 2;
        f->ppow = (t_zpow *)realloc(f->ppow, size * sizeof(t_zpow));
		for (int i = f->count; i < size; ++i)
			mpz_init(f->ppow[i].p);
        f->size = size;
    }
	mpz_set(f->ppow[count].p, pp.p);
	f->ppow[count].e = pp.e;
}

uint try_simple_fact(uint n, uint d, t_fact *f) {
	uint e = 0;
	while ((n % d) == 0) {
		n /= d;
		++e;
	}
	if (e) {
		t_ppow pp;
		pp.p = d;
		pp.e = e;
		add_fact(f, pp);
	}
	return n;
}

void simple_fact(uint n, t_fact *f) {
	uint d = 3;
	if (n > 1)
		n = try_simple_fact(n, 2, f);
	while (n > 1) {
		n = try_simple_fact(n, d, f);
		d += 2;
	}
	return;
}

uint simple_tau(t_fact *f) {
	uint t = 1;
	for (uint i = 0; i < f->count; ++i)
		t *= f->ppow[i].e + 1;
	return t;
}

uint simple_valuation(ulong n, ulong p) {
    uint v = 0;
    while ((n % p) == 0) {
        ++v;
        n /= p;
    }
    return v;
}

uint simple_prime_count(ulong n) {
    mpz_t zn, zc;
    mpz_init_set_ui(zn, n);
    mpz_init(zc);
    prime_count(zc, zn);
    uint c = mpz_get_ui(zc);
    mpz_clear(zn);
    mpz_clear(zc);
    return c;
}

uint tiny_gcd(uint a, uint b) {
    if (a > b)
        return tiny_gcd(b, a);
    if (a == 0)
        return b;
    return tiny_gcd(b % a, a);
}

ulong simple_gcd(ulong a, ulong b) {
    if (a > b)
        return simple_gcd(b, a);
    if (a == 0)
        return b;
    return simple_gcd(b % a, a);
}
