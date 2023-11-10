#include "unit.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

/* Search recursively for representations of p/q as a sum of <depth> distinct
   unit fractions. */

uint maxdepth;

/* a prime p, used for factorizing */
typedef struct prime_s {
    uint p;             /* the prime itself */
    uint p_squared;     /* p^2, or 0 if p^2 > MAXUINT */
    mpz_t zp_squared;   /* p^2 */
} prime_t;

prime_t* primes;    /* ordered list of known primes */
uint primecount;    /* number of entries in <primes> */
uint primesize;     /* malloced size of <primes> */

/* a prime power p^k, used for representing a factorization */
typedef struct fac_s {
    uint p;
    uint k;
} fac_t;

/* Structure recording information about a rational p/q at a given depth in
 * the recursion.
 * Note that we always set rat[ri].r = rat[ri-1].r - rat[ri].qmin: both
 * qmin and max in rat[ri] refer to the range of unit fractions to be
 * trial-subtracted from rat[ri-1].r. We leave rat[0].qmin = 1/0 to get the
 * correct limits.
 */
typedef struct rat_s {
    mpq_t r;        /* the rational itself, p/q */
    mpq_t qmin;     /* 1/min, the last unit fraction to be trial-subtracted */
    mpz_t max;      /* final den. for trial subtraction */
    fac_t* f;       /* factors of q */
    uint fcount;    /* number of factors in f */
    uint fsize;     /* malloced size of f */
    mpz_t* d;       /* divisors of q^2 */
    uint dcount;    /* number of divisors in d */
    uint dsize;     /* malloced size of d */
} rat_t;
rat_t *rat;

#define RI(ri) rat[ri].r
#define QI(ri) mpq_denref(RI(ri))
#define PI(ri) mpq_numref(RI(ri))
#define MINI(ri) mpq_denref(rat[ri].qmin)
#define MAXI(ri) rat[ri].max

mpz_t factor_n, factor_q, factor_r;
mpz_t f2_mod, f2_min, ztemp;

void resize_primes(uint size) {
    primes = (prime_t *)realloc(primes, size * sizeof(prime_t));
    for (uint i = primesize; i < size; ++i)
        ZINIT(primes[i].zp_squared);
    primesize = size;
}

void resize_fac(uint ri, uint size) {
    rat[ri].f = (fac_t *)realloc(rat[ri].f, size * sizeof(fac_t));
    rat[ri].fsize = size;
}

void resize_div(uint ri, uint size) {
    rat[ri].d = (mpz_t *)realloc(rat[ri].d, size * sizeof(mpz_t));
    for (uint i = rat[ri].dsize; i < size; ++i)
        ZINIT(rat[ri].d[i]);
    rat[ri].dsize = size;
}

/* print details of the structure rat[ri] to stderr */
void diagnose(uint ri) {
    rat_t* r = &rat[ri];
    gmp_fprintf(stderr, "rat[%d]: r = %Qd; qmin = %Qd; max = %Zd\n",
            ri, r->r, r->qmin, r->max);
    gmp_fprintf(stderr, "  f[size=%u]:", r->fsize);
    for (uint i = 0; i < r->fcount; ++i) {
        if (i)
            fprintf(stderr, " .");
        gmp_fprintf(stderr, " %d^%d", r->f[i].p, r->f[i].k);
    }
    gmp_fprintf(stderr, "\n  d[size=%u]:", r->dsize);
    for (uint i =0; i < r->dcount; ++i)
        gmp_fprintf(stderr, " %Zu", r->d[i]);
    gmp_fprintf(stderr, "\n");
}

void init_unit(uint max) {
    maxdepth = max;
    rat = malloc(maxdepth * sizeof(rat_t));
    for (uint i = 0; i < maxdepth; ++i) {
        QINIT(rat[i].r);
        QINIT(rat[i].qmin);
        mpz_set_ui(mpq_numref(rat[i].qmin), 1);    /* never changes */
        ZINIT(rat[i].max);
        rat[i].fsize = 0;
        rat[i].f = NULL;
        resize_fac(i, 10);
        rat[i].dsize = 0;
        rat[i].d = NULL;
        resize_div(i, 1024);
    }
    mpz_set_ui(mpq_denref(rat[0].qmin), 0);     /* base, never changes */
    ZINIT(factor_n);
    ZINIT(factor_q);
    ZINIT(factor_r);
    ZINIT(f2_mod);
    ZINIT(f2_min);
    ZINIT(ztemp);
    primes = NULL;
    resize_primes(1024);
    primes[0].p = 2;
    primes[0].p_squared = 4;
    mpz_set_ui(primes[0].zp_squared, 4);
    primes[1].p = 3;
    primes[1].p_squared = 9;
    mpz_set_ui(primes[1].zp_squared, 9);
    primecount = 2;
}

void done_unit(void) {
    for (uint i = 0; i < maxdepth; ++i) {
        QCLEAR(rat[i].r);
        QCLEAR(rat[i].qmin);
        ZCLEAR(rat[i].max);
        for (uint j = 0; j < rat[i].dsize; ++j)
            ZCLEAR(rat[i].d[j]);
        free(rat[i].f);
    }
    free(rat);
    ZCLEAR(factor_n);
    ZCLEAR(factor_q);
    ZCLEAR(factor_r);
    ZCLEAR(f2_mod);
    ZCLEAR(f2_min);
    ZCLEAR(ztemp);
    for (uint i = 0; i < primesize; ++i)
        ZCLEAR(primes[i].zp_squared);
    free(primes);
}

/* boolean gcd test; assumes 0 <= p < q */
bool have_common_divisor(uint p, uint q) {
    while (p > 0) {
        uint t = q % p;
        q = p;
        p = t;
    }
    return (q == 1) ? 0 : 1;
}

/* find the next prime not in primes[], extend the list, and return the
 * index of the new prime */
uint nextprime(void) {
    uint p = primes[primecount - 1].p;
    bool isprime;
    if (primecount == primesize)
        resize_primes(primesize * 3 / 2);
    while (p > 1) {
        p += 2;
        isprime = 1;
        for (uint i = 1; primes[i].p_squared <= p; ++i) {
            if ((p % primes[i].p) == 0) {
                isprime = 0;
                break;
            }
        }
        if (isprime) {
            primes[primecount].p = p;
            if (p > 65535) {
                primes[primecount].p_squared = 0;
                mpz_set_ui(primes[primecount].zp_squared, p);
                mpz_mul_ui(primes[primecount].zp_squared,
                        primes[primecount].zp_squared, p);
            } else {
                primes[primecount].p_squared = p * p;
                mpz_set_ui(primes[primecount].zp_squared,
                        primes[primecount].p_squared);
            }
            return primecount++;
        }
    }
    fprintf(stderr, "Overflow: prime > MAXUINT required\n");
    exit(-1);
}

/* given factor_n = QI(ri), if p^k divides factor_n:
 *   factor_n := factor_n / p^k
 *   append p^k to rat[ri].f[]
 * return TRUE if factor_n == 1
 */
bool try_div(uint ri, uint p) {
    uint power;

    mpz_fdiv_qr_ui(factor_q, factor_r, factor_n, p);
    if (mpz_cmp_ui(factor_r, 0) != 0)
        return 0;    /* p is not a factor */

    power = 1;
    while (1) {
        mpz_set(factor_n, factor_q);
        mpz_fdiv_qr_ui(factor_q, factor_r, factor_n, p);
        if (mpz_cmp_ui(factor_r, 0) != 0)
            break;
        ++power;
    }
    if (rat[ri].fcount >= rat[ri].fsize)
        resize_fac(ri, rat[ri].fsize + 8);
    rat[ri].f[rat[ri].fcount].p = p;
    rat[ri].f[rat[ri].fcount].k = power;
    ++rat[ri].fcount;
    return (mpz_cmp_ui(factor_n, 1) > 0) ? 0 : 1;
}

/* set rat[ri].f to the factors of QI(ri) */
void factorize(uint ri) {
    uint p;
    mpz_set(factor_n, QI(ri));
    rat[ri].fcount = 0;
    if (mpz_cmp_ui(factor_n, 1) == 0)
        return;
    if (ri > 0) {
        /* given r[ri] = r[ri-1] - 1/something, try first dividing by the prime
         * factors of QI(ri-1) */
        rat_t* prev = &rat[ri - 1];
        for (uint i = 0; i < prev->fcount; ++i)
            if (try_div(ri, prev->f[i].p))
                return;
    }
    for (uint i = 0; 1; ++i) {
        if (i >= primecount)
            i = nextprime();
        if (try_div(ri, primes[i].p))
            return;
        if (mpz_cmp(primes[i].zp_squared, factor_n) > 0) {
            p = mpz_get_ui(factor_n);
            if (p == 0) {
                gmp_fprintf(stderr, "Overflow: prime %Zu > MAXUINT\n",
                        factor_n);
                exit(-1);
            }
            try_div(ri, p);
            return;
        }
    }
}

/* Set rat[ri].d to the  divisors of QI(ri)^2.
 * This assumes factorize(ri) has previously been called.
 */
void divisors(uint ri) {
    uint count = 1, prev, k, next;
    uint p;
    for (uint i = 0; i < rat[ri].fcount; ++i)
        count *= 2 * rat[ri].f[i].k + 1;
    if (count > rat[ri].dsize)
        resize_div(ri, count * 2);

    mpz_set_ui(rat[ri].d[0], 1);
    prev = 1;
    for (uint i = 0; i < rat[ri].fcount; ++i) {
        p = rat[ri].f[i].p;
        k = rat[ri].f[i].k * 2;
        next = prev * k;
        for (uint j = 0; j < next; ++j)
            mpz_mul_ui(rat[ri].d[j + prev], rat[ri].d[j], p);
        prev += next;
    }
    rat[ri].dcount = count;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators > m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q.
 * This assumes divisors(ri) has previously been called.
 */
bool find_s2(uint ri) {
    uint end_div = rat[ri].dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (uint i = 0; i < end_div; ++i) {
        if (mpz_cmp(f2_min, rat[ri].d[i]) >= 0)
            continue;
        if (mpz_cmp(QI(ri), rat[ri].d[i]) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, rat[ri].d[i], PI(ri)))
            return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators >= m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q <= d < q.
 * This assumes divisors(ri) has previously been called.
 */
bool find_m2(uint ri) {
    uint end_div = rat[ri].dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (uint i = 0; i < end_div; ++i) {
        if (mpz_cmp(f2_min, rat[ri].d[i]) > 0)
            continue;
        if (mpz_cmp(QI(ri), rat[ri].d[i]) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, rat[ri].d[i], PI(ri)))
            return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth> distinct
 * unit fractions with denominators > m = MINI(ri).
 * Requires depth >= 3.
 */
bool find_sn(uint ri, uint depth) {
    uint rj = ri + 1;
    if (mpz_cmp_ui(PI(ri), 1) == 0)
        return 1;
    mpz_mul_ui(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    if (mpz_cmp(MINI(ri), MINI(rj)) >= 0)
        mpz_add_ui(MINI(rj), MINI(ri), 1);
    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_sub(RI(rj), RI(ri), rat[rj].qmin);
        factorize(rj);
        divisors(rj);
        if ((depth == 3) ? find_s2(rj) : find_sn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth>
 * unit fractions, not necessarily distinct, with denominators >= m = MINI(ri).
 * Requires depth >= 3.
 */
bool find_mn(uint ri, uint depth) {
    uint rj = ri + 1;
    if (mpz_cmp_ui(PI(ri), 1) == 0)
        return 1;
    mpz_mul_ui(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    if (mpz_cmp(MINI(ri), MINI(rj)) > 0)
        mpz_set(MINI(rj), MINI(ri));
    if (mpz_sgn(MINI(rj)) == 0)
        mpz_set_ui(MINI(rj), 1);
    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_sub(RI(rj), RI(ri), rat[rj].qmin);
        factorize(rj);
        divisors(rj);
        if ((depth == 3) ? find_m2(rj) : find_mn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return TRUE if p/q can be represented as fewer than c distinct unit
 * fractions.
 */
bool better_set(mpq_t q, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(q), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), q);
    mpq_get_den(QI(0), q);
    factorize(0);
    if (find_s2(0))
        return 1;
    if (c == 3)
        return 0;
    if (find_sn(0, c - 1))
        return 1;
    return 0;
}

/* Return TRUE if p/q can be represented as fewer than c unit fractions
 * (not necessarily distinct).
 */
bool better_multi(mpq_t q, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(q), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), q);
    mpq_get_den(QI(0), q);
    factorize(0);
    if (find_m2(0))
        return 1;
    if (c == 3)
        return 0;
    if (find_mn(0, c - 1))
        return 1;
    return 0;
}

uint find_set(mpq_t q) {
    if (mpq_sgn(q) == 0)
        return 0;
    if (mpz_cmp_ui(mpq_numref(q), 1) == 0)
        return 1;
    mpq_get_num(PI(0), q);
    mpq_get_den(QI(0), q);
    factorize(0);
    if (find_s2(0))
        return 2;
    for (uint c = 3; 1; ++c) {
        if (find_sn(0, c))
            return c;
    }
    /* not reached */
}

uint mpq_int(mpq_t q) {
    if (mpz_cmp(mpq_numref(q), mpq_denref(q)) < 0)
        return 0;
    mpz_fdiv_q(ztemp, mpq_numref(q), mpq_denref(q));
    return mpz_get_ui(ztemp);
}

uint find_multi(mpq_t q) {
    if (mpq_sgn(q) == 0)
        return 0;
    uint unit = mpq_int(q);
    if (unit) {
        mpz_fdiv_r(mpq_numref(q), mpq_numref(q), mpq_denref(q));
        if (mpq_sgn(q) == 0)
            return unit;
    }
    if (mpz_cmp_ui(mpq_numref(q), 1) == 0)
        return unit + 1;
    mpq_get_num(PI(0), q);
    mpq_get_den(QI(0), q);
    factorize(0);
    if (find_m2(0))
        return unit + 2;
    for (uint c = 3; 1; ++c) {
        if (find_mn(0, c))
            return unit + c;
    }
    /* not reached */
}
