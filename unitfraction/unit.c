#include "unit.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

/* Search recursively for representations of p/q as a sum of <depth> distinct
   unit fractions. */

uint maxdepth;

/* a prime power p^k, used for representing a factorization */
typedef struct fac_s {
    uint p;
    uint k;
} fac_t;

fac_t* f;       /* factors of q */
uint fcount;    /* number of factors in f */
uint fsize;     /* malloced size of f */

mpz_t* divs = NULL;     /* divisors of some q^2 */
uint dcount = 0;        /* number of divisors */
uint dsize = 0;         /* malloced size of divs */

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
} rat_t;
rat_t *rat;

#define RI(ri) rat[ri].r
#define QI(ri) mpq_denref(RI(ri))
#define PI(ri) mpq_numref(RI(ri))
#define MINI(ri) mpq_denref(rat[ri].qmin)
#define MAXI(ri) rat[ri].max

mpz_t f2_mod, f2_min, ztemp;
mpq_t qtemp;

double t0 = 0;
struct rusage rusage_buf;
static inline double seconds(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000
            - t0;
}

void resize_fac(uint size) {
    f = (fac_t *)realloc(f, size * sizeof(fac_t));
    fsize = size;
}

void resize_div(uint size) {
    divs = (mpz_t *)realloc(divs, size * sizeof(mpz_t));
    for (uint i = dsize; i < size; ++i)
        ZINIT(divs[i]);
    dsize = size;
}

/* print details of the structure rat[ri] to stderr */
void diagnose(uint ri) {
    rat_t* r = &rat[ri];
    gmp_fprintf(stderr, "rat[%d]: r = %Qd; qmin = %Qd; max = %Zd\n",
            ri, r->r, r->qmin, r->max);
    gmp_fprintf(stderr, "\n");
}

void init_unit(uint max) {
    t0 = seconds();
    maxdepth = max;
    rat = malloc(maxdepth * sizeof(rat_t));
    for (uint i = 0; i < maxdepth; ++i) {
        QINIT(rat[i].r);
        QINIT(rat[i].qmin);
        mpz_set_ui(mpq_numref(rat[i].qmin), 1);    /* never changes */
        ZINIT(rat[i].max);
    }
    mpz_set_ui(mpq_denref(rat[0].qmin), 0);     /* base, never changes */
    fsize = 0;
    f = NULL;
    resize_fac(10);
    dsize = 0;
    resize_div(1024);
    ZINIT(f2_mod);
    ZINIT(f2_min);
    ZINIT(ztemp);
    QINIT(qtemp);
}

void done_unit(void) {
    for (uint i = 0; i < maxdepth; ++i) {
        QCLEAR(rat[i].r);
        QCLEAR(rat[i].qmin);
        ZCLEAR(rat[i].max);
    }
    free(rat);
    free(f);
    for (uint i = 0; i < dsize; ++i)
        ZCLEAR(divs[i]);
    free(divs);
    ZCLEAR(f2_mod);
    ZCLEAR(f2_min);
    ZCLEAR(ztemp);
    QCLEAR(qtemp);
}

/* set rat[ri].f to the factors of QI(ri) */
void factorize(uint ri) {
    uint ni = 0;
    mpz_set(ztemp, QI(ri));
    if (mpz_even_p(ztemp)) {
        f[ni].p = 2;
        f[ni].k = 0;
        while (mpz_even_p(ztemp)) {
            ++f[ni].k;
            mpz_divexact_ui(ztemp, ztemp, 2);
        }
        ++ni;
    }
    uint p = 3;
    while (mpz_cmp_ui(ztemp, 1) > 0) {
        if (mpz_divisible_ui_p(ztemp, p)) {
            if (ni + 1 >= fsize)
                resize_fac(fsize + 8);
            f[ni].p = p;
            f[ni].k = 0;
            while (mpz_divisible_ui_p(ztemp, p)) {
                ++f[ni].k;
                mpz_divexact_ui(ztemp, ztemp, p);
            }
            ++ni;
        }
        p += 2;
    }
    fcount = ni;
    return;
}

/* Set divs to the divisors of QI(ri)^2.
 */
void divisors(uint ri) {
    factorize(ri);
    uint count = 1, prev, k, next;
    uint p;
    for (uint i = 0; i < fcount; ++i)
        count *= 2 * f[i].k + 1;
    if (count > dsize)
        resize_div(count * 2);

    mpz_set_ui(divs[0], 1);
    prev = 1;
    for (uint i = 0; i < fcount; ++i) {
        p = f[i].p;
        k = f[i].k * 2;
        next = prev * k;
        for (uint j = 0; j < next; ++j)
            mpz_mul_ui(divs[j + prev], divs[j], p);
        prev += next;
    }
    dcount = count;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators > m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q.
 */
bool find_s2(uint ri) {
    divisors(ri);
    uint end_div = dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (uint i = 0; i < end_div; ++i) {
        if (mpz_cmp(f2_min, divs[i]) >= 0)
            continue;
        if (mpz_cmp(QI(ri), divs[i]) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, divs[i], PI(ri)))
            return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators >= m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q <= d < q.
 */
bool find_m2(uint ri) {
    divisors(ri);
    uint end_div = dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (uint i = 0; i < end_div; ++i) {
        if (mpz_cmp(f2_min, divs[i]) > 0)
            continue;
        if (mpz_cmp(QI(ri), divs[i]) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, divs[i], PI(ri)))
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
        if ((depth == 3) ? find_m2(rj) : find_mn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return TRUE if p/q can be represented as fewer than c distinct unit
 * fractions.
 */
bool better_set(mpq_t r, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
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
bool better_multi(mpq_t r, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_m2(0))
        return 1;
    if (c == 3)
        return 0;
    if (find_mn(0, c - 1))
        return 1;
    return 0;
}

uint find_set(mpq_t r) {
    if (mpq_sgn(r) == 0)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_s2(0))
        return 2;
    for (uint c = 3; 1; ++c) {
        if (find_sn(0, c))
            return c;
    }
    /* not reached */
}

/* return floor(r), assumed to fit in uint */
uint mpq_int(mpq_t r) {
    if (mpz_cmp(mpq_numref(r), mpq_denref(r)) < 0)
        return 0;
    mpz_fdiv_q(ztemp, mpq_numref(r), mpq_denref(r));
    return mpz_get_ui(ztemp);
}

uint find_multi(mpq_t r) {
    if (mpq_sgn(r) == 0)
        return 0;
    uint unit = mpq_int(r);
    if (unit) {
        mpz_fdiv_r(mpq_numref(r), mpq_numref(r), mpq_denref(r));
        if (mpq_sgn(r) == 0)
            return unit;
    }
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return unit + 1;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_m2(0))
        return unit + 2;
    for (uint c = 3; 1; ++c) {
        if (find_mn(0, c))
            return unit + c;
    }
    /* not reached */
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * square unit fractions with denominators > (m = MINI(ri))^2.
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q such that (q + d)/p and (q + q^2/d)/p
 * are both perfect squares.
 */
bool find_square_s2(uint ri) {
    divisors(ri);
    uint end_div = dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_mul(f2_min, f2_min, MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (uint i = 0; i < end_div; ++i) {
        if (mpz_cmp(f2_min, divs[i]) >= 0)
            continue;
        if (mpz_cmp(QI(ri), divs[i]) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (!mpz_congruent_p(f2_mod, divs[i], PI(ri)))
            continue;
        mpz_add(ztemp, QI(ri), divs[i]);
        mpz_divexact(ztemp, ztemp, PI(ri));
        if (!mpz_perfect_square_p(ztemp))
            continue;
        mpz_mul(ztemp, QI(ri), QI(ri));
        mpz_divexact(ztemp, ztemp, divs[i]);
        mpz_add(ztemp, ztemp, QI(ri));
        mpz_divexact(ztemp, ztemp, PI(ri));
        if (!mpz_perfect_square_p(ztemp))
            continue;
gmp_printf("success with d=%Zi for q=%Qi\n", divs[i], RI(ri));
        return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth> distinct
 * square unit fractions with denominators > (m = MINI(ri))^2.
 * Requires depth >= 3.
 */
bool find_square_sn(uint ri, uint depth) {
    uint rj = ri + 1;
/*
    if (mpz_cmp_ui(PI(ri), 1) == 0
        && mpz_perfect_square_p(QI(ri))
        ... and QI(ri) > MINI(ri)^2
    )
        return 1;
*/

    /* max := floor(sqrt(q * depth / p)) */
    mpz_mul_ui(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_sqrt(MAXI(rj), MAXI(rj));

    /* min := max(ceil(sqrt(q / p)), prevmin + 1) */
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    mpz_sqrtrem(MINI(rj), ztemp, MINI(rj));
    if (mpz_sgn(ztemp) > 0)
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    if (mpz_cmp(MINI(rj), MINI(ri)) <= 0)
        mpz_add_ui(MINI(rj), MINI(ri), 1);

    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_mul(qtemp, rat[rj].qmin, rat[rj].qmin);
        mpq_sub(RI(rj), RI(ri), qtemp);
        if ((depth == 3) ? find_square_s2(rj) : find_square_sn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

uint find_square_set(mpq_t r, uint max_depth) {
    if (mpq_sgn(r) == 0)
        return 0;
    gmp_printf("%Qu: not 0 (%.2fs)\n", r, seconds());
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0
        && mpz_perfect_square_p(mpq_denref(r))
    )
        return 1;
    gmp_printf("%Qu: not 1 (%.2fs)\n", r, seconds());
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_square_s2(0))
        return 2;
    gmp_printf("%Qu: not 2 (%.2fs)\n", r, seconds());
    for (uint c = 3; c <= max_depth; ++c) {
        if (find_square_sn(0, c)) {
            for (uint i = 0; i <= c - 2; ++i)
                gmp_printf("%Zu (%Qd); ", MINI(i), RI(i));
            printf(" (%.2fs)\n", seconds());
            return c;
        }
        gmp_printf("%Qu: not %u (%.2fs)\n", r, c, seconds());
    }
    return 0;
    /* not reached */
}

