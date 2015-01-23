#include <gmp.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/times.h>

/* Search recursively for representations of p/q as a sum of <depth> distinct
   unit fractions. */

/* handy macros in case we need to debug our GMP memory use */
#define ZINIT(z) mpz_init(z)
#define ZCLEAR(z) mpz_clear(z)
#define QINIT(q) mpq_init(q)
#define QCLEAR(q) mpq_clear(q)

/* we are searching for depth = 8, so 16 should be more than we'll ever need */
#define MAXDEPTH 16

long clock_tick;    /* ticks per second */

/* a prime p, used for factorizing */
typedef struct prime_s {
    unsigned int p;            /* the prime itself */
    unsigned int p_squared;    /* p^2, or 0 if p^2 > MAXUINT */
    mpz_t zp_squared;        /* p^2 */
} prime_t;

prime_t* primes;    /* ordered list of known primes */
int primecount;        /* number of entries in <primes> */
int primesize;        /* malloced size of <primes> */

/* a prime power p^k, used for representing a factorization */
typedef struct fac_s {
    unsigned int p;
    int k;
} fac_t;

/* Structure recording information about a rational p/q at a given depth in
 * the recursion.
 * Note that we always set rat[ri].r = rat[ri-1].r - rat[ri].qmin: both
 * qmin and max in rat[ri] refer to the range of unit fractions to be
 * trial-subtracted from rat[ri-1].r. We leave rat[0].qmin = 1/0 to get the
 * correct limits.
 */
typedef struct rat_s {
    mpq_t r;    /* the rational itself, p/q */
    mpq_t qmin;    /* 1/min, the last unit fraction to be trial-subtracted */
    mpz_t max;    /* final den. for trial subtraction */
    fac_t* f;    /* factors of q */
    int fcount;    /* number of factors in f */
    int fsize;    /* malloced size of f */
    mpz_t* d;    /* divisors of q^2 */
    int dcount;    /* number of divisors in d */
    int dsize;    /* malloced size of d */
} rat_t;
rat_t rat[MAXDEPTH];

#define RI(ri) rat[ri].r
#define QI(ri) mpq_denref(RI(ri))
#define PI(ri) mpq_numref(RI(ri))
#define MINI(ri) mpq_denref(rat[ri].qmin)
#define MAXI(ri) rat[ri].max

int vec[MAXDEPTH];    /* record required depths for a given <q> */

mpz_t factor_n, factor_q, factor_r;
mpz_t f2_mod, f2_min;

void resize_primes(int size) {
    int i;
    primes = (prime_t *)realloc(primes, size * sizeof(prime_t));
    for (i = primesize; i < size; ++i)
        ZINIT(primes[i].zp_squared);
    primesize = size;
}

void resize_fac(int ri, int size) {
    rat[ri].f = (fac_t *)realloc(rat[ri].f, size * sizeof(fac_t));
    rat[ri].fsize = size;
}

void resize_div(int ri, int size) {
    int i;
    rat[ri].d = (mpz_t *)realloc(rat[ri].d, size * sizeof(mpz_t));
    for (i = rat[ri].dsize; i < size; ++i)
        ZINIT(rat[ri].d[i]);
    rat[ri].dsize = size;
}

/* print details of the structure rat[ri] to stderr */
void diagnose(int ri) {
    rat_t* r = &rat[ri];
    int i;
    gmp_fprintf(stderr, "rat[%d]: r = %Qd; qmin = %Qd; max = %Zd\n",
            ri, r->r, r->qmin, r->max);
    gmp_fprintf(stderr, "  f[size=%u]:", r->fsize);
    for (i = 0; i < r->fcount; ++i) {
        if (i) fprintf(stderr, " .");
        gmp_fprintf(stderr, " %d^%d", r->f[i].p, r->f[i].k);
    }
    gmp_fprintf(stderr, "\n  d[size=%u]:", r->dsize);
    for (i =0; i < r->dcount; ++i)
        gmp_fprintf(stderr, " %Zu", r->d[i]);
    gmp_fprintf(stderr, "\n");
}

void init_gmp(void) {
    int i;
    for (i = 0; i < MAXDEPTH; ++i) {
        QINIT(rat[i].r);
        QINIT(rat[i].qmin);
        mpz_set_si(mpq_numref(rat[i].qmin), 1);    /* never changes */
        ZINIT(rat[i].max);
        rat[i].fsize = 0;
        rat[i].f = NULL;
        resize_fac(i, 10);
        rat[i].dsize = 0;
        rat[i].d = NULL;
        resize_div(i, 1024);
    }
    ZINIT(factor_n);
    ZINIT(factor_q);
    ZINIT(factor_r);
    ZINIT(f2_mod);
    ZINIT(f2_min);
    primes = NULL;
    resize_primes(1024);
    primes[0].p = 2;
    primes[0].p_squared = 4;
    mpz_set_si(primes[0].zp_squared, 4);
    primes[1].p = 3;
    primes[1].p_squared = 9;
    mpz_set_si(primes[1].zp_squared, 9);
    primecount = 2;
}

void clear_gmp(void) {
    int i, j;
    for (i = 0; i < MAXDEPTH; ++i) {
        QCLEAR(rat[i].r);
        QCLEAR(rat[i].qmin);
        ZCLEAR(rat[i].max);
        for (j = 0; j < rat[i].dsize; ++j)
            ZCLEAR(rat[i].d[j]);
    }
    ZCLEAR(factor_n);
    ZCLEAR(factor_q);
    ZCLEAR(factor_r);
    ZCLEAR(f2_mod);
    ZCLEAR(f2_min);
    for (i = 0; i < primesize; ++i)
        ZCLEAR(primes[i].zp_squared);
    free(primes);
}

/* time used so far (CPU seconds) */
double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

/* boolean gcd test; assumes 0 <= p < q */
int have_common_divisor(int p, int q) {
    while (p > 0) {
        int t = q % p;
        q = p;
        p = t;
    }
    return (q == 1) ? 0 : 1;
}

/* find the next prime not in primes[], extend the list, and return the
 * index of the new prime */
int nextprime(void) {
    unsigned int p = primes[primecount - 1].p;
    int i, isprime;
    if (primecount == primesize)
        resize_primes(primesize * 3 / 2);
    while (p > 1) {
        p += 2;
        isprime = 1;
        for (i = 1; primes[i].p_squared <= p; ++i) {
            if ((p % primes[i].p) == 0) {
                isprime = 0;
                break;
            }
        }
        if (isprime) {
            primes[primecount].p = p;
            if (p >= 65536) {
                primes[primecount].p_squared = 0;
                mpz_set_si(primes[primecount].zp_squared, p);
                mpz_mul_si(primes[primecount].zp_squared,
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
 * return true if factor_n == 1
 */
int try_div(int ri, unsigned int p) {
    int power;

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
void factorize(int ri) {
    int i, p;
    mpz_set(factor_n, QI(ri));
    rat[ri].fcount = 0;
    if (mpz_cmp_si(factor_n, 1) == 0)
        return;
    if (ri > 0) {
        /* given r[ri] = r[ri-1] - 1/something, try first dividing by the prime
         * factors of QI(ri-1) */
        rat_t* prev = &rat[ri - 1];
        for (i = 0; i < prev->fcount; ++i)
            if (try_div(ri, prev->f[i].p))
                return;
    }
    for (i = 0; 1; ++i) {
        if (i >= primecount)
            i = nextprime();
        if (try_div(ri, primes[i].p))
            return;
        if (mpz_cmp(primes[i].zp_squared, factor_n) > 0) {
            p = mpz_get_si(factor_n);
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
void divisors(int ri) {
    int count = 1, i, j, prev, k, next;
    unsigned int p;
    for (i = 0; i < rat[ri].fcount; ++i)
        count *= 2 * rat[ri].f[i].k + 1;
    if (count > rat[ri].dsize)
        resize_div(ri, count * 2);

    mpz_set_ui(rat[ri].d[0], 1);
    prev = 1;
    for (i = 0; i < rat[ri].fcount; ++i) {
        p = rat[ri].f[i].p;
        k = rat[ri].f[i].k * 2;
        next = prev * k;
        for (j = 0; j < next; ++j)
            mpz_mul_ui(rat[ri].d[j + prev], rat[ri].d[j], p);
        prev += next;
    }
    rat[ri].dcount = count;
}

/* Return true if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators greater than m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q.
 * This assumes divisors(ri) has previously been called.
 */
int find_2(int ri) {
    int i, end_div = rat[ri].dcount;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    for (i = 0; i < end_div; ++i) {
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

/* Return true if r = RI(ri) can be expressed as the sum of <depth> distinct
 * unit fractions with denominators greater than m = MINI(ri).
 * Requires depth >= 3, and that find_n(ri, depth - 1) is false.
 */
int find_n(int ri, int depth) {
    int rj = ri + 1;
    mpz_mul_si(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    if (mpz_cmp(MINI(ri), MINI(rj)) >= 0)
        mpz_add_ui(MINI(rj), MINI(ri), 1);
    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_sub(RI(rj), RI(ri), rat[rj].qmin);
        factorize(rj);
        divisors(rj);
        if ((depth == 3) ? find_2(rj) : find_n(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return the number of distinct unit fractions required to express p/q.
 * Assumes QI(0) = q has been set, and the factors and divisors found.
 */
int find_best(int p, int q) {
    int depth;

    if (p == 1) return 1;
    mpz_set_si(PI(0), p);

    if (find_2(0)) return 2;
    depth = 3;
    while (1) {
        if (find_n(0, depth)) return depth;
        ++depth;
    }
}

int main(int argc, char** argv) {
    int p, q, qend, depth, result, i;
    int arg = 1;
    int opt_u = 0, opt_p = 0;

    while (arg < argc && argv[arg][0] == '-') {
        char* s = argv[arg++];
        if (strcmp(s, "--") == 0)
            break;
        if (strcmp(s, "-q") == 0) {
            opt_u = 0;
            opt_p = 0;
            continue;
        }
        if (strcmp(s, "-u") == 0) {
            opt_u = 1;
            opt_p = 0;
            continue;
        }
        if (strcmp(s, "-p") == 0) {
            opt_u = 1;
            opt_p = 1;
            continue;
        }
        fprintf(stderr, "Unknown option '%s'\n", s);
        argc = -1;  /* force usage message */
        break;
    }

    if (argc - arg != 3) {
        fprintf(stderr, "Usage: %s [ -q <start> <end> <minsize>\n", argv[0]);
        return(-1);
    }

    init_gmp();
    q = atoi(argv[arg++]) - 1;
    qend = atoi(argv[arg++]);
    depth = atoi(argv[arg++]);
    mpz_set_si(MINI(0), 0);    /* min_0 never changes */

    clock_tick = sysconf(_SC_CLK_TCK);
    setvbuf(stdout, NULL, _IOLBF, 0);

    while (++q <= qend) {
        if (!opt_u)
            memset(vec, 0, (depth + 1) * sizeof(int));
        mpz_set_si(QI(0), q);
        factorize(0);
        if (opt_p && (rat[0].fcount != 1 || rat[0].f[0].k != 1)) {
            /* not a prime, skip if "-p" */
            continue;
        }
        divisors(0);
        for (p = (opt_u ? q - 1 : 1); p < q; ++p) {
            if (!opt_u && have_common_divisor(p, q))
                continue;
            result = find_best(p, q);
            if (result > depth) {
                printf("*** f(%d) = %d / %d ***\n", result, p, q);
                if (!opt_u) {
                    if (depth > MAXDEPTH) {
                        fprintf(stderr, "%d exceeds max depth, stopping\n", depth);
                        return(0);
                    }
                    for (i = depth; i < result; ++i)
                        vec[i + 1] = 0;
                }
                depth = result;
            }
            if (!opt_u)
                ++vec[result];
        }
        if (opt_u) {
            printf("%d: %d [%.2fs]\n", q, result, timing());
        } else {
            printf("%d:", q);
            for (i = 1; i <= depth; ++i)
                printf(" %d", vec[i]);
            printf(" [%.2fs]\n", timing());
        }
    }
    clear_gmp();
    return 0;
}
