#include <gmp.h>
#include "ptypes.h"

#include "coul.h"
#include "coultau.h"
#include "factor.h"
#include "primality.h"
#include "prime_iterator.h"
#include "utility.h"
#include "pbrent63.h"
#include "squfof126.h"
#include "ecm.h"
#include "tinyqs.h"
#include "simpqs.h"

#include <time.h>
#include <string.h>
#include <errno.h>

t_tm *taum = NULL;
uint taum_alloc = 0;
uint taum_size;
uint test_rough = 0;
mpz_t tmf, tmf2, tmp_lim;
#define SIMPQS_SIZE 66
mpz_t simpqs_array[SIMPQS_SIZE];

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
     _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

struct timespec cg_tp0;
struct timespec cg_tp1;
#define GIG 1000000000
static inline ulong cgdiff(struct timespec *t0) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp1);
    return (cg_tp1.tv_sec - t0->tv_sec) * GIG
            + cg_tp1.tv_nsec - t0->tv_nsec;
}
#ifdef VERBOSE
#define dz(format, ...) do { \
    gmp_printf("(%ld) " format "\n", cgdiff(&cg_tp0), ## __VA_ARGS__); \
} while (0)
#else
#define dz(...) 1
#endif

#ifdef VERBOSE
static inline bool ct_prime(mpz_t n) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = _GMP_is_prob_prime(n);
    gmp_printf("(%ld) p: %Zd %u\n", cgdiff(&cg_tp0), n, r ? 1 : 0);
    return r;
}
static inline ulong ct_power(mpz_t n) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    ulong r = power_factor(n, n);
    gmp_printf("(%ld) pow: %Zd^%lu\n", cgdiff(&cg_tp0), n, r);
    return r;
}
static inline bool ct_ecm(mpz_t n, mpz_t f, ulong b1, ulong curves) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = _GMP_ECM_FACTOR(n, f, b1, curves);
    gmp_printf("(%ld) ecm: %Zu [%lu, %lu] %u", cgdiff(&cg_tp0), n, b1, curves, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_pminus1(mpz_t n, mpz_t f, ulong b1, ulong b2) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = _GMP_pminus1_factor(n, f, b1, b2);
    gmp_printf("(%ld) p-1: %Zu [%lu, %lu] %u", cgdiff(&cg_tp0), n, b1, b2, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_tinyqs(mpz_t n, mpz_t f) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = tinyqs(n, f);
    gmp_printf("(%ld) tqs: %Zu %u", cgdiff(&cg_tp0), n, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_simpqs(mpz_t n, mpz_t *fa) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    int r = _GMP_simpqs(n, fa);
    gmp_printf("(%ld) sqs: %Zu %d", cgdiff(&cg_tp0), n, r);
    if (r)
        gmp_printf(" %Zu", fa[0]);
    gmp_printf("\n");
    return r;
}
static inline bool ct_holf(mpz_t n, mpz_t f, ulong rounds) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = _GMP_holf_factor(n, f, rounds);
    gmp_printf("(%ld) hlf: %Zu [%lu] %d", cgdiff(&cg_tp0), n, rounds, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_squfof(mpz_t n, mpz_t f, ulong rounds) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = squfof126(n, f, rounds);
    gmp_printf("(%ld) sqf: %Zu [%lu] %d", cgdiff(&cg_tp0), n, rounds, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_brent63(mpz_t n, mpz_t f, ulong rounds) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = pbrent63(n, f, rounds);
    gmp_printf("(%ld) b63: %Zu [%lu] %d", cgdiff(&cg_tp0), n, rounds, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
static inline bool ct_brent(mpz_t n, mpz_t f, ulong a, ulong rounds) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = _GMP_pbrent_factor(n, f, a, rounds);
    gmp_printf("(%ld) sqf: %Zu [%lu, %lu] %d", cgdiff(&cg_tp0), n, a, rounds, r);
    if (r)
        gmp_printf(" %Zu", f);
    gmp_printf("\n");
    return r;
}
extern int fs_trial(factor_state* fs);
static inline bool ct_trial(factor_state *fs) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
    bool r = fs_trial(fs);
    gmp_printf("(%ld) div: %Zu %d", cgdiff(&cg_tp0), fs->n, r);
    if (r)
        gmp_printf(" %Zu", fs->f);
    gmp_printf("\n");
    return r;
}
#else
#   define ct_prime(n) _GMP_is_prob_prime(n)
#   define ct_power(n) power_factor(n, n)
#   define ct_ecm(n, f, b1, curves) _GMP_ECM_FACTOR(n, f, b1, curves)
#   define ct_pminus1(n, f, b1, b2) _GMP_pminus1_factor(n, f, b1, b2)
#   define ct_tinyqs(n, f) tinyqs(n, f)
#   define ct_simpqs(n, fa) _GMP_simpqs(n, fa)
#   define ct_holf(n, f, rounds) _GMP_holf_factor(n, f, rounds)
#   define ct_squfof(n, f, rounds) squfof126(n, f, rounds)
#   define ct_brent63(n, f, rounds) pbrent63(n, f, rounds)
#   define ct_brent(n, f, a, rounds) _GMP_pbrent_factor(n, f, a, rounds)
#   define ct_trial(fs) fs_trial(fs)
#endif

#define NPRIMES_SMALL 2000
/* MPUG declares this static, so we must copy it */
static unsigned short primes_small[NPRIMES_SMALL];

void init_tmfbl(void);
void done_tmfbl(void);
void init_tau(uint rough) {
    UV pn;
    test_rough = rough;
    PRIME_ITERATOR(iter);
    primes_small[0] = 0;
    primes_small[1] = 2;
    for (pn = 2; pn < NPRIMES_SMALL; pn++)
        primes_small[pn] = prime_iterator_next(&iter);
    prime_iterator_destroy(&iter);
    mpz_init(tmf);
    mpz_init(tmf2);
    mpz_init(tmp_lim);
    init_tmfbl();
    for (uint i = 0; i < SIMPQS_SIZE; ++i)
        mpz_init(simpqs_array[i]);
}

void done_tau(void) {
    for (uint i = 0; i < SIMPQS_SIZE; ++i)
        mpz_clear(simpqs_array[i]);
    done_tmfbl();
    if (taum_alloc)
        for (uint i = 0; i < taum_alloc; ++i)
            mpz_clear(taum[i].n);
    free(taum);
    mpz_clear(tmf);
    mpz_clear(tmf2);
    mpz_clear(tmp_lim);
}

void alloc_taum(uint size) {
    if (size > taum_alloc) {
        taum = (t_tm *)realloc(taum, size * sizeof(t_tm));
        for (uint i = taum_alloc; i < size; ++i)
            mpz_init(taum[i].n);
        taum_alloc = size;
    }
}

void fs_init(factor_state* fs) {
    fs->state = FS_INIT;
    mpz_init(fs->n);
    mpz_init(fs->f);
    fs->ef = 1;
    fs->sp = 0;
    fs->log = get_verbose_level();
    fs->ntofac = 0;
    fs->tlim = 0;
    return;
}

void fs_clear(factor_state* fs) {
    mpz_clear(fs->n);
    mpz_clear(fs->f);
    for (int i = 0; i < fs->ntofac; ++i)
        mpz_clear(fs->tofac_stack[i]);
    return;
}

int fs_trial(factor_state* fs) {
    UV tlim = fs->tlim;
    UV sp = fs->sp;
    UV un = fs->un;
    UV lim, p;

    if (sp == 0) {
        int e2 = 0;
        while (mpz_even_p(fs->n)) {
            mpz_divexact_ui(fs->n, fs->n, 2);
            ++e2;
        }
        sp = 2;
        un = (mpz_cmp_ui(fs->n, 2 * tlim) >= 0) ? 2 * tlim : mpz_get_ui(fs->n);
        if (e2) {
            mpz_set_ui(fs->f, 2);
            fs->e = e2;
            goto found_trial;
        }
    }

    lim = (tlim < un) ? tlim : un;
    for (p = primes_small[sp]; p * p < lim; p = primes_small[++sp]) {
        int ep = 0;
        while (mpz_divisible_ui_p(fs->n, p)) {
            mpz_divexact_ui(fs->n, fs->n, p);
            ++ep;
        }
        if (ep) {
            mpz_set_ui(fs->f, p);
            fs->e = ep;
            un = (mpz_cmp_ui(fs->n, 2 * tlim) > 0)
                    ? 2 * tlim : mpz_get_ui(fs->n);
            goto found_trial;
        }
    }

    if (un < p * p) {
        mpz_set(fs->f, fs->n);
        fs->e = 1;
        mpz_set_ui(fs->n, 1);
        goto found_trial;
    }
    return 0;

found_trial:
    fs->sp = sp;
    fs->un = un;
    return 1;
}

static int _fs_remove(factor_state* fs) {
    int e = 0;
    e += mpz_remove(fs->n, fs->n, fs->f);
    for (int i = 0; i < fs->ntofac; ++i)
        e += mpz_remove(fs->tofac_stack[i], fs->tofac_stack[i], fs->f);
    return e;
}

/* Try to find one more prime factor and its exponent. Returns true if
 * it found one, else factorization is complete.
 */
int factor_one(factor_state* fs) {
    UV nbits = mpz_sizeinbase(fs->n, 2);
    UV B1;

fs_retry:
    if (mpz_cmp_ui(fs->n, 1) == 0) {
        if (fs->ntofac == 0) {
            fs->state = FS_TERM;
            return 0; /* no more factors to find */
        }
        --fs->ntofac;
        mpz_set(fs->n, fs->tofac_stack[fs->ntofac]);
        mpz_clear(fs->tofac_stack[fs->ntofac]);
        fs->state = FS_LARGE;
        goto fs_retry;
    }

    switch (fs->state) {
    default:
        croak("Unknown state\n");
    case FS_TERM:
        return 0;
    case FS_INIT:
        if (mpz_cmp_ui(fs->n, 0) == 0) {
            fs->state = FS_TERM;
            mpz_set_ui(fs->f, 0);
            fs->e = 1;
            return 1;
        }
        fs->sp = 0;
        fs->tlim = (nbits > 80) ? 4001 * 4001 : 16001 * 16001;
        fs->state = FS_TRIAL;
    case FS_TRIAL:
        if (ct_trial(fs))
            return 1;
        fs->state = FS_POWER;
    case FS_POWER:
        fs->ef = ct_power(fs->n);
        if (!fs->ef)
            fs->ef = 1;
        fs->state = FS_LARGE;
/*
 * This set of operations is meant to provide good performance for
 * "random" numbers as input. Hence we stack lots of effort up front
 * looking for small factors: prho and pbrent are ~ O(f^1/2) where
 * f is the smallest factor. SQUFOF is O(N^1/4), so arguably not
 * any better. p-1 and ECM are quite useful for pulling out small
 * factors (6-20 digits).
 *
 * Factoring a 778-digit number consisting of 101 8-digit factors
 * should complete in under 3 seconds. Factoring numbers consisting
 * of many 12-digit or 14-digit primes should take under 10 seconds.
 *
 * For state > 2, all returned exponents should be multiplied by fs->ef,
 * and tofac_stack must be checked (see label fs_retry above).
 */
    case FS_LARGE:
        if (mpz_cmp_ui(fs->n, fs->tlim) <= 0 || ct_prime(fs->n)) {
            mpz_set(fs->f, fs->n);
            fs->e = fs->ef * _fs_remove(fs);
            return 1;
        }

        if (nbits <= 63)
            if (ct_brent63(fs->n, fs->f, 400000))
                goto found_factor;
        if (nbits >= 65 && nbits <= 126) {
            if (ct_pminus1(fs->n, fs->f, 5000, 5000))
                goto found_factor;
            if (ct_tinyqs(fs->n, fs->f))
                goto found_factor;
        }
        /* It's possible the previous calls failed or weren't available */
        if (nbits <= 53)
            if (ct_squfof(fs->n, fs->f, 400000))
                goto found_factor;
        if (nbits <= 77) {
            int sb1 = (nbits < 58) ? 1
                : (nbits < 63) ? 2
                : (nbits < 72) ? 4
                : 10;
            if (ct_pminus1(fs->n, fs->f, sb1 * 1000, sb1 * 10000))
                goto found_factor;
            if (ct_squfof(fs->n, fs->f, 1000000))
                goto found_factor;
        }
        /* recheck power? */
        if (ct_pminus1(fs->n, fs->f, 20000, 200000))
            goto found_factor;

        /* Small ECM to find small factors */
        if (ct_ecm(fs->n, fs->f, 200, 4))
            goto found_factor;
        if (ct_ecm(fs->n, fs->f, 600, 20))
            goto found_factor;
        if (ct_ecm(fs->n, fs->f, 2000, 10))
            goto found_factor;

        /* Small p-1 */
        if (nbits < 100 || nbits >= 160)
            if (ct_pminus1(fs->n, fs->f, 200000, 3000000))
                goto found_factor;

        /* Set ECM parameters that have a good chance of success */
        UV curves;
        if (nbits < 100) {
            B1 = 5000;
            curves = 20;
        } else if (nbits < 128) {
            /* FIXME: surely this and the next case should be 'curves = 20'?
             * there was a comment on each "go to QS" - is the intent to do
             * a quick hit here, then rely on QS for more progress?
             */
            B1 = 10000;
            curves = 2;
        } else if (nbits < 160) {
            B1 = 20000;
            curves = 2;
        } else if (nbits < 192) {
            B1 = 30000;
            curves = 20;
        } else if (nbits < 224) {
            B1 = 40000;
            curves = 40;
        } else if (nbits < 256) {
            B1 = 80000;
            curves = 40;
        } else if (nbits < 512) {
            B1 = 160000;
            curves = 80;
        } else {
            B1 = 320000;
            curves = 160;
        }
        if (curves > 0 && ct_ecm(fs->n, fs->f, B1, curves))
            goto found_factor;

        /* QS (30+ digits). Fantastic if it is a semiprime, but can be
         * slow and a memory hog if not (compared to ECM). Restrict to
         * reasonable size numbers (< 91 digits). Because of the way it
         * works, it will generate (possibly) multiple factors for the same
         * amount of work. Go to some trouble to use them.
         */
        if (nbits >= 90 && nbits < 300) {
            mpz_t farray[66];
            int i, qs_nfactors;
            for (i = 0; i < 66; i++)
                mpz_init(farray[i]);
            qs_nfactors = ct_simpqs(fs->n, farray);
            mpz_set(fs->f, farray[0]);
            if (qs_nfactors > 2) {
                /* We found multiple factors */
                for (i = 2; i < qs_nfactors; i++) {
                    if (fs->log)
                        gmp_printf("SIMPQS found extra factor %Zd\n",
                                farray[i]);
                    if (fs->ntofac >= MAX_FACTORS - 1)
                        croak("Too many factors\n");
                    mpz_init_set(fs->tofac_stack[fs->ntofac], farray[i]);
                    ++fs->ntofac;
                    mpz_divexact(fs->n, fs->n, farray[i]);
                }
                /* f = farray[0], n = farray[1], farray[2..] pushed */
            }
            for (i = 0; i < 66; i++)
                mpz_clear(farray[i]);
            if (qs_nfactors > 1) {
                if (fs->log)
                    gmp_printf("SIMPQS found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }

        if (ct_ecm(fs->n, fs->f, 2 * B1, 20))
            goto found_factor;
        if (ct_brent(fs->n, fs->f, 1, 1024 * 1024))
            goto found_factor;
        if (ct_ecm(fs->n, fs->f, 4 * B1, 20))
            goto found_factor;
        if (ct_ecm(fs->n, fs->f, 8 * B1, 20))
            goto found_factor;
        /* HOLF in case it's a near-ratio-of-perfect-square */
        if (ct_holf(fs->n, fs->f, 1024 * 1024))
            goto found_factor;
        /* Large p-1 with stage 2: B2 = 20 * B1 */
        if (ct_pminus1(fs->n, fs->f, 5000000, 5000000 * 20))
            goto found_factor;
        if (ct_ecm(fs->n, fs->f, 32 * B1, 40))
            goto found_factor;
#if 0
        if (ct_brent(fs->n, fs->f, 2, 512 * 1024 * 1024))
            goto found_factor;
#endif

        /* Our method of last resort: ECM with high bmax and many curves*/
        if (fs->log)
            gmp_printf("starting large ECM on %Zd\n", fs->n);
        B1 *= 8;
        for (UV i = 0; i < 10; B1 *= 2, i++) {
            if (ct_ecm(fs->n, fs->f, B1, 100)) {
                if (!mpz_divisible_p(fs->n, fs->f)
                    || mpz_cmp_ui(fs->f, 1) == 0
                    || mpz_cmp(fs->f, fs->n) == 0
                ) {
                    gmp_printf("n = %Zd    f = %Zd\n", fs->n, fs->f);
                    croak("Incorrect factoring");
                }
                goto found_factor;
            }
        }

        /* FIXME: What to do with composites we can't factor?
         * Push them as "C#####" ? For now, just push them as if we factored.
         */
        if (fs->log)
            gmp_printf("gave up on %Zd\n", fs->n);
        goto found_factor;
    }
found_factor:
    if (!ct_prime(fs->f)) {
        mpz_init(fs->tofac_stack[fs->ntofac]);
        mpz_divexact(fs->tofac_stack[fs->ntofac], fs->n, fs->f);
        ++fs->ntofac;
        mpz_set(fs->n, fs->f);
        goto fs_retry;
    }
    fs->e = fs->ef * _fs_remove(fs);
    return 1;
}

/* Return true if what's left to factor is a prime */
static int _scanp(factor_state* fs) {
    int result = 0;
    if (mpz_cmp_ui(fs->n, 1) != 0) {
        if (ct_prime(fs->n))
            result = 1;
        else
            return 0;
    }
    for (int i = 0; i < fs->ntofac; ++i) {
        if (mpz_cmp_ui(fs->tofac_stack[i], 1) == 0)
            continue;
        if (result)
            return 0;
        if (ct_prime(fs->tofac_stack[i]))
            result = 1;
        else
            return 0;
    }
    return result;
}

/* Return true if tau(n^x) == k */
int is_taux(mpz_t n, uint32_t k, uint32_t x) {
    int cmp = mpz_cmp_ui(n, 1);
    int result = 0;
    factor_state fs;

#ifdef VERBOSE
    gmp_printf("is_taux t=%u (%u) %Zu^%u\n", k, mpz_sizeinbase(n, 2), n, x);
#endif
    if (cmp < 0)
        return 0;
    if (cmp == 0)
        return k == 1 || x == 0;
    if (k == 1 || x == 0)
        return 0;
    if (k == x + 1)
        return ct_prime(n) ? 1 : 0;

    fs_init(&fs);
    mpz_set(fs.n, n);
    while (1) {
        if ((k & 1) && (x & 1)) {
            int e = ct_power(fs.n);
            if (e == 0 || e & 1 || e > k)
                break;
            /* we actually need e divisible by gcd(map $_ - 1, divisors(k)) */
            x *= e;
        }
        if (k & 1) {
            while (1) {
                if (k == x + 1) {
                    result = _scanp(&fs);
                    break;
                }
                if (!factor_one(&fs))
                    break;
                if (k % (fs.e * x + 1))
                    break;
                k /= fs.e * x + 1;
                if (k == 1) {
                    result = (mpz_cmp_ui(fs.n, 1) == 0);
                    for (int i = 0; i < fs.ntofac; ++i)
                        result &= (mpz_cmp_ui(fs.tofac_stack[i], 1) == 0);
                    break;
                }
            }
            break;
        }
        if (k == x + 1) {
            result = _scanp(&fs);
            break;
        }
        if (!factor_one(&fs))
            break;
        if (k % (fs.e * x + 1))
            break;
        k /= fs.e * x + 1;
        if (k == 1) {
            result = (mpz_cmp_ui(fs.n, 1) == 0);
            for (int i = 0; i < fs.ntofac; ++i)
                result &= (mpz_cmp_ui(fs.tofac_stack[i], 1) == 0);
            break;
        }
        /* if fs.state == FS_LARGE, we need n >= max_trial_p ^ min_exp, where
         * max_trial_p = sqrt(fs.tlim), min_exp = sum(map $_ - 1, factor(k))
         */
    }
    fs_clear(&fs);
    return result;
}

/* We need to invoke simpqs, but find a single prime factor to return;
 * we accept that means we may do extra work to find back additional
 * factors.
 */
bool do_simpqs(mpz_t n, mpz_t f) {
    int qs = ct_simpqs(n, simpqs_array);

    /* if not factorized, it's a fail */
    if (qs < 2)
        return 0;

    /* look for a prime among the factors */
    for (uint i = 0; i < qs; ++i)
        if (ct_prime(simpqs_array[i])) {
            mpz_set(f, simpqs_array[i]);
            return 1;
        }

    /* all factors composite: pick the first one and find a factor there */
    factor_state fs;
    fs_init(&fs);
    mpz_set(fs.n, simpqs_array[0]);
    if (!factor_one(&fs)) {
        gmp_fprintf(stderr,
            "factor_one(%Zd) failed to find factor in simpqs composite\n",
            fs.n
        );
        exit(1);
    }
    mpz_set(f, fs.f);
    fs_clear(&fs);
    return 1;
}

static inline bool prep_abort(t_tm *tm, bool result) {
    if (result)
        tm->state = 0;  /* done */
    return result;
}

bool tau_multi_prep(uint i) {
    t_tm *tm = &taum[i];
    uint t = tm->t;
    uint nbits = mpz_sizeinbase(tm->n, 2);
    tm->state = 1;  /* init */

#ifdef VERBOSE
    gmp_printf("tau_multi_prep t=%u (%u) %Zu\n", t, nbits, tm->n);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cg_tp0);
#endif
    if (t == 1) {
        dz("div: t=1");
        return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);
    }

    /* do FS_TRIAL stage directly */
    int ep = 0;
    while (mpz_even_p(tm->n)) {
        mpz_divexact_ui(tm->n, tm->n, 2);
        ++ep;
    }
    if (ep) {
        if ((t % (ep + 1)) != 0) {
            dz("div: %u ~| t=%u", ep + 1, t);
            return 0;
        }
        t /= ep + 1;
        if (t == 1) {
            dz("div: t=1");
            return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);
        }
        nbits -= ep;
        ep = 0;
    }

    UV p;
    UV sp = 2;
    UV tlim = (nbits > 80) ? 4001 * 4001 : 16001 * 16001;
    if (test_rough && t >= test_rough) {
        uint roughness = divisors[t].sumpm;
        mpz_root(tmp_lim, tm->n, roughness);
        if (mpz_fits_uint_p(tmp_lim)) {
            ulong lim = mpz_get_ui(tmp_lim);
            tlim = lim * lim;
        }
        /* else what? */
    }

    UV un = mpz_cmp_ui(tm->n, 2 * tlim) >= 0
        ? 2 * tlim
        : mpz_get_ui(tm->n);
    UV lim = (tlim < un) ? tlim : un;
    PRIME_ITERATOR(iter);
    while (1) {
        p = prime_iterator_next(&iter);
        if (p * p >= lim)
            break;
        while (mpz_divisible_ui_p(tm->n, p)) {
            mpz_divexact_ui(tm->n, tm->n, p);
            ++ep;
        }
        if (ep) {
            if ((t % (ep + 1)) != 0) {
                dz("div: %u ~| t=%u", ep + 1, t);
                prime_iterator_destroy(&iter);
                return 0;
            }
            t /= ep + 1;
            if (t == 1) {
                dz("div: t=1");
                prime_iterator_destroy(&iter);
                return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);
            } else if (t == 2) {
                dz("div: t=2");
                prime_iterator_destroy(&iter);
                return prep_abort(tm, ct_prime(tm->n));
            } else if (mpz_cmp_ui(tm->n, 1) == 0) {
                dz("div: n=1");
                prime_iterator_destroy(&iter);
                return 0;
            }
            ep = 0;
            if (test_rough && t >= test_rough) {
                uint roughness = divisors[t].sumpm;
                mpz_root(tmp_lim, tm->n, roughness);
                if (mpz_fits_uint_p(tmp_lim)) {
                    ulong lim = mpz_get_ui(tmp_lim);
                    tlim = lim * lim;
                }
                /* else what? */
            }
            un = mpz_cmp_ui(tm->n, 2 * tlim) > 0
                ? 2 * tlim
                : mpz_get_ui(tm->n);
            lim = (tlim < un) ? tlim : un;
        }
    }
    prime_iterator_destroy(&iter);

    if (un < p * p) {
        dz("div: tail is prime");
        return prep_abort(tm, t == ((mpz_cmp_ui(tm->n, 1) == 0) ? 1 : 2));
    } else if (mpz_cmp_ui(tm->n, 1) == 0) {
        dz("div: n = 1");
        return prep_abort(tm, t == 1);
    } else if (t == 2) {
        dz("div: t == 2");
        return prep_abort(tm, ct_prime(tm->n));
    }
    if (test_rough && t >= test_rough) {
        mpz_ui_pow_ui(tmp_lim, p, divisors[t].sumpm);
        if (mpz_cmp(tm->n, tmp_lim) < 0) {
            dz("div: rough[%u]", divisors[t].sumpm);
            return 0;
        }
    }
    dz("div: done, t=%u", t);
    if (ct_prime(tm->n))
        return prep_abort(tm, t == 2);
    tm->e = ct_power(tm->n);
    if (!tm->e)
        tm->e = 1;
    if ((t & 1) && (tm->e & 1))
        return 0;
    if (!(t & 1) && !(tm->e & 1))
        return 0;
    tm->t = t;
    tm->tlim = tlim;
    return 1;
}

bool tmf_2(t_tm *tm) { return ct_brent63(tm->n, tmf, 400000); }
bool tmf_3(t_tm *tm) { return ct_pminus1(tm->n, tmf, 5000, 5000); }
bool tmf_4(t_tm *tm) { return ct_tinyqs(tm->n, tmf); }
bool tmf_5(t_tm *tm) { return ct_squfof(tm->n, tmf, 400000); }
bool tmf_6(t_tm *tm) { return ct_pminus1(tm->n, tmf, 1000, 10000); }
bool tmf_7(t_tm *tm) { return ct_pminus1(tm->n, tmf, 2000, 20000); }
bool tmf_8(t_tm *tm) { return ct_pminus1(tm->n, tmf, 4000, 40000); }
bool tmf_9(t_tm *tm) { return ct_pminus1(tm->n, tmf, 10000, 100000); }
bool tmf_10(t_tm *tm) { return ct_squfof(tm->n, tmf, 1000000); }
bool tmf_11(t_tm *tm) { return ct_pminus1(tm->n, tmf, 20000, 200000); }
bool tmf_12(t_tm *tm) { return ct_ecm(tm->n, tmf, 200, 4); }
bool tmf_13(t_tm *tm) { return ct_ecm(tm->n, tmf, 600, 20); }
bool tmf_14(t_tm *tm) { return ct_ecm(tm->n, tmf, 2000, 10); }
bool tmf_15(t_tm *tm) { return ct_pminus1(tm->n, tmf, 200000, 3000000); }
bool tmf_16(t_tm *tm) { tm->B1 = 5000; return ct_ecm(tm->n, tmf, tm->B1, 20); }
/* FIXME: surely this and the next case should have curves = 20?
 * There was a comment on each "go to QS" - is the intent to do
 * a quick hit here, then rely on QS for more progress?
 */
bool tmf_17(t_tm *tm) { tm->B1 = 10000; return ct_ecm(tm->n, tmf, tm->B1, 2); }
bool tmf_18(t_tm *tm) { tm->B1 = 20000; return ct_ecm(tm->n, tmf, tm->B1, 2); }
bool tmf_19(t_tm *tm) { tm->B1 = 30000; return ct_ecm(tm->n, tmf, tm->B1, 20); }
bool tmf_20(t_tm *tm) { tm->B1 = 40000; return ct_ecm(tm->n, tmf, tm->B1, 40); }
bool tmf_21(t_tm *tm) { tm->B1 = 80000; return ct_ecm(tm->n, tmf, tm->B1, 40); }
bool tmf_22(t_tm *tm) { tm->B1 = 160000; return ct_ecm(tm->n, tmf, tm->B1, 80); }
bool tmf_23(t_tm *tm) { tm->B1 = 320000; return ct_ecm(tm->n, tmf, tm->B1, 160); }
bool tmf_24(t_tm *tm) { return do_simpqs(tm->n, tmf); }
bool tmf_25(t_tm *tm) { return ct_ecm(tm->n, tmf, 2 * tm->B1, 20); }
bool tmf_26(t_tm *tm) { return ct_brent(tm->n, tmf, 1, 1 << 20); }
bool tmf_27(t_tm *tm) { return ct_ecm(tm->n, tmf, 4 * tm->B1, 20); }
bool tmf_28(t_tm *tm) { return ct_ecm(tm->n, tmf, 8 * tm->B1, 20); }
bool tmf_29(t_tm *tm) { return ct_holf(tm->n, tmf, 1 << 20); }
bool tmf_30(t_tm *tm) { return ct_pminus1(tm->n, tmf, 5000000, 5000000 * 20); }
bool tmf_31(t_tm *tm) { return ct_ecm(tm->n, tmf, 32 * tm->B1, 40); }
/* last resort tests */
bool tmf_32(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 4, 100); }
bool tmf_33(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 5, 100); }
bool tmf_34(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 6, 100); }
bool tmf_35(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 7, 100); }
bool tmf_36(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 8, 100); }
bool tmf_37(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 9, 100); }
bool tmf_38(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 10, 100); }
bool tmf_39(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 11, 100); }
bool tmf_40(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 12, 100); }
bool tmf_41(t_tm *tm) { return ct_ecm(tm->n, tmf, tm->B1 << 13, 100); }

typedef bool (*t_tmf)(t_tm *tm);
const t_tmf tmfa[] = {
    NULL, NULL, &tmf_2, &tmf_3, &tmf_4, &tmf_5, &tmf_6, &tmf_7,
    &tmf_8, &tmf_9, &tmf_10, &tmf_11, &tmf_12, &tmf_13, &tmf_14, &tmf_15,
    &tmf_16, &tmf_17, &tmf_18, &tmf_19, &tmf_20, &tmf_21, &tmf_22, &tmf_23,
    &tmf_24, &tmf_25, &tmf_26, &tmf_27, &tmf_28, &tmf_29, &tmf_30, &tmf_31,
    &tmf_32, &tmf_33, &tmf_34, &tmf_35, &tmf_36, &tmf_37, &tmf_38, &tmf_39,
    &tmf_40, &tmf_41
};
#define TM_TERM 0
/* leave room for possible power check at tmf_1() */
#define TM_INIT 2
#define TM_MAX (sizeof(tmfa) / sizeof(t_tmf))
typedef struct s_tmf_bits {
    uint maxlen;
    ulong tmf_bits;
} t_tmf_bits;

t_tmf_bits tmfb[] = {
    { 53, 0b11111110000000011111110000100100},
    { 58, 0b11111110000000011111110001000100},
    { 63, 0b11111110000000011111110010000100},
    { 64, 0b11111110000000011111110100000000},
    { 72, 0b11111110000000011111111100011000},
    { 77, 0b11111110000000011111111000011000},
    { 89, 0b11111110000000011111100000011000},
    { 99, 0b11111111000000011111100000011000},
    {126, 0b11111111000000100111100000011000},
    {127, 0b11111111000000100111100000000000},
    {159, 0b11111111000001000111100000000000},
    {191, 0b11111111000010001111100000000000},
    {223, 0b11111111000100001111100000000000},
    {255, 0b11111111001000001111100000000000},
    {299, 0b11111111010000001111100000000000},
    {511, 0b11111110010000001111100000000000},
    {  0, 0b11111110100000001111100000000000}
};
#define TMFB_MAX (sizeof(tmfb) / sizeof(t_tmf_bits))
ulong *tmfbl = NULL;
uint tmfb_maxb;
ulong tmfb_lim;

void init_tmfbl(void) {
    tmfb_lim = tmfb[TMFB_MAX - 1].tmf_bits;
    tmfb_maxb = tmfb[TMFB_MAX - 2].maxlen;
    tmfbl = (ulong *)malloc((tmfb_maxb + 1) * sizeof(ulong));
    uint i = 0;
    for (uint j = 0; j <= TMFB_MAX - 2; ++j) {
        uint lim = tmfb[j].maxlen;
        ulong bits = tmfb[j].tmf_bits;
        for (; i <= lim; ++i)
            tmfbl[i] = bits;
    }
}

void done_tmfbl(void) {
    free(tmfbl);
}

static inline ulong _find_tmfb(uint size) {
    if (size <= tmfb_maxb)
        return tmfbl[size];
    return tmfb_lim;
}

/* see also other_comparator() in coul.c */
int taum_comparator(const void *va, const void *vb) {
    t_tm *tma = (t_tm *)va;
    t_tm *tmb = (t_tm *)vb;
    uint at2 = tma->t ^ (tma->t - 1);
    uint bt2 = tmb->t ^ (tmb->t - 1);
    if (at2 < bt2)
        return -1;
    if (at2 > bt2)
        return 1;
    if (tma->t < tmb->t)
        return -1;
    if (tma->t > tmb->t)
        return 1;
    return mpz_cmp(tma->n, tmb->n);
}

mpz_t *tm_factor(t_tm *tm) {
    if (ct_prime(tmf))
        return &tmf;
    /* we have a composite factor */
    mpz_divexact(tmf2, tm->n, tmf);
    if (ct_prime(tmf2))
        return &tmf2;
    /* .. leaving a composite residue */
    factor_state fs;
    fs_init(&fs);
    fs.tlim = tm->tlim;
    if (mpz_cmp(tmf, tmf2) < 0)
        mpz_set(fs.n, tmf);
    else
        mpz_set(fs.n, tmf2);
    fs.state = FS_POWER;
    if (!factor_one(&fs)) {
        gmp_fprintf(stderr, "no factor found for composite %Zd\n", fs.n);
        exit(1);
    }
    mpz_set(tmf, fs.f);
    fs_clear(&fs);
    return &tmf;
}

bool tau_multi_run(uint count) {
    uint i = 0;
    /* Shuffle the entries that did not complete by trial division to
     * the front. Find size and thus the associated tmfb entry for each. */
    for (uint j = 0; j < count; ++j) {
        if (taum[j].state == 0)
            continue;
        if (i < j) {
            mpz_swap(taum[i].n, taum[j].n);
            taum[i].t = taum[j].t;
            taum[i].e = taum[j].e;
        }
        taum[i].state = 2;
        taum[i].bits = _find_tmfb(mpz_sizeinbase(taum[i].n, 2));
        ++i;
    }
    if (i == 0)
        return 1;
    count = i;

    qsort(taum, count, sizeof(t_tm), &taum_comparator);

    uint next_i;
  tmr_retry:
    for (uint i = TM_INIT; i < TM_MAX; i = next_i) {
        next_i = i + 1;
        for (uint j = 0; j < count; ++j) {
            t_tm *tm = &taum[j];
          tmr_redo:
            if (tm->state > i)
                continue;
            if (!(tm->bits & (1UL << i)))
                continue;
            if (!(*tmfa[i])(tm)) {
                tm->state = i + 1;
                continue;
            }

            /* make sure we have a _prime_ factor */
            mpz_t *f = tm_factor(tm);
            uint e = 0;
            while (mpz_divisible_p(tm->n, *f)) {
                ++e;
                mpz_divexact(tm->n, tm->n, *f);
            }
            if (e == 0) {
                gmp_fprintf(stderr,
                    "state %d found non-divisible factor %Zd for %Zd\n",
                    tm->state, *f, tm->n
                );
                exit(1);
            }
            e = e * tm->e + 1;
            if (tm->t % e)
                return 0;
            tm->t /= e;
            if (tm->t == 1) {
                if (mpz_cmp_ui(tm->n, 1) != 0)
                    return 0;
                goto tmr_splice;
            } else if (tm->t == 2) {
                if (!ct_prime(tm->n))
                    return 0;
                goto tmr_splice;
            } else if (mpz_cmp_ui(tm->n, 1) == 0)
                return 0;
            else if (tm->t & 1) {
                /* odd tau should be easy, do immediate full check */
                if (!is_taux(tm->n, 1, tm->t))
                    return 0;
                goto tmr_splice;
            } else if (ct_prime(tm->n))
                return 0;
            else if ((tm->t & 1) && (tm->e & 1)) {
                e = ct_power(tm->n);
                if (e == 0 || e & 1 || e > tm->t)
                    return 0;
                tm->e *= e;
            }
            tm->state = TM_INIT;
            next_i = TM_INIT;
            tm->bits = _find_tmfb(mpz_sizeinbase(tm->n, 2));
            goto tmr_retry;

          tmr_splice:
            --count;
            if (count == 0)
                return 1;
            if (j < count) {
                mpz_swap(taum[j].n, taum[count].n);
                taum[j].t = taum[count].t;
                taum[j].state = taum[count].state;
                taum[j].e = taum[count].e;
                taum[j].bits = taum[count].bits;
                taum[j].B1 = taum[count].B1;
                taum[j].tlim = taum[count].tlim;
            }
            goto tmr_redo;
        }
    }
    gmp_fprintf(stderr,
        "no factors found for %u non-primes starting %Zd (%u)\n",
        count, taum[0].n, taum[0].t
    );
    exit(1);
}

bool tau_single_try(uint i) {
    t_tm *tm = &taum[i];
    for (uint i = TM_INIT; i < TM_MAX; ++i) {
        if (!(tm->bits & (1UL << i)))
            continue;
        if (!(*tmfa[i])(tm)) {
            tm->state = i + 1;
            continue;
        }

        mpz_set(tm->n, tmf);
        return 1;
    }
    return 0;
}
