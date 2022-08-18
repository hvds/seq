#include <gmp.h>
#include "ptypes.h"

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

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
     _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

#define NPRIMES_SMALL 2000
/* MPUG declares this static, so we must copy it */
static unsigned short primes_small[NPRIMES_SMALL];
void init_tau(void) {
    UV pn;
    PRIME_ITERATOR(iter);
    primes_small[0] = 0;
    primes_small[1] = 2;
    for (pn = 2; pn < NPRIMES_SMALL; pn++) {
        primes_small[pn] = prime_iterator_next(&iter);
    }
    prime_iterator_destroy(&iter);
}

void fs_init(factor_state* fs)
{
    fs->state = FS_INIT;
    mpz_init(fs->n);
    mpz_init(fs->f);
    fs->ef = 1;
    fs->sp = 0;
    fs->log = get_verbose_level();
    fs->ntofac = 0;
    return;
}

void fs_clear(factor_state* fs)
{
    mpz_clear(fs->n);
    mpz_clear(fs->f);
    for (int i = 0; i < fs->ntofac; ++i)
        mpz_clear(fs->tofac_stack[i]);
    return;
}

int fs_trial(factor_state* fs)
{
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
            un = (mpz_cmp_ui(fs->n, 2 * tlim) > 0) ? 2 * tlim : mpz_get_ui(fs->n);
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

static int _fs_remove(factor_state* fs)
{
    int e = 0;
    e += mpz_remove(fs->n, fs->n, fs->f);
    for (int i = 0; i < fs->ntofac; ++i)
        e += mpz_remove(fs->tofac_stack[i], fs->tofac_stack[i], fs->f);
    return e;
}

/* Try to find one more prime factor and its exponent. Returns true if
 * it found one, else factorization is complete.
 */
int factor_one(factor_state* fs)
{
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
        if (fs_trial(fs))
            return 1;
        fs->state = FS_POWER;
    case FS_POWER:
        fs->ef = power_factor(fs->n, fs->f);
        if (fs->ef) {
            mpz_set(fs->n, fs->f);
        } else {
            fs->ef = 1;
        }
        fs->state = FS_LARGE;
/*
 * This set of operations is meant to provide good performance for
 * "random" numbers as input.    Hence we stack lots of effort up front
 * looking for small factors: prho and pbrent are ~ O(f^1/2) where
 * f is the smallest factor.    SQUFOF is O(N^1/4), so arguably not
 * any better.    p-1 and ECM are quite useful for pulling out small
 * factors (6-20 digits).
 *
 * Factoring a 778-digit number consisting of 101 8-digit factors
 * should complete in under 3 seconds.    Factoring numbers consisting
 * of many 12-digit or 14-digit primes should take under 10 seconds.
 *
 * For state > 2, all returned exponents should be multiplied by fs->ef,
 * and tofac_stack must be checked (see label fs_retry above).
 */
    case FS_LARGE:
        if (mpz_cmp_ui(fs->n, fs->tlim) <= 0 || _GMP_is_prob_prime(fs->n)) {
            mpz_set(fs->f, fs->n);
            fs->e = fs->ef * _fs_remove(fs);
            return 1;
        }

        if (nbits <= 63) {
            if (pbrent63(fs->n, fs->f, 400000)) {
                if (fs->log) gmp_printf("UV Rho-Brent found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        if (nbits >= 65 && nbits <= 126) {
            if (_GMP_pminus1_factor(fs->n, fs->f, 5000, 5000)) {
                if (fs->log) gmp_printf("p-1 (%dk) found factor %Zd\n", 5000, fs->f);
                goto found_factor;
            }
            if (tinyqs(fs->n, fs->f)) {
                if (fs->log) gmp_printf("tinyqs found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        /* It's possible the previous calls failed or weren't available */
        if (nbits <= 53) {
            if (squfof126(fs->n, fs->f, 400000)) {
                if (fs->log) gmp_printf("UV SQUFOF126 found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        if (nbits <= 77) {
            int sb1 = (nbits < 58) ?    1
                            : (nbits < 63) ?    2
                            : (nbits < 72) ?    4
                                                                : 10;
            if (_GMP_pminus1_factor(fs->n, fs->f, sb1*1000, sb1*10000)) {
                if (fs->log) gmp_printf("p-1 (%dk) found factor %Zd\n", sb1, fs->f);
                goto found_factor;
            }
            if (squfof126(fs->n, fs->f, 1000000)) {
                if (fs->log) gmp_printf("SQUFOF126 found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        /* recheck power? */
        if (_GMP_pminus1_factor(fs->n, fs->f, 20000, 200000)) {
            if (fs->log) gmp_printf("p-1 (20k) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Small ECM to find small factors */
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 200, 4)) {
            if (fs->log) gmp_printf("tiny ecm (200) found factor %Zd\n", fs->f);
            goto found_factor;
        }
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 600, 20)) {
            if (fs->log) gmp_printf("tiny ecm (600) found factor %Zd\n", fs->f);
            goto found_factor;
        }
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 2000, 10)) {
            if (fs->log) gmp_printf("tiny ecm (2000) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Small p-1 */
        if (nbits < 100 || nbits >= 160) {
            if (_GMP_pminus1_factor(fs->n, fs->f, 200000, 3000000)) {
                if (fs->log) gmp_printf("p-1 (200k) found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }

        /* Set ECM parameters that have a good chance of success */
        {
            UV curves;
            if            (nbits < 100){ B1 =     5000; curves =    20; }
            else if (nbits < 128){ B1 =    10000; curves =     2; } /* go to QS */
            else if (nbits < 160){ B1 =    20000; curves =     2; } /* go to QS */
            else if (nbits < 192){ B1 =    30000; curves =    20; }
            else if (nbits < 224){ B1 =    40000; curves =    40; }
            else if (nbits < 256){ B1 =    80000; curves =    40; }
            else if (nbits < 512){ B1 = 160000; curves =    80; }
            else                                 { B1 = 320000; curves = 160; }
            if (curves > 0 && _GMP_ECM_FACTOR(fs->n, fs->f, B1, curves)) {
                if (fs->log) gmp_printf("small ecm (%luk,%lu) found factor %Zd\n", B1/1000, curves, fs->f);
                goto found_factor;
            }
        }

        /* QS (30+ digits).    Fantastic if it is a semiprime, but can be
         * slow and a memory hog if not (compared to ECM).    Restrict to
         * reasonable size numbers (< 91 digits).    Because of the way it
         * works, it will generate (possibly) multiple factors for the same
         * amount of work.    Go to some trouble to use them. */
        if (nbits >= 90 && nbits < 300) {
            mpz_t farray[66];
            int i, qs_nfactors;
            for (i = 0; i < 66; i++)
                mpz_init(farray[i]);
            qs_nfactors = _GMP_simpqs(fs->n, farray);
            mpz_set(fs->f, farray[0]);
            if (qs_nfactors > 2) {
                /* We found multiple factors */
                for (i = 2; i < qs_nfactors; i++) {
                    if (fs->log) gmp_printf("SIMPQS found extra factor %Zd\n", farray[i]);
                    if (fs->ntofac >= MAX_FACTORS-1) croak("Too many factors\n");
                    mpz_init_set(fs->tofac_stack[fs->ntofac], farray[i]);
                    ++fs->ntofac;
                    mpz_divexact(fs->n, fs->n, farray[i]);
                }
                /* f = farray[0], n = farray[1], farray[2..] pushed */
            }
            for (i = 0; i < 66; i++)
                mpz_clear(farray[i]);
            if (qs_nfactors > 1) {
                if (fs->log) gmp_printf("SIMPQS found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 2*B1, 20)) {
            if (fs->log) gmp_printf("ecm (%luk,20) found factor %Zd\n", 2*B1/1000, fs->f);
            goto found_factor;
        }

        if (_GMP_pbrent_factor(fs->n, fs->f, 1, 1*1024*1024)) {
            if (fs->log) gmp_printf("pbrent (1,1M) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 4*B1, 20)) {
            if (fs->log) gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 4*B1, fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 8*B1, 20)) {
            if (fs->log) gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 8*B1, fs->f);
            goto found_factor;
        }

        /* HOLF in case it's a near-ratio-of-perfect-square */
        if (_GMP_holf_factor(fs->n, fs->f, 1*1024*1024)) {
            if (fs->log) gmp_printf("holf found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Large p-1 with stage 2: B2 = 20*B1 */
        if (_GMP_pminus1_factor(fs->n, fs->f, 5000000, 5000000*20)) {
            if (fs->log) gmp_printf("p-1 (5M) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 32*B1, 40)) {
            if (fs->log) gmp_printf("ecm (%luk,40) ecm found factor %Zd\n", 32*B1, fs->f);
            goto found_factor;
        }

        /*
        if (_GMP_pbrent_factor(fs->n, fs->f, 2, 512*1024*1024)) {
            if (fs->log) gmp_printf("pbrent (2,512M) found factor %Zd\n", fs->f);
            goto found_factor;
        }
        */

        /* Our method of last resort: ECM with high bmax and many curves*/
        if (fs->log) gmp_printf("starting large ECM on %Zd\n", fs->n);
        B1 *= 8;
        for (UV i = 0; i < 10; B1 *= 2, i++) {
            if (_GMP_ECM_FACTOR(fs->n, fs->f, B1, 100)) {
                if (!mpz_divisible_p(fs->n, fs->f)
                        || mpz_cmp_ui(fs->f, 1) == 0
                        || mpz_cmp(fs->f, fs->n) == 0
                ) {
                    gmp_printf("n = %Zd    f = %Zd\n", fs->n, fs->f);
                    croak("Incorrect factoring");
                }
                if (fs->log) gmp_printf("ecm (%luk,100) ecm found factor %Zd\n", B1, fs->f);
                goto found_factor;
            }
        }

        /* TODO: What to do with composites we can't factor?
         *             Push them as "C#####" ?
         *             For now, just push them as if we factored.
         */
        if (fs->log) gmp_printf("gave up on %Zd\n", fs->n);
        goto found_factor;
    }
found_factor:
    if (!_GMP_is_prob_prime(fs->f)) {
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
static int _scanp(factor_state* fs)
{
    int result = 0;
    if (mpz_cmp_ui(fs->n, 1) != 0) {
        if (_GMP_is_prime(fs->n))
            result = 1;
        else
            return 0;
    }
    for (int i = 0; i < fs->ntofac; ++i) {
        if (mpz_cmp_ui(fs->tofac_stack[i], 1) == 0)
            continue;
        if (result) return 0;
        if (_GMP_is_prime(fs->tofac_stack[i])) {
            result = 1;
        } else
            return 0;
    }
    return result;
}

/* Return true if tau(n^x) == k */
int is_taux(mpz_t n, uint32_t k, uint32_t x) {
    int cmp = mpz_cmp_ui(n, 1);
    int result = 0;
    factor_state fs;

    if (cmp < 0) return 0;
    if (cmp == 0) return k == 1 || x == 0;
    if (k == 1 || x == 0) return 0;
    if (k == x + 1) return _GMP_is_prime(n) ? 1 : 0;

    fs_init(&fs);
    mpz_set(fs.n, n);
    while (1) {
        if ((k & 1) && (x & 1)) {
            int e = power_factor(fs.n, fs.f);
            if (e == 0 || e & 1 || e > k)
                break;
            mpz_set(fs.n, fs.f);
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
