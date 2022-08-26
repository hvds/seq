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

t_tm *taum = NULL;
uint taum_alloc = 0;
uint taum_size;
mpz_t tmf;
#define SIMPQS_SIZE 66
mpz_t simpqs_array[SIMPQS_SIZE];

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
     _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

#define NPRIMES_SMALL 2000
/* MPUG declares this static, so we must copy it */
static unsigned short primes_small[NPRIMES_SMALL];

void init_tmfbl(void);
void done_tmfbl(void);
void init_tau(void) {
    UV pn;
    PRIME_ITERATOR(iter);
    primes_small[0] = 0;
    primes_small[1] = 2;
    for (pn = 2; pn < NPRIMES_SMALL; pn++)
        primes_small[pn] = prime_iterator_next(&iter);
    prime_iterator_destroy(&iter);
    mpz_init(tmf);
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
        if (mpz_cmp_ui(fs->n, fs->tlim) <= 0 || _GMP_is_prob_prime(fs->n)) {
            mpz_set(fs->f, fs->n);
            fs->e = fs->ef * _fs_remove(fs);
            return 1;
        }

        if (nbits <= 63) {
            if (pbrent63(fs->n, fs->f, 400000)) {
                if (fs->log)
                    gmp_printf("UV Rho-Brent found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        if (nbits >= 65 && nbits <= 126) {
            if (_GMP_pminus1_factor(fs->n, fs->f, 5000, 5000)) {
                if (fs->log)
                    gmp_printf("p-1 (%dk) found factor %Zd\n", 5000, fs->f);
                goto found_factor;
            }
            if (tinyqs(fs->n, fs->f)) {
                if (fs->log)
                    gmp_printf("tinyqs found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        /* It's possible the previous calls failed or weren't available */
        if (nbits <= 53) {
            if (squfof126(fs->n, fs->f, 400000)) {
                if (fs->log)
                    gmp_printf("UV SQUFOF126 found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        if (nbits <= 77) {
            int sb1 = (nbits < 58) ? 1
                : (nbits < 63) ? 2
                : (nbits < 72) ? 4
                : 10;
            if (_GMP_pminus1_factor(fs->n, fs->f, sb1 * 1000, sb1 * 10000)) {
                if (fs->log)
                    gmp_printf("p-1 (%dk) found factor %Zd\n", sb1, fs->f);
                goto found_factor;
            }
            if (squfof126(fs->n, fs->f, 1000000)) {
                if (fs->log)
                    gmp_printf("SQUFOF126 found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }
        /* recheck power? */
        if (_GMP_pminus1_factor(fs->n, fs->f, 20000, 200000)) {
            if (fs->log)
                gmp_printf("p-1 (20k) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Small ECM to find small factors */
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 200, 4)) {
            if (fs->log)
                gmp_printf("tiny ecm (200) found factor %Zd\n", fs->f);
            goto found_factor;
        }
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 600, 20)) {
            if (fs->log)
                gmp_printf("tiny ecm (600) found factor %Zd\n", fs->f);
            goto found_factor;
        }
        if (_GMP_ECM_FACTOR(fs->n, fs->f, 2000, 10)) {
            if (fs->log)
                gmp_printf("tiny ecm (2000) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Small p-1 */
        if (nbits < 100 || nbits >= 160) {
            if (_GMP_pminus1_factor(fs->n, fs->f, 200000, 3000000)) {
                if (fs->log)
                    gmp_printf("p-1 (200k) found factor %Zd\n", fs->f);
                goto found_factor;
            }
        }

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
        if (curves > 0 && _GMP_ECM_FACTOR(fs->n, fs->f, B1, curves)) {
            if (fs->log)
                gmp_printf("small ecm (%luk,%lu) found factor %Zd\n",
                        B1 / 1000, curves, fs->f);
            goto found_factor;
        }

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
            qs_nfactors = _GMP_simpqs(fs->n, farray);
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

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 2 * B1, 20)) {
            if (fs->log)
                gmp_printf("ecm (%luk,20) found factor %Zd\n",
                        2 * B1 / 1000, fs->f);
            goto found_factor;
        }

        if (_GMP_pbrent_factor(fs->n, fs->f, 1, 1024 * 1024)) {
            if (fs->log)
                gmp_printf("pbrent (1,1M) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 4 * B1, 20)) {
            if (fs->log)
                gmp_printf("ecm (%luk,20) ecm found factor %Zd\n",
                        4 * B1, fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 8 * B1, 20)) {
            if (fs->log)
                gmp_printf("ecm (%luk,20) ecm found factor %Zd\n",
                        8 * B1, fs->f);
            goto found_factor;
        }

        /* HOLF in case it's a near-ratio-of-perfect-square */
        if (_GMP_holf_factor(fs->n, fs->f, 1024 * 1024)) {
            if (fs->log)
                gmp_printf("holf found factor %Zd\n", fs->f);
            goto found_factor;
        }

        /* Large p-1 with stage 2: B2 = 20 * B1 */
        if (_GMP_pminus1_factor(fs->n, fs->f, 5000000, 5000000 * 20)) {
            if (fs->log)
                gmp_printf("p-1 (5M) found factor %Zd\n", fs->f);
            goto found_factor;
        }

        if (_GMP_ECM_FACTOR(fs->n, fs->f, 32 * B1, 40)) {
            if (fs->log)
                gmp_printf("ecm (%luk,40) ecm found factor %Zd\n",
                        32 * B1, fs->f);
            goto found_factor;
        }

#if 0
        if (_GMP_pbrent_factor(fs->n, fs->f, 2, 512 * 1024 * 1024)) {
            if (fs->log)
                gmp_printf("pbrent (2,512M) found factor %Zd\n", fs->f);
            goto found_factor;
        }
#endif

        /* Our method of last resort: ECM with high bmax and many curves*/
        if (fs->log)
            gmp_printf("starting large ECM on %Zd\n", fs->n);
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
                if (fs->log)
                    gmp_printf("ecm (%luk,100) ecm found factor %Zd\n",
                            B1, fs->f);
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
static int _scanp(factor_state* fs) {
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
        if (result)
            return 0;
        if (_GMP_is_prime(fs->tofac_stack[i]))
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

    if (cmp < 0)
        return 0;
    if (cmp == 0)
        return k == 1 || x == 0;
    if (k == 1 || x == 0)
        return 0;
    if (k == x + 1)
        return _GMP_is_prime(n) ? 1 : 0;

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

/* We need to invoke simpqs, but find a single prime factor to return;
 * we accept that means we may do extra work to find back additional
 * factors.
 */
bool do_GMP_simpqs(mpz_t n, mpz_t f) {
    int qs = _GMP_simpqs(n, simpqs_array);

    /* if not factorized, it's a fail */
    if (qs < 2)
        return 0;

    /* look for a prime among the factors */
    for (uint i = 0; i < qs; ++i)
        if (_GMP_is_prob_prime(simpqs_array[i])) {
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

    if (t == 1)
        return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);

    /* do FS_TRIAL stage directly */
    int ep = 0;
    while (mpz_even_p(tm->n)) {
        mpz_divexact_ui(tm->n, tm->n, 2);
        ++ep;
    }
    if (ep) {
        if ((t % (ep + 1)) != 0)
            return 0;
        t /= ep + 1;
        if (t == 1)
            return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);
        nbits -= ep;
        ep = 0;
    }

    UV p;
    UV sp = 2;
    UV tlim = (nbits > 80) ? 4001 * 4001 : 16001 * 16001;
    UV un = mpz_cmp_ui(tm->n, 2 * tlim) >= 0
        ? 2 * tlim
        : mpz_get_ui(tm->n);
    UV lim = (tlim < un) ? tlim : un;
    for (p = primes_small[sp]; p * p < lim; p = primes_small[++sp]) {
        while (mpz_divisible_ui_p(tm->n, p)) {
            mpz_divexact_ui(tm->n, tm->n, p);
            ++ep;
        }
        if (ep) {
            if ((t % (ep + 1)) != 0)
                return 0;
            t /= ep + 1;
            if (t == 1)
                return prep_abort(tm, mpz_cmp_ui(tm->n, 1) == 0);
            if (t == 2)
                return prep_abort(tm, _GMP_is_prob_prime(tm->n));
            if (mpz_cmp_ui(tm->n, 1) == 0)
                return 0;
            ep = 0;
            un = mpz_cmp_ui(tm->n, 2 * tlim) > 0
                ? 2 * tlim
                : mpz_get_ui(tm->n);
            lim = (tlim < un) ? tlim : un;
        }
    }

    if (un < p * p)
        return prep_abort(tm, t == (mpz_cmp_ui(tm->n, 1) == 0 ? 1 : 2));
    if (mpz_cmp_ui(tm->n, 1) == 0)
        return prep_abort(tm, t == 1);
    if (t == 2)
        return prep_abort(tm, _GMP_is_prob_prime(tm->n));
    if (_GMP_is_prob_prime(tm->n))
        return prep_abort(tm, t == 2);
    tm->e = power_factor(tm->n, tm->n);
    if (!tm->e)
        tm->e = 1;
    if ((t & 1) && (tm->e & 1))
        return 0;
    if (!(t & 1) && !(tm->e & 1))
        return 0;
    tm->t = t;
    return 1;
}

bool tmf_2(t_tm *tm) { return pbrent63(tm->n, tmf, 400000); }
bool tmf_3(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 5000, 5000); }
bool tmf_4(t_tm *tm) { return tinyqs(tm->n, tmf); }
bool tmf_5(t_tm *tm) { return squfof126(tm->n, tmf, 400000); }
bool tmf_6(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 1000, 10000); }
bool tmf_7(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 2000, 20000); }
bool tmf_8(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 4000, 40000); }
bool tmf_9(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 10000, 100000); }
bool tmf_10(t_tm *tm) { return squfof126(tm->n, tmf, 1000000); }
bool tmf_11(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 20000, 200000); }
bool tmf_12(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 200, 4); }
bool tmf_13(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 600, 20); }
bool tmf_14(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 2000, 10); }
bool tmf_15(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 200000, 3000000); }
bool tmf_16(t_tm *tm) {
    tm->B1 = 5000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 20);
}
/* FIXME: surely this and the next case should have curves = 20?
 * There was a comment on each "go to QS" - is the intent to do
 * a quick hit here, then rely on QS for more progress?
 */
bool tmf_17(t_tm *tm) {
    tm->B1 = 10000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 2);
}
bool tmf_18(t_tm *tm) {
    tm->B1 = 20000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 2);
}
bool tmf_19(t_tm *tm) {
    tm->B1 = 30000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 20);
}
bool tmf_20(t_tm *tm) {
    tm->B1 = 40000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 40);
}
bool tmf_21(t_tm *tm) {
    tm->B1 = 80000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 40);
}
bool tmf_22(t_tm *tm) {
    tm->B1 = 160000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 80);
}
bool tmf_23(t_tm *tm) {
    tm->B1 = 320000;
    return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1, 160);
}
bool tmf_24(t_tm *tm) { return do_GMP_simpqs(tm->n, tmf); }
bool tmf_25(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 2 * tm->B1, 20); }
bool tmf_26(t_tm *tm) { return _GMP_pbrent_factor(tm->n, tmf, 1, 1 << 20); }
bool tmf_27(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 4 * tm->B1, 20); }
bool tmf_28(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 8 * tm->B1, 20); }
bool tmf_29(t_tm *tm) { return _GMP_holf_factor(tm->n, tmf, 1 << 20); }
bool tmf_30(t_tm *tm) { return _GMP_pminus1_factor(tm->n, tmf, 5000000, 5000000 * 20); }
bool tmf_31(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, 32 * tm->B1, 40); }
/* last resort tests */
bool tmf_32(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 4, 100); }
bool tmf_33(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 5, 100); }
bool tmf_34(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 6, 100); }
bool tmf_35(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 7, 100); }
bool tmf_36(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 8, 100); }
bool tmf_37(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 9, 100); }
bool tmf_38(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 10, 100); }
bool tmf_39(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 11, 100); }
bool tmf_40(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 12, 100); }
bool tmf_41(t_tm *tm) { return _GMP_ECM_FACTOR(tm->n, tmf, tm->B1 << 13, 100); }

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
            if (!(tm->bits & (1 << i)))
                continue;
            if (!(*tmfa[i])(tm)) {
                tm->state = i + 1;
                continue;
            }
            uint e = 0;
            while (mpz_divisible_p(tm->n, tmf)) {
                ++e;
                mpz_divexact(tm->n, tm->n, tmf);
            }
            if (e == 0) {
                gmp_fprintf(stderr,
                    "state %d found non-divisible factor %Zd for %Zd\n",
                    tm->state, tmf, tm->n
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
                if (!_GMP_is_prob_prime(tm->n))
                    return 0;
                goto tmr_splice;
            } else if (mpz_cmp_ui(tm->n, 1) == 0)
                return 0;
            else if (tm->t & 1) {
                /* odd tau should be easy, do immediate full check */
                if (!is_taux(tm->n, 1, tm->t))
                    return 0;
                goto tmr_splice;
            } else if (_GMP_is_prob_prime(tm->n))
                return 0;
            else if ((tm->t & 1) && (tm->e & 1)) {
                e = power_factor(tm->n, tm->n);
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
                taum[j].e = taum[count].e;
                taum[j].state = taum[count].state;
                taum[j].bits = taum[count].bits;
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
