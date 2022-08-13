#include <stdlib.h>
#include <string.h>
#include "rootmod.h"
#include "coulfact.h"
#include "utility.h"

extern void fs_init(factor_state *fs);
extern void fs_clear(factor_state *fs);

typedef struct s_results {
    uint size;
    uint count;
    mpz_t *r;
} t_results;

typedef enum {
    rm_base,
    arm_scratch,
    armkp_base, armkp_new, armkp_stash,
    armc_r1,
    armpp_copy,
    armppr_copy,
    asmf_r1, asmf_r2,

    E_RESULTS_MAX
} e_results;

t_results ra[E_RESULTS_MAX];

typedef struct s_lpow {
    ulong p;
    uint e;
} t_lpow;
t_lpow *rm_nf = NULL;
uint rm_nf_size = 0;
t_ppow *rm_kf = NULL;
uint rm_kf_size = 0;

typedef enum {
    rm_r, rm_p,
    arm_a,
    armkp_m, armkp_px,
    armpp_px,
    armc_n, armc_inv,
    arzpp_px,
    armppr_a, armppr_z, armppr_px, armppr_px2,
    earmpp_t, earmpp_t1, earmpp_t2, earmpp_g,
    armp_a, armp_t,
    tsp_A, tsp_B, tsp_T, tsp_y, tsp_z,
    si_d, si_m,

    E_RMSTASH_MAX
} e_rmstash;
mpz_t rm_stash[E_RMSTASH_MAX];
static inline mpz_t *ZP(e_rmstash e) { return &rm_stash[e]; }
#define Z(e) *ZP(e)

void resize_results(t_results *rp, uint size) {
    if (rp->size < size) {
        if (size < rp->size + 16)
            size = rp->size + 16;
        rp->r = (mpz_t *)realloc(rp->r, size * sizeof(mpz_t));
        for (uint i = rp->size; i < size; ++i)
            mpz_init(rp->r[i]);
        rp->size = size;
    }
}

void resize_nf(uint size) {
    if (size > rm_nf_size) {
        if (size < rm_nf_size + 16)
            size = rm_nf_size + 16;
        rm_nf = (t_lpow *)realloc(rm_nf, size * sizeof(t_lpow));
        rm_nf_size = size;
    }
}

void resize_kf(uint size) {
    if (size > rm_kf_size) {
        if (size < rm_kf_size + 16)
            size = rm_kf_size + 16;
        rm_kf = (t_ppow *)realloc(rm_kf, size * sizeof(t_ppow));
        rm_kf_size = size;
    }
}

void _swap_r(e_results e) {
    t_results *rp = &ra[rm_base];
    t_results *rq = &ra[e];

    t_results temp = *rq;
    *rq = *rp;
    *rp = temp;
}

void _swapz_r(e_results e) {
    _swap_r(e);
    ra[rm_base].count = 0;
}

void init_rootmod(void) {
    memset(ra, 0, E_RESULTS_MAX * sizeof(t_results));
    resize_results(&ra[rm_base], 16);
    for (e_rmstash e = 0; e < E_RMSTASH_MAX; ++e)
        mpz_init(Z(e));
    resize_nf(16);
    resize_kf(16);
}

void done_rootmod(void) {
    free(rm_kf);
    free(rm_nf);
    for (e_rmstash e = 0; e < E_RMSTASH_MAX; ++e)
        mpz_clear(Z(e));
    for (uint i = 0; i < E_RESULTS_MAX; ++i) {
        t_results *rp = &ra[i];
        if (rp->size) {
            for (uint j = 0; j < rp->size; ++j)
                mpz_clear(rp->r[j]);
            free(rp->r);
            rp->size = 0;
        }
    }
}

int _mpz_comparator(const void *va, const void *vb) {
    return mpz_cmp(*(mpz_t *)va, *(mpz_t *)vb);
}

void save_result(t_results *rp, mpz_t r) {
    uint i = rp->count++;
    resize_results(rp, i + 1);
    mpz_set(rp->r[i], r);
}

void save_base(mpz_t r) {
    save_result(&ra[rm_base], r);
}

uint valuation(mpz_t result, mpz_t base, ulong p) {
    mpz_set(result, base);
    uint e = 0;
    while (mpz_divisible_ui_p(result, p)) {
        ++e;
        mpz_divexact_ui(result, result, p);
    }
    return e;
}

/* FIXME: write this without GMP */
ulong simple_invert(ulong d, ulong m) {
    mpz_set_ui(Z(si_d), d);
    mpz_set_ui(Z(si_m), m);
    if (mpz_invert(Z(si_d), Z(si_d), Z(si_m)))
        return mpz_get_ui(Z(si_d));
    return 0;
}

/* Given coprime n1, n2 and results arrays r1 (constant rm_base) and r2,
 * such that elements of r1 and r2 are kth roots of some a (mod n1) and
 * (mod n2) respectively, updates rm_base to be a new list of the kth
 * roots of a mod (n1 * n2).
 * CHECKME: is this implementing CRT? why is it not calling chinese()?
 */
void _allrootmod_cprod(e_results e2, mpz_t n1, mpz_t n2) {
    t_results *r = &ra[rm_base];
    t_results *r1 = &ra[armc_r1];
    t_results *r2 = &ra[e2];
    _swapz_r(armc_r1);
    resize_results(r, r1->count * r2->count);

    mpz_mul(Z(armc_n), n1, n2);
    if (!mpz_invert(Z(armc_inv), n1, n2)) {
        gmp_fprintf(stderr, "_allrootmod_cprod(%Zu, %Zu) has no inverse\n",
                n1, n2);
        exit(1);
    }
    for (uint i1 = 0; i1 < r1->count; ++i1) {
        mpz_t *z1 = &r1->r[i1];
        for (uint i2 = 0; i2 < r2->count; ++i2) {
            mpz_t *z2 = &r2->r[i2];
            /* save z1 + n1 * ((inv * (z2 - z1)) % n2)) % n */
            mpz_sub(Z(rm_r), *z2, *z1);
            mpz_mul(Z(rm_r), Z(rm_r), Z(armc_inv));
            mpz_mod(Z(rm_r), Z(rm_r), n2);
            mpz_mul(Z(rm_r), Z(rm_r), n1);
            mpz_mod(Z(rm_r), Z(rm_r), Z(armc_n));
            mpz_add(Z(rm_r), Z(rm_r), *z1);
            save_base(Z(rm_r));
        }
    }
    return;
}

/* "Tonelli-Shanks kth roots alternate version"
 * Given k prime, a coprime to p returns a root in rm_r, and zeta in tsp_z.
 */
void _ts_prime(mpz_t a, uint k, ulong p) {
    uint e = 0;
    ulong r = p - 1;
    while ((r % k) == 0) {
        ++e;
        r /= k;
    }

    ulong ke = (p - 1) / r;
    mpz_powm_ui(Z(rm_r), a, simple_invert(k % r, r), Z(rm_p));
    mpz_invert(Z(tsp_B), a, Z(rm_p));
    ulong ainv = mpz_get_ui(Z(tsp_B));
    mpz_powm_ui(Z(tsp_B), Z(rm_r), k, Z(rm_p));
    mpz_mul_ui(Z(tsp_B), Z(tsp_B), ainv);
    mpz_mod_ui(Z(tsp_B), Z(tsp_B), p);
    mpz_set_ui(Z(tsp_T), 2);
    mpz_set_ui(Z(tsp_y), 1);

    while (mpz_cmp_ui(Z(tsp_y), 1) == 0) {
        mpz_powm_ui(Z(tsp_z), Z(tsp_T), r, Z(rm_p));
        mpz_powm_ui(Z(tsp_y), Z(tsp_z), ke / k, Z(rm_p));
        mpz_add_ui(Z(tsp_T), Z(tsp_T), 1);
    }

    while (ke != k) {
        ke /= k;
        mpz_set(Z(tsp_T), Z(tsp_z));
        mpz_powm_ui(Z(tsp_z), Z(tsp_z), k, Z(rm_p));
        mpz_powm_ui(Z(tsp_A), Z(tsp_B), ke / k, Z(rm_p));
        while (mpz_cmp_ui(Z(tsp_A), 1) != 0) {
            mpz_mul(Z(rm_r), Z(rm_r), Z(tsp_T));
            mpz_mod_ui(Z(rm_r), Z(rm_r), p);
            mpz_mul(Z(tsp_B), Z(tsp_B), Z(tsp_z));
            mpz_mod_ui(Z(tsp_B), Z(tsp_B), p);
            mpz_mul(Z(tsp_A), Z(tsp_A), Z(tsp_y));
            mpz_mod_ui(Z(tsp_A), Z(tsp_A), p);
        }
    }
    return;
}

/* Helper function for _allrootmod_prime_power_r() to return
 * s + (a - s^k) / ks^{k - 1} (mod p^e).
 * Intermediate calculations may need to be done mod p^{e+1}, provided
 * by caller.
 */
void _eval_armpp(mpz_t s, mpz_t a, uint k, mpz_t px, mpz_t px2) {
    mpz_powm_ui(Z(earmpp_t), s, k - 1, px2);

    mpz_set(Z(earmpp_t1), a);
    mpz_submul(Z(earmpp_t1), Z(earmpp_t), s);
    mpz_mod(Z(earmpp_t1), Z(earmpp_t1), px2);

    mpz_mul_ui(Z(earmpp_t2), Z(earmpp_t), k);
    mpz_mod(Z(earmpp_t2), Z(earmpp_t2), px2);

    mpz_gcd(Z(earmpp_g), Z(earmpp_t1), Z(earmpp_t2));
    mpz_divexact(Z(earmpp_t1), Z(earmpp_t1), Z(earmpp_g));
    mpz_divexact(Z(earmpp_t2), Z(earmpp_t2), Z(earmpp_g));

    if (!mpz_invert(Z(earmpp_t2), Z(earmpp_t2), px)) {
        gmp_fprintf(stderr, "_eval_armpp() no inverse %Zu %% %Zu\n",
                Z(earmpp_t2), px);
        exit(1);
    }
    mpz_mul(Z(rm_r), Z(earmpp_t1), Z(earmpp_t2));
    /* mpz_mod(Z(rm_r), Z(rm_r), px); */
    mpz_add(Z(rm_r), Z(rm_r), s);
    mpz_mod(Z(rm_r), Z(rm_r), px);
    return;
}

/* Save the kth roots of a (mod p) given k prime.
 */
void _allrootmod_prime(mpz_t za, uint k, ulong p) {
    ulong a = mpz_fdiv_r_ui(Z(armp_a), za, p);
    if (p == 2 || a == 0) {
        mpz_set_ui(Z(rm_r), a);
        save_base(Z(rm_r));
        return;
    }

    ulong pm = p - 1;
    ulong g = simple_gcd(k, pm);

    /* If co-prime, there is exactly one root. */
    if (g == 1) {
        ulong inverse = simple_invert(k, pm);
        mpz_powm_ui(Z(rm_r), Z(armp_a), inverse, Z(rm_p));
        save_base(Z(rm_r));
        return;
    }

    /* Check generalized Euler's criterion. */
    mpz_powm_ui(Z(rm_r), Z(armp_a), pm / g, Z(rm_p));
    if (mpz_cmp_ui(Z(rm_r), 1) != 0)
        return;

    /* Special case p=3 for performance. */
    if (p == 3) {
        mpz_set_ui(Z(rm_r), 1);
        save_base(Z(rm_r));
        mpz_set_ui(Z(rm_r), 2);
        save_base(Z(rm_r));
        return;
    }

    /* Call a Tonelli-Shanks solver that also allows us to get all the roots.
     * Puts result in rm_r and zeta in tsp_z. */
    _ts_prime(Z(armp_a), k, p);
    if (mpz_sgn(Z(tsp_z)) == 0) {
        fprintf(stderr, "ts_prime(%lu, %u, %lu) gave zeta=0\n", a, k, p);
        exit(1);
    }
    mpz_powm_ui(Z(armp_t), Z(rm_r), k, Z(rm_p));
    if (mpz_cmp_ui(Z(armp_t), a) != 0) {
        fprintf(stderr, "ts_prime(%lu, %u, %lu) gave bad result\n", a, k, p);
        exit(1);
    }

    t_results *rp = &ra[rm_base];
    save_base(Z(rm_r));
    ulong r = mpz_get_ui(Z(rm_r));
    while (1) {
        mpz_mul(Z(rm_r), Z(rm_r), Z(tsp_z));
        mpz_mod_ui(Z(rm_r), Z(rm_r), p);
        if (mpz_cmp_ui(Z(rm_r), r) == 0)
            break;
        if (rp->count == k) {
            fprintf(stderr,
                "excess roots found for _allrootmod_prime(%lu, %u, %lu)",
                a, k, p
            );
            exit(1);
        }
        save_base(Z(rm_r));
    }
    return;
}

/* Save the kth roots of a (mod p^e) given k prime, a coprime to p.
 */
void _allrootmod_prime_power_r(mpz_t a, uint k, ulong p, uint e) {
    if (e == 1) {
        _allrootmod_prime(a, k, p);
        return;
    }
    uint e2 = (p > 2 || e < 5) ? ((e + 1) >> 1) : ((e + 3) >> 1);
    /* Recurse until e = 1, then walk back up the stack */
    _allrootmod_prime_power_r(a, k, p, e2);

    t_results *rp = &ra[rm_base];
    t_results *r2 = &ra[armppr_copy];
    if (rp->count == 0)
        return;
    _swapz_r(armppr_copy);

    mpz_pow_ui(Z(armppr_px), Z(rm_p), e);
    if (k != p) {
        for (uint ri = 0; ri < r2->count; ++ri) {
            _eval_armpp(r2->r[ri], a, k, Z(armppr_px), Z(armppr_px));
            save_base(Z(rm_r));
        }
        return;
    }

    mpz_mod(Z(armppr_a), a, Z(armppr_px));
    mpz_pow_ui(Z(armppr_px2), Z(rm_p), e + 1);
    for (uint ri = 0; ri < r2->count; ++ri) {
        _eval_armpp(r2->r[ri], a, k, Z(armppr_px), Z(armppr_px2));
        /* check if it's a solution */
        mpz_powm_ui(Z(armppr_z), Z(rm_r), k, Z(armppr_px));
        if (mpz_cmp(Z(armppr_z), Z(armppr_a)) == 0)
            save_base(Z(rm_r));
    }
    /* now deduplicate */
    if (rp->count == 0)
        return;
    _swapz_r(armppr_copy);
    qsort(r2->r, r2->count, sizeof(mpz_t), &_mpz_comparator);
    mpz_t *prev = NULL;
    for (uint ri = 0; ri < r2->count; ++ri)
        if (prev == NULL || mpz_cmp(*prev, r2->r[ri]) != 0) {
            prev = &r2->r[ri];
            save_base(*prev);
        }
    return;
}

/* Save the kth roots of 0 (mod p^e = px) given k prime.
 */
void _allrootzero_prime_power(uint k, ulong p, uint e, mpz_t px) {
    uint t = (e - 1) / k + 1;
    mpz_pow_ui(Z(arzpp_px), Z(rm_p), e - t);
    if (!mpz_fits_uint_p(Z(arzpp_px))) {
        fprintf(stderr, "_allrootzero_prime_power() overflow %lu^%u\n",
                p, e - t);
        exit(1);
    }
    uint r = mpz_get_ui(Z(arzpp_px));
    mpz_pow_ui(Z(arzpp_px), Z(rm_p), t);
    for (uint i = 0; i < r; ++i) {
        mpz_mul_ui(Z(rm_r), Z(arzpp_px), i);
        mpz_mod(Z(rm_r), Z(rm_r), px);
        save_base(Z(rm_r));
    }
    return;
}

/* Save the kth roots of a (mod p^e = px) given k prime.
 */
void _allrootmod_prime_power(mpz_t a, uint k, ulong p, uint e, mpz_t px) {
    if (e == 1) {
        _allrootmod_prime(a, k, p);
        return;
    }

    uint v = valuation(Z(rm_r), a, p);
    if (v == 0) {
        _allrootmod_prime_power_r(a, k, p, e);
        return;
    }
    if (v >= e) {
        _allrootzero_prime_power(k, p, e, px);
        return;
    }
    if (v % k)
        return;
    uint m = v / k;

    /* We now know that p^{mk} divides a leaving a' coprime to p.
     * We now want roots R = { r_i } of a' mod p^{e - mk}.
     * Solutions will then be of the form r_i p^m + jp^{e - m(k - 1)} (mod p^e)
     * for all r_i in R and all j: 0 <= j < p^{m(k - 1)}.
     */

    /* Note: we rely on armpp_r copying rm_r before overwriting it */
    _allrootmod_prime_power_r(Z(rm_r), k, p, e - m * k);

    t_results *rp = &ra[rm_base];
    t_results *r2 = &ra[armpp_copy];
    if (rp->count == 0)
        return;
    _swapz_r(armpp_copy);

    mpz_pow_ui(Z(armpp_px), Z(rm_p), m * (k - 1));
    if (!mpz_fits_uint_p(Z(armpp_px))) {
        fprintf(stderr, "_allrootmod_prime_power() overflow %lu^%u\n",
                p, m * (k - 1));
        exit(1);
    }
    uint range = mpz_get_ui(Z(armpp_px));
    resize_results(rp, range * r2->count);

    mpz_pow_ui(Z(armpp_px), Z(rm_p), m);
    for (uint ri = 0; ri < r2->count; ++ri) {
        mpz_t *r = &r2->r[ri];
        mpz_mul(*r, *r, Z(armpp_px));
        mpz_mod(*r, *r, px);
    }

    mpz_pow_ui(Z(armpp_px), Z(rm_p), e - m * (k - 1));
    for (uint ri = 0; ri < r2->count; ++ri) {
        mpz_t *r = &r2->r[ri];
        for (uint i = 0; i < range; ++i) {
            mpz_mul_ui(Z(rm_r), Z(armpp_px), i);
            mpz_add(Z(rm_r), Z(rm_r), *r);
            mpz_mod(Z(rm_r), Z(rm_r), px);
            save_base(Z(rm_r));
        }
    }
    return;
}

/* Save (append) kth roots of a (mod n) given k prime, and array nf with
 * size nfc of factors of n.
 */
void _allrootmod_kprime(mpz_t a, uint k, mpz_t n, t_lpow *nf, uint nfc) {
#if 0
    /* I think the original perl source added this because it had
     * allsqrtmodfact() already written, but I don't think it gains enough
     * to justify translating. */
    if (k == 2) {
        _allsqrtmodfact(a, n, nf, nfc);
        return;
    }
#endif
    /* save the results we should append to */
    t_results *r = &ra[rm_base];
    t_results *stash = &ra[armkp_stash];
    t_results *base = &ra[armkp_base];
    _swapz_r(armkp_stash);

    mpz_set_ui(Z(armkp_m), 1);
    for (uint nfi = 0; nfi < nfc; ++nfi) {
        /* loop: given previous results at r, find new results and
         * given them, leaving the augmented list at r; when nfi == 0
         * there are no previous results */
        if (nfi > 0)
            _swapz_r(armkp_base);

        ulong p = nf[nfi].p;
        uint e = nf[nfi].e;
        mpz_set_ui(Z(rm_p), p);
        if (e == 1) {
            mpz_set_ui(Z(armkp_px), p);
            _allrootmod_prime(a, k, p);
        } else {
            mpz_pow_ui(Z(armkp_px), Z(rm_p), e);
            _allrootmod_prime_power(a, k, p, e, Z(armkp_px));
        }
        if (r->count == 0)
            goto armkp_abort;
        if (nfi > 0)
            _allrootmod_cprod(armkp_base, Z(armkp_px), Z(armkp_m));
        mpz_mul(Z(armkp_m), Z(armkp_m), Z(armkp_px));
    }
    resize_results(stash, stash->count + r->count);
    for (uint i = 0; i < r->count; ++i)
        save_result(stash, r->r[i]);
  armkp_abort:
    _swap_r(armkp_stash);
    return;
}

/* Given positive even integer k and positive mpz_t a, n, where n is
 * a product of powers of primes fitting in ulong, constructs a list
 * of kth roots x of a (mod n) having 0 <= x < n.
 *
 * Returns the number of roots found, and writes the location of the
 * (singleton) array of roots (sorted ascending) into *result.
 *
 * It is the caller's responsibility to avoid calling allrootmod() again
 * before they have finished looking at the results.
 *
 */
uint allrootmod(mpz_t **result, mpz_t a, uint k, mpz_t n) {
    t_results *rp = &ra[rm_base];
    rp->count = 0;

    mpz_mod(Z(arm_a), a, n);
    if (mpz_cmp_ui(n, 2) < 0) {
        save_base(Z(arm_a));
        goto arm_done;
    }

    uint nfc = 0;
    factor_state fs;
    fs_init(&fs);
    mpz_set(fs.n, n);
    while (factor_one(&fs)) {
        resize_nf(nfc + 1);
        rm_nf[nfc].p = mpz_get_ui(fs.f);
        rm_nf[nfc].e = fs.e;
        ++nfc;
    }
    if (fs.state != FS_TERM) {
        gmp_fprintf(stderr, "In allrootmod failed to factorize n=%Zu\n", fs.n);
        exit(1);
    }
    fs_clear(&fs);

    /* now similarly factorize k */
    uint kfc = 0;
    fs_init(&fs); 
    mpz_set_ui(fs.n, k);
    while (factor_one(&fs)) {
        resize_kf(kfc + 1);
        rm_kf[kfc].p = mpz_get_ui(fs.f);
        rm_kf[kfc].e = fs.e;
        ++kfc;
    }
    fs_clear(&fs);

    if (kfc == 1 && rm_kf[0].e == 1) {
        _allrootmod_kprime(Z(arm_a), k, n, rm_nf, nfc);
        goto arm_done;
    }

    uint ki = 0;
    uint ke = 0;
    t_results *rq = &ra[arm_scratch];
    while (1) {
        _swapz_r(arm_scratch);
        for (uint ri = 0; ri < rq->count; ++ri)
            _allrootmod_kprime(rq->r[ri], rm_kf[ki].p, n, rm_nf, nfc);
        ++ke;
        if (ke >= rm_kf[ki].e) {
            ++ki;
            ke = 0;
            if (ki >= kfc)
                break;
        }
    }

  arm_done:
    if (rp->count > 1)
        qsort(rp->r, rp->count, sizeof(mpz_t), &_mpz_comparator);
    *result = rp->r;
    return rp->count;
}

#undef Z
