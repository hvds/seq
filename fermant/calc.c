#include <gmp.h>
#include <stdio.h>

#include "calc.h"
#include "int.h"
#include "frag.h"
#include "path.h"
#include "source.h"
#include "diag.h"

typedef unsigned char pow_t;
typedef unsigned char lcsize_t;

typedef struct {
    pow_t pow;
    lincomz_t lcp[0];
} term_t;

typedef struct {
    uint count;
    uint alloc;
    term_t term[0];
} mul_t;
typedef uint constid_t;
typedef uint mulid_t;

typedef struct {
    constid_t c;
    mulid_t m;
} cmul_t;

typedef struct {
    uint count;
    uint alloc;
    cmul_t cmul[0];
} expr_t;
typedef uint exprid_t;
/* we actually only need 2, but debugging with 3 gives some provenance */
#define NUM_EXPRS 3

mpq_t *path_total = NULL;
lincom_t *path_lc = NULL;
mpq_t constq;       /* expr_const */
mpq_t cdsiq;        /* const_div_si */
mpq_t ddq, dd2q;    /* do_distrib */
mpq_t ieq;          /* inteval */
mpq_t totalq;       /* exported, set up by report_total */

term_t *terms = NULL;
uint nterms = 0;
uint sizeterms = 0;
mul_t **muls = NULL;
uint nmuls = 0;
uint sizemuls = 0;
mpq_t *consts = NULL;
uint nconsts = 0;
uint sizeconsts = 0;
expr_t *exprs[NUM_EXPRS];

static inline void mpq_set_sisi(mpq_t q, int num, int den) {
    /* mpq_set_si accepts only _unsigned_ den */
    if (den < 0)
        mpq_set_si(q, -num, -den);
    else
        mpq_set_si(q, num, den);
}

static inline lincom_t *calc_path_lc(uint pi) {
    return (lincom_t *)(add_p(path_lc, (size_t)pi * (size_t)lc_size()));
}

static inline uint term_size(void) {
    return sizeof(pow_t) + sizeof(lcsize_t) + lc_size();
}

static inline pow_t term_pow(term_t *tp) {
    return tp->pow;
}

static inline void term_pow_set(term_t *tp, pow_t pow) {
    tp->pow = pow;
}

static inline lincom_t *term_lc(term_t *tp) {
    return (lincom_t *)&tp->lcp[0];
}

static inline void term_copy(term_t *tdst, term_t *tsrc) {
    term_pow_set(tdst, term_pow(tsrc));
    lc_copy(term_lc(tdst), term_lc(tsrc));
}

static inline int term_cmp(term_t *ta, term_t *tb, uint vmax) {
    if (ta->pow != tb->pow)
        return (int)ta->pow - (int)tb->pow;
    return lc_cmp(term_lc(ta), term_lc(tb), vmax);
}

static inline uint mul_size(uint terms) {
    return sizeof(mul_t) + terms * term_size();
}

static inline void reset_muls(void) {
    nmuls = 0;
}

static inline void resize_muls(uint extra) {
    if (nmuls + extra <= sizemuls)
        return;
    uint newsize = 3 * sizemuls / 2;
    while (nmuls + extra > newsize)
        newsize += 100;
    muls = realloc(muls, newsize * sizeof(mul_t *));
    for (mulid_t mi = sizemuls; mi < newsize; ++mi) {
        muls[mi] = calloc(1, mul_size(10));
        muls[mi]->alloc = 10;
    }
    sizemuls = newsize;
}

static inline mulid_t new_mul(void) {
    resize_muls(1);
    return nmuls++;
}

static inline void free_mul(mul_t *mp) {
    free(mp);
}

static inline mul_t *mul_p(mulid_t mi) {
    return muls[mi];
}

static inline void mul_resize(mulid_t mi, uint count) {
    if (mul_p(mi)->alloc >= count)
        return;
    uint newsize = mul_p(mi)->alloc * 2;
    if (newsize < count)
        newsize = count;
    mul_t *mp = realloc(mul_p(mi), mul_size(newsize));
    mp->alloc = newsize;
    muls[mi] = mp;
}

static inline uint mul_count(mulid_t mi) {
    return mul_p(mi)->count;
}

static inline term_t *mul_term(mulid_t mi, uint off) {
    return (term_t *)add_p(&mul_p(mi)->term[0],
            (size_t)off * (size_t)term_size());
}

static inline lincom_t *mul_lc(mulid_t mi, uint off) {
    return term_lc(mul_term(mi, off));
}

static inline void mul_count_set(mulid_t mi, uint count) {
    mul_resize(mi, count);
    mul_p(mi)->count = count;
}

static inline void mul_copy(mulid_t mdst, mulid_t msrc) {
    uint count = mul_count(msrc);
    mul_count_set(mdst, count);
    memcpy(mul_term(mdst, 0), mul_term(msrc, 0), count * term_size());
}

static inline void mul_remove(mulid_t mi, uint off) {
    uint last = --mul_p(mi)->count;
    if (off < last)
        term_copy(mul_term(mi, off), mul_term(mi, last));
}

/* must be called only with non-const lc */
static inline uint mul_mul_lc(mulid_t mi, lincom_t *lc, uint vmax, pow_t p) {
    for (uint vi = 1; vi <= vmax; ++vi)
        if (lc_get(lc, vi) != 0)
            goto non_const;
    assert(0);
    
  non_const:
    for (uint ti = 0; ti < mul_count(mi); ++ti) {
        term_t *tp = mul_term(mi, ti);
        if (lc_cmp(term_lc(tp), lc, vmax) != 0)
            continue;
        term_pow_set(tp, term_pow(tp) + p);
        return ti;
    }
    uint tn = mul_count(mi);
    mul_count_set(mi, tn + 1);
    term_t *tp = mul_term(mi, tn);
    term_pow_set(tp, p);
    lc_copy(term_lc(tp), lc);
    return tn;
}

/* this could be easier if mul_t terms were sorted, but we expect the
 * average number of terms to be between 1 and 2, so it seems unlikely
 * to be worth the effort.
 */
static inline int mul_match(mulid_t mi, mulid_t mj, uint vmax) {
    uint count = mul_count(mi);
    if (mul_count(mj) != count)
        return 0;
    for (uint i = 0; i < count; ++i) {
        term_t *ti = mul_term(mi, i);
        for (uint j = 0; j < count; ++j) {
            term_t *tj = mul_term(mj, j);
            if (term_cmp(ti, tj, vmax) == 0)
                goto matched_term;
        }
        return 0;
      matched_term:
        ;
    }
    return 1;
}

static inline mpq_t *const_mpq(constid_t cmi) {
    return &consts[cmi];
}

static inline void reset_consts(void) {
    nconsts = 0;
}

static inline void resize_consts(uint extra) {
    if (nconsts + extra <= sizeconsts)
        return;
    uint newsize = 3 * sizeconsts / 2;
    while (nconsts + extra > newsize)
        newsize += 100;
    consts = realloc(consts, newsize * sizeof(mpq_t));
    for (constid_t ci = sizeconsts; ci < newsize; ++ci)
        mpq_init(*const_mpq(ci));
    sizeconsts = newsize;
}

static inline constid_t new_const(void) {
    resize_consts(1);
    constid_t ci = nconsts++;
    mpq_set_ui(*const_mpq(ci), 1, 1);
    return ci;
}

static inline void const_div_si(constid_t cmi, int div) {
    mpq_set_sisi(cdsiq, 1, div);
    mpq_mul(*const_mpq(cmi), *const_mpq(cmi), cdsiq);
}

static inline constid_t cmul_const(cmul_t *cmp) {
    return cmp->c;
}

static inline void cmul_const_set_ui(cmul_t *cmp, uint p, uint q) {
    mpq_set_ui(*const_mpq(cmp->c), p, q);
}

static inline mulid_t cmul_mul(cmul_t *cmp) {
    return cmp->m;
}

int cmul_comparator(const void *va, const void *vb) {
    const cmul_t *cma = (const cmul_t *)va;
    const cmul_t *cmb = (const cmul_t *)vb;
    int count = mul_count(cma->m);
    if (count != mul_count(cmb->m))
        return count - mul_count(cmb->m);
    return memcmp(mul_term(cma->m, 0), mul_term(cmb->m, 0),
            count * term_size());
}

static inline uint expr_size(uint count) {
    return sizeof(expr_t) + count * sizeof(cmul_t);
}

static inline expr_t *expr_p(exprid_t ei) {
    return exprs[ei];
}

expr_t *new_expr(uint alloc) {
    expr_t *e = malloc(sizeof(expr_t) + alloc * sizeof(cmul_t));
    e->count = 0;
    e->alloc = alloc;
    return e;
}

void free_expr(exprid_t ei) {
    free(expr_p(ei));
}

static inline void expr_resize(exprid_t ei, uint count) {
    if (expr_p(ei)->alloc >= count)
        return;
    uint newsize = expr_p(ei)->alloc * 2;
    if (newsize < count)
        newsize = count;
    expr_t *ep = realloc(expr_p(ei), expr_size(newsize));
    ep->alloc = newsize;
    exprs[ei] = ep;
}

static inline uint expr_count(exprid_t ei) {
    return expr_p(ei)->count;
}

static inline cmul_t *expr_cmul(exprid_t ei, uint off) {
    return &(expr_p(ei)->cmul[off]);
}

static inline void expr_count_set(exprid_t ei, uint count) {
    uint oc = expr_count(ei);
    expr_resize(ei, count);
    expr_p(ei)->count = count;
    while (oc < count) {
        expr_cmul(ei, oc)->m = new_mul();
        expr_cmul(ei, oc)->c = new_const();
        ++oc;
    }
}

static inline mulid_t expr_mul(exprid_t ei, uint off) {
    return expr_p(ei)->cmul[off].m;
}

static inline mpq_t *expr_const(exprid_t ei) {
    if (expr_count(ei)) {
        assert(expr_count(ei) == 1);
        mpq_set(constq, *const_mpq(cmul_const(expr_cmul(ei, 0))));
        assert(mul_count(expr_mul(ei, 0)) == 0);
    } else {
        mpq_set_ui(constq, 0, 1);
    }
    return &constq;
}

static inline void expr_remove(exprid_t ei, uint off) {
    uint last = --expr_p(ei)->count;
    if (off < last)
        *expr_cmul(ei, off) = *expr_cmul(ei, last);
}

static inline uint expr_match(exprid_t ei, mulid_t mi, uint vmax) {
    uint count = expr_count(ei);
    for (uint i = 0; i < count; ++i)
        if (mul_match(cmul_mul(expr_cmul(ei, i)), mi, vmax))
            return i;
    return count;
}

static inline void expr_add(exprid_t ei, mpq_t c, mulid_t mi, uint vmax) {
    if (mpq_sgn(c) == 0)
        return;

    uint off = expr_match(ei, mi, vmax);
    if (off < expr_count(ei)) {
        mpq_t *q = const_mpq(cmul_const(expr_cmul(ei, off)));
        mpq_add(*q, *q, c);
    } else {
        expr_count_set(ei, off + 1);
        cmul_t *cmp = expr_cmul(ei, off);
        cmp->c = new_const();
        mpq_set(*const_mpq(cmp->c), c);
        cmp->m = new_mul();
        mul_copy(cmp->m, mi);
    }
}

uint mul_disp(char *buf, uint bufsize, mulid_t mi, uint vi) {
    uint pos = 0;
    uint count = mul_count(mi);
    char vbuf[lc_dumpsize()];
    for (uint i = 0; i < count; ++i) {
        term_t *tp = mul_term(mi, i);
        pow_t p = term_pow(tp);
        uint vpos = lc_disp(vbuf, sizeof(vbuf), term_lc(tp), vi);
        if (vpos > 1)
            pos += snprintf(buf + pos, bufsize - pos, "(%s)", vbuf);
        else
            pos += snprintf(buf + pos, bufsize - pos, "%s", vbuf);
        if (p != 1)
            pos += snprintf(buf + pos, bufsize - pos, "^%u", p);
    }
    return pos;
}

uint cmul_disp(char *buf, uint bufsize, cmul_t *cmp, uint vi) {
    mpq_t *c = const_mpq(cmul_const(cmp));
    if (mpz_cmp_ui(mpq_numref(*c), 0) == 0)
        return snprintf(buf, bufsize, "0");
    uint pos = gmp_snprintf(buf, bufsize, "%Qd", *c);
    if (pos < bufsize)
        pos += snprintf(buf + pos, bufsize - pos, " ");
    return pos + mul_disp(buf + pos, bufsize - pos, cmul_mul(cmp), vi);
}

#define EXPR_DUMPSIZE 4096
void expr_dump(exprid_t ei, uint vi) {
    char buf[EXPR_DUMPSIZE];
    for (uint i = 0; i < expr_count(ei); ++i) {
        cmul_t *cmp = expr_cmul(ei, i);
        cmul_disp(buf, sizeof(buf), cmp, vi);
        if (i && buf[0] != '-')
            fprintf(stderr, "+%s", buf);
        else
            fprintf(stderr, "%s", buf);
    }
    fprintf(stderr, "\n");
}

/* take num_paths to ensure paths have been initialized */
void init_calc(uint num_paths) {
    path_total = malloc(num_paths * sizeof(mpq_t));
    path_lc = calloc(num_paths, lc_size());
    for (uint i = 0; i < num_paths; ++i) {
        mpq_init(path_total[i]);
        path_t p = path_p(i);
        lincom_t *lc = calc_path_lc(i);
        while (p) {
            uint vi = path_first(p) + 1;
            lc_set(lc, vi, 1);
            /* FIXME: path.[ch] should own this logic */
            p &= ~(1 << (vi - 1));
        }
    }
    for (uint i = 0; i < NUM_EXPRS; ++i)
        exprs[i] = new_expr(20);
    mpq_init(constq);
    mpq_init(cdsiq);
    mpq_init(ddq);
    mpq_init(dd2q);
    mpq_init(ieq);
    mpq_init(totalq);
}

void done_calc(void) {
    for (uint i = 0; i < npaths; ++i)
        mpq_clear(path_total[i]);
    free(path_total);
    free(path_lc);
    for (exprid_t i = 0; i < NUM_EXPRS; ++i)
        free_expr(i);
    for (constid_t ci = 0; ci < sizeconsts; ++ci)
        mpq_clear(*const_mpq(ci));
    for (mulid_t mi = 0; mi < sizemuls; ++mi)
        free(mul_p(mi));
    free(muls);
    free(consts);
    mpq_clear(constq);
    mpq_clear(cdsiq);
    mpq_clear(ddq);
    mpq_clear(dd2q);
    mpq_clear(ieq);
    mpq_clear(totalq);
}

exprid_t init_expr(path_t pi) {
    reset_muls();
    reset_consts();
    exprid_t ei = 0;
    expr_count_set(ei, 1);
    cmul_t *cmp = expr_cmul(ei, 0);
    cmp->c = new_const();
    cmul_const_set_ui(cmp, 1, 1);
    cmp->m = new_mul();
    mulid_t mi = cmul_mul(cmp);
    mul_count_set(mi, 1);
    lc_copy(mul_lc(mi, 0), calc_path_lc(pi));
    term_pow_set(mul_term(mi, 0), 1);
    return ei;
}

void do_distrib(exprid_t e0, uint off, uint ti, uint tj, uint vi) {
    LINCOM_ALLOC(lck);
    constid_t ci = cmul_const(expr_cmul(e0, off));
    mulid_t mi = cmul_mul(expr_cmul(e0, off));
    expr_remove(e0, off);
    pow_t pj = term_pow(mul_term(mi, tj));
    lincom_t *lci = term_lc(mul_term(mi, ti));
    lincom_t *lcj = term_lc(mul_term(mi, tj));
    /* given l_i, l_j have terms a v_i, b v_i with g = gcd(a, b), we first
     * denormalise l_j to a/g l_j, then split it as
     * b/g l_i + [c l_k = (a/g l_j - b/g l_i)]; distributing that gives terms
     * (b/g)^p c^{pj-p} l_i^p l_k^{p_j-p} C(p_j, p) with an additional constant
     * factor (a/g)^{-p_j}.
     */
#if 0

we may end up with the 'l_k' co-term coinciding with an existing term, so
we should a) allocate room for l_k on the stack, and generate it into there;
b) _always_ mul_remove tj; c) multiply l_k^pj back in, recording tk as the
result; then proceed as before moving powers from tk to ti.

#endif
    int a = lc_get(lci, vi);
    int b = lc_get(lcj, vi);
    uint g = ugcd(abs(a), abs(b));
    int aog = a / g;
    int bog = b / g;
    if (aog == 1)
        lc_copy(lck, lcj);
    else {
        for (uint vj = 0; vj <= vi; ++vj)
            /* FIXME: bounds check */
            lc_set(lck, vj, aog * lc_get(lcj, vj));
        mpq_set_sisi(ddq, 1, aog);
        mpz_pow_ui(mpq_denref(ddq), mpq_denref(ddq), pj);
        mpq_mul(*const_mpq(ci), *const_mpq(ci), ddq);
    }
    for (uint vj = 0; vj <= vi; ++vj)
        lc_set(lck, vj, lc_get(lck, vj) - bog * lc_get(lci, vj));
    int c = lc_norm(lck, vi);
    int is_const = lc_is_const(lck, vi);

    mul_remove(mi, tj);
    /* if ti was last term, it will have moved to fill the gap for tj */
    if (ti == mul_count(mi))
        ti = tj;
    uint tk;
    if (!is_const)
        tk = mul_mul_lc(mi, lck, vi, pj);

    mpq_set_si(ddq, c, 1);
    mpz_pow_ui(mpq_numref(ddq), mpq_numref(ddq), pj);
    mpq_mul(*const_mpq(ci), *const_mpq(ci), ddq);

    /* expr_add() within the loop may trigger reallocations, so must not
     * rely on pointers here, only id types */
    mpq_set_sisi(ddq, bog, c);
    for (uint p = 0; p <= pj; ++p) {
        mpq_set_si(dd2q, icomb(pj, p), 1);
        mpq_mul(dd2q, dd2q, *const_mpq(ci));
        expr_add(e0, dd2q, mi, vi);
        if (p < pj) {
            term_pow_set(mul_term(mi, ti), term_pow(mul_term(mi, ti)) + 1);
            if (!is_const) {
                term_pow_set(mul_term(mi, tk), term_pow(mul_term(mi, tk)) - 1);
                /* can only happen on last loop iteration */
                if (term_pow(mul_term(mi, tk)) == 0)
                    mul_remove(mi, tk);
            }
            mpq_mul(*const_mpq(ci), *const_mpq(ci), ddq);
        }
    }
}

/* If there is any mul_t in the expr with more than one term dependent
 * on vi, distribute out one such term iteratively until every mul_t has
 * at most one dependent term (and is therefore integrable).
 * Beware that expr_count() can be increased or decreased by do_distrib().
 */
void find_distrib(exprid_t e0, uint vi) {
    for (uint off = 0; off < expr_count(e0); ++off) {
        mulid_t mi = expr_mul(e0, off);
        uint count = 0;
        uint ta, tb;
        for (uint ti = 0; ti < mul_count(mi); ++ti) {
            if (!lc_get(mul_lc(mi, ti), vi))
                continue;
            if (count == 0)
                ta = ti;
            if (++count > 1) {
                tb = ti;
                break;
            }
        }
        if (count > 1) {
            if (term_cmp(mul_term(mi, ta), mul_term(mi, tb), vi) >= 0)
                do_distrib(e0, off, ta, tb, vi);
            else
                do_distrib(e0, off, tb, ta, vi);
            --off;  /* redo */
            continue;
        }
    }
}

/* each mul_t is directly integrable, we can even do it in place */
void safe_integrate(exprid_t e0, uint vi) {
    /* make a lincom consisting of just vi */
    /* FIXME: have a cache of these */
    char vlc[lc_size()];
    lincom_t *lc = (lincom_t *)&vlc[0];
    memset(lc, 0, lc_size());
    lc_set(lc, vi, 1);

    uint ecount = expr_count(e0);
    for (uint off = 0; off < ecount; ++off) {
        mulid_t mi = expr_mul(e0, off);
        for (uint ti = 0; ti < mul_count(mi); ++ti) {
            int coeff = (int)lc_get(mul_lc(mi, ti), vi);
            if (!coeff)
                continue;
            /* this is the dependent term */
            pow_t pow = term_pow(mul_term(mi, ti)) + 1;
            term_pow_set(mul_term(mi, ti), pow);
            const_div_si(cmul_const(expr_cmul(e0, off)), coeff * (int)pow);
            goto integrate_done;
        }
        /* independent of vi, so just multiply by vi */
        mul_mul_lc(mi, lc, vi, 1);
      integrate_done:
        ;
    }
}

static inline void eval_lim(
    limit_t *lim, int dir,
    mulid_t msrc, mulid_t mdst, uint targ,
    mpq_t *cq,
    lincom_t *lcsrc, lincom_t *lcdst,
    int coeff, pow_t p,
    uint vi, cmul_t *cmp,
    exprid_t e
) {
    mul_copy(mdst, msrc);
    mul_remove(mdst, targ);
    int clim = coeff * limitp_num(lim);
    int csrc = limitp_den(lim);
    mpq_set_si(*cq, 1, csrc);
    lincom_t *lclim = limitp_lc(lim);
    /* given (x+yv)^a, we substitute p/q z for v by calculating
     * (qx+pyz) and adjusting external const by 1/q^a */
    for (uint j = 0; j < vi; ++j)
        lc_set(lcdst, j, csrc * lc_get(lcsrc, j) + clim * lc_get(lclim, j));
    int g = lc_norm(lcdst, vi - 1);
    if (!lc_is_const(lcdst, vi - 1))
        mul_mul_lc(mdst, lcdst, vi - 1, p);
    if (g != 1) {
        mpz_mul_si(mpq_numref(*cq), mpq_numref(*cq), g);
        mpq_canonicalize(*cq);
    }
    if (p != 1) {
        mpz_pow_ui(mpq_numref(*cq), mpq_numref(*cq), p);
        mpz_pow_ui(mpq_denref(*cq), mpq_denref(*cq), p);
    }
    mpq_mul(*cq, *cq, *const_mpq(cmul_const(cmp)));
    if (dir < 0)
        mpq_neg(*cq, *cq);
    expr_add(e, *cq, mdst, vi - 1);
}

exprid_t inteval(exprid_t e0, range_t *rp, uint vi) {
    uint ecount = expr_count(e0);
    exprid_t e1 = (e0 + 1) % NUM_EXPRS;
    expr_count_set(e1, 0);
    mulid_t m1 = new_mul();
    char vlc[lc_size()];
    lincom_t *lc1 = (lincom_t *)&vlc[0];
    for (uint i = 0; i < ecount; ++i) {
        cmul_t *cmp = expr_cmul(e0, i);
        mulid_t m0 = cmul_mul(cmp);
        uint mcount = mul_count(m0);
        uint seen = mcount;
        for (uint j = 0; j < mcount; ++j) {
            if (lc_get(term_lc(mul_term(m0, j)), vi) == 0)
                continue;
            seen = j;
            break;
        }
        /* there must be exactly one term dependent on vi */
        /* FIXME: we could cache this when we integrate (or combine
         * inteval with integrate) */
        assert(seen < mcount);

        term_t *t0 = mul_term(m0, seen);
        lincom_t *lc0 = term_lc(t0);
        pow_t p = term_pow(t0);
        int coeff = lc_get(lc0, vi);

        eval_lim(range_high(rp), 1,
                m0, m1, seen, &ieq, lc0, lc1, coeff, p, vi, cmp, e1);
        eval_lim(range_low(rp), -1,
                m0, m1, seen, &ieq, lc0, lc1, coeff, p, vi, cmp, e1);
    }
    return e1;
}

exprid_t do_integrate(exprid_t e0, range_t *rp, uint vi) {
    find_distrib(e0, vi);
    if (debug_integrate) {
        fprintf(stderr, "  distrib %c to ", 'a' + vi - 1);
        expr_dump(e0, vi);
    }
    safe_integrate(e0, vi);
    if (debug_integrate) {
        fprintf(stderr, "  integrate %c to ", 'a' + vi - 1);
        expr_dump(e0, vi);
    }
    exprid_t e1 = inteval(e0, rp, vi);
    if (debug_integrate) {
        fprintf(stderr, "  inteval %c to ", 'a' + vi - 1);
        expr_dump(e1, vi - 1);
    }
    return e1;
}

static inline void alloc_result(exprid_t e, path_t pi) {
    mpq_t *c = expr_const(e);
    mpq_add(path_total[pi], path_total[pi], *c);
}

void report_total(void) {
    for (uint i = 0; i < npaths; ++i) {
        report("000 path[%d] = %Qd\n", i, path_total[i]);
        mpq_add(totalq, totalq, path_total[i]);
    }
}

void integrate(fid_t fi) {
    pathset_t ps = frag_ps(fi);
    if (pathset_count(ps) != 1)
        fail("panic: ps=%#x for fi %u\n", ps, fi);
    path_t pi = pathset_first(ps);

    char buf[frag_dumpsize()];
    if (debug_integrate) {
        frag_disp(buf, sizeof(buf), fi);
        fprintf(stderr, "integrate %s\n", buf);
    }

    exprid_t e = init_expr(pi);
    for (uint vi = nv; vi > 0; --vi) {
        if (debug_integrate)
            expr_dump(e, vi);
        e = do_integrate(e, frag_range(fi, vi), vi);
    }
    if (debug_integrate) {
        fprintf(stderr, "result ");
        expr_dump(e, 0);
    }
    alloc_result(e, pi);
}

void integrate_path(uint pi) {
    for (uint ri = 0; ri < nresolve; ++ri) {
        if (resolves[ri].pi != pi && resolves[ri].pj != pi)
            continue;
        if (!integrate_open(ri, pi))
            continue;
        uint count = 0;
        while (1) {
            reset_frags();
            ++count;
            fid_t fi = new_frag();
            if (!read_frag(fi))
                break;
            if (need_diag) {
                need_diag = 0;
                diag("int(%i) %u: %u", pi, ri, count);
            }
            integrate(fi);
            if (need_log) {
                integrate_checkpoint();
                need_log = 0;
            }
        }
        integrate_close(ri, pi);
    }
}
