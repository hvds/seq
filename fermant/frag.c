#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "frag.h"
#include "int.h"
#include "path.h"
#include "num.h"
#include "lincom.h"
#include "limit.h"
#include "diag.h"

extern inline uint range_size(void);
extern inline uint frag_size(void);
extern inline frag_t *frag_p(fid_t f);
extern inline void resize_frags(uint extra);
extern inline void free_frag(frag_t *fp);
extern inline fid_t new_frag(void);
extern inline fid_t frag_dup(fid_t fi);
extern inline range_t *frag_range(fid_t f, uint vi);

crange_t FRC[2];
frag_t **frags = NULL;
uint nfrags = 0;
uint sizefrags = 0;
uint next_fragid = 0;

uint frag_dumpsize(void) {
    /* "frag 999(998) 0x18: [0, a+b]; [...]...\n" */
    assert(sizeof(fid_t) <= 4);
    return 30 + (npaths + 3)/4 + (limit_dumpsize() * 2 + 6) * nv + 2;
}

uint fragp_disp(char *buf, uint bufsize, frag_t *fp) {
    uint pos = 0;
    pos += snprintf(&buf[pos], bufsize - pos, "frag %u(%u) 0x%x: ",
            fp->fid, fp->parent, fp->ps);
    for (uint vi = 1; vi <= nv; ++vi) {
        pos += snprintf(&buf[pos], bufsize - pos, "[");
        pos += limit_disp(&buf[pos], bufsize - pos, fp->r[vi-1].low);
        pos += snprintf(&buf[pos], bufsize - pos, ", ");
        pos += limit_disp(&buf[pos], bufsize - pos, fp->r[vi-1].high);
        pos += snprintf(&buf[pos], bufsize - pos, "]");
        if (vi < nv)
            pos += snprintf(&buf[pos], bufsize - pos, "; ");
    }
    return pos;
}

uint frag_disp(char *buf, uint bufsize, fid_t fi) {
    return fragp_disp(buf, bufsize, frag_p(fi));
}

void frag_dump(fid_t fi) {
    char buf[frag_dumpsize()];
    frag_disp(buf, sizeof(buf), fi);
    fprintf(stderr, "%s\n", buf);
}

void done_frags(void) {
    for (uint fi = 0; fi < nfrags; ++fi)
        free_frag(frags[fi]);
    free(frags);
    free(limits);
}

void init_frags(void) {
    resize_limits(0);           /* malloc for fixed limits */
    memset(limit_lc(LIM0), 0, lc_size());
    limit_num_set(LIM0, 0);
    limit_den_set(LIM0, 1);
    memset(limit_lc(LIM1), 0, lc_size());
    lc_set(limit_lc(LIM1), 0, 1);
    limit_num_set(LIM1, 1);
    limit_den_set(LIM1, 1);

    fid_t fi = new_frag();
    frag_p(fi)->fid = fi;
    frag_p(fi)->parent = fi;
    frag_p(fi)->ps = all_paths();
    for (uint vi = 1; vi <= nv; ++vi) {
        frag_range(fi, vi)->low = LIM0;
        frag_range(fi, vi)->high = LIM1;
    }
}

int denorm_addmul(fid_t fi, int *c, int q, uint vmax, int dir) {
    int cmax = c[vmax];
    assert(cmax != 0);
    limitid_t lmax = ((dir < 0) ^ (cmax < 0))
        ? frag_range(fi, vmax)->low
        : frag_range(fi, vmax)->high;
    cmax *= limit_num(lmax);
    int qmax = limit_den(lmax);
    if (qmax > 1) {
        int g = igcd(cmax, qmax);
        cmax /= g;
        qmax /= g;
    }
    for (uint vi = 0; vi < vmax; ++vi)
        c[vi] = c[vi] * qmax + lc_get(limit_lc(lmax), vi) * cmax;
    return q * qmax;
}

void find_range_dir(crange_t *targ, fid_t fi, int *cs, int q, uint vmax, int dir) {
    int c[vmax + 1];
    memcpy(&c[0], &cs[0], (vmax + 1) * sizeof(int));
    while (vmax > 0) {
        if (c[vmax])
            q = denorm_addmul(fi, c, q, vmax, dir);
        --vmax;
    }
    int g = igcd(c[0], q);
    targ->num = c[0] / g;
    targ->den = q / g;
}

/* find const min and max for given expression according to limits of fi,
 * storing results in FRC[0, 1].
 */
void find_range(fid_t fi, int *c, int q, int vmax) {
    find_range_dir(&FRC[0], fi, c, q, vmax, -1);
    find_range_dir(&FRC[1], fi, c, q, vmax, 1);
}

/* Find a useful point to split frag fi pivoting on cs[]; split fi
 * to (fi, fn) and return fn.
 * We know that the pivot expression p takes values straddling 0.
 * If p = q(x + cy) for some variable y, we have min(x) + min(cy) < 0 and
 * max(x) + max(cy) > 0.
 * If we find that min(x) + max(cy) >= 0, and that max(x) + min(cy) <= 0,
 * then y = x/c is a valid splitting point. If not, we must recurse deeper
 * to constrain p - q max(cy) and/or p - q min(cy).
 * 
 * FIXME: we care about the sign of the ranges found, not absolute values -
 * consider denormalizing to lc * sign(num), and working only with that.
 * That may mean no more than discarding q, qmax.
 * FIXME: add overflow checks
 */
fid_t find_split(fid_t fi, int *cs, int q, uint vmax) {
    while (vmax && cs[vmax] == 0)
        --vmax;
    assert(vmax > 0);

    int c[vmax];

    /* p/q */
    memcpy(&c[0], &cs[0], (vmax + 1) * sizeof(int));
    /* (p/q - max(cy))/q_max */
    int qmax = denorm_addmul(fi, c, q, vmax, 1);
    find_range_dir(&FRC[0], fi, c, qmax, vmax - 1, -1);
    if (FRC[0].num < 0)
        return find_split(fi, c, qmax, vmax - 1);

    /* try r[vi] max */
    memcpy(&c[0], &cs[0], (vmax + 1) * sizeof(int));
    qmax = denorm_addmul(fi, c, q, vmax, -1);
    find_range_dir(&FRC[1], fi, c, qmax, vmax - 1, 1);
    if (FRC[1].num > 0)
        return find_split(fi, c, qmax, vmax - 1);

    /* we want to split vi at -li_{vi-1} / c */
    limitid_t ln = new_limit();
    limit_set_norm(ln, cs, -cs[vmax], vmax - 1);

    fid_t fn = frag_dup(fi);
    if (debug_split) {
        char buf1[limit_dumpsize()], buf2[limit_dumpsize()],
                buf3[limit_dumpsize()];
        limit_disp(buf1, sizeof(buf1), frag_range(fi, vmax)->low);
        limit_disp(buf2, sizeof(buf2), frag_range(fi, vmax)->high);
        limit_disp(buf3, sizeof(buf3), ln);
        fprintf(stderr, "split %c[%s, %s] at %s to %u, %u\n",
                'a' + vmax - 1, buf1, buf2, buf3, fi, fn);
    }
    assert(limit_cmp(ln, frag_range(fi, vmax)->low, vmax - 1) != 0);
    assert(limit_cmp(ln, frag_range(fi, vmax)->high, vmax - 1) != 0);
    frag_range(fi, vmax)->high = ln;
    frag_range(fn, vmax)->low = ln;
    return fn;
}

void split_one(fid_t fi, int *c, uint vmax, uint pi, uint pj) {
    find_range(fi, c, 1, vmax);
    int plow = FRC[0].num;
    int qlow = FRC[0].den;
    int phigh = FRC[1].num;
    int qhigh = FRC[1].den;
    assert(plow * qhigh < qlow * phigh);
    if (phigh <= 0) {
        frag_p(fi)->ps &= ~(1 << pj);
        return;
    }
    if (plow >= 0) {
        frag_p(fi)->ps &= ~(1 << pi);
        return;
    }
    /* pivot straddles zero, so find a split */
    fid_t fj = find_split(fi, c, 1, vmax);
    split_one(fi, c, vmax, pi, pj);
    split_one(fj, c, vmax, pi, pj);
}

static inline uint fls(uint x) {
    uint c = __builtin_clz(x);
    return 8 * sizeof(x) - c - 1;
}

void split_all_for(uint pi, uint pj) {
    int c[nv + 1];
    pathset_t ps = (1 << pi) | (1 << pj);
    path_t ppi = path_p(pi);
    path_t ppj = path_p(pj);
    uint vmax = fls(ppi ^ ppj) + 1;
    path_t common = ppi & ppj;
    ppi &= ~common;
    ppj &= ~common;
    c[0] = 0;
    for (uint i = 0; i < nv; ++i) {
        uint var = i + 1;
        uint bit = 1 << i;
        c[var] = (ppi & bit) ? 1 : (ppj & bit) ? -1 : 0;
    }

    for (fid_t fi = 0; fi < nfrags; ++fi) {
        if ((fi % 100) == 0)
            diag("split %u/%u", fi, nfrags);
        if (debug_split) {
            char buf[frag_dumpsize()];
            frag_disp(buf, sizeof(buf), fi);
            fprintf(stderr, "try split %s\n", buf);
        }
        if ((frag_p(fi)->ps & ps) != ps) {
            if (debug_split)
                fprintf(stderr, ".. does not match\n");
            continue;
        }
        split_one(fi, &c[0], vmax, pi, pj);
    }
    diag("");
}
