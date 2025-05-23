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
#include "source.h"
#include "diag.h"

extern inline void reset_frags(void);
extern inline uint range_size(void);
extern inline limit_t *range_low(range_t *rp);
extern inline limit_t *range_high(range_t *rp);
extern inline void range_low_set(range_t *rp, limit_t *li);
extern inline void range_high_set(range_t *rp, limit_t *li);
extern inline uint frag_size(void);
extern inline frag_t *frag_p(fid_t f);
extern inline void resize_frags(uint extra);
extern inline fid_t new_frag(void);
extern inline fid_t frag_dup(fid_t fi);
extern inline fid_t frag_id(fid_t fi);
extern inline fid_t frag_parent(fid_t fi);
extern inline pathset_t frag_ps(fid_t fi);
extern inline void frag_ps_set(fid_t fi, pathset_t ps);
extern inline range_t *frag_range(fid_t f, uint vi);

#ifdef DEBUG
#   define dassert(x) assert(x)
#else
#   define dassert(x)
#endif

crange_t FRC[2];
frag_t *frags = NULL;
uint nfrags = 0;
uint sizefrags = 0;
fid_t nextfragid = 0;

uint frag_dumpsize(void) {
    /* "frag 999(998) 0x18: [0, a+b]; [...]...\n" */
    dassert(sizeof(fid_t) <= 4);
    return 5+22+25 + (npaths + 3)/4 + (limit_dumpsize() * 2 + 6) * nv + 2;
}

uint frag_disp(char *buf, uint bufsize, fid_t fi) {
    uint pos = 0;
    pos += snprintf(&buf[pos], bufsize - pos, "frag %d(%d) 0x%x: ",
            (int)frag_id(fi), (int)frag_parent(fi), frag_ps(fi));
    for (uint vi = 1; vi <= nv; ++vi) {
        pos += snprintf(&buf[pos], bufsize - pos, "[");
        pos += limitp_disp(&buf[pos], bufsize - pos,
                range_low(frag_range(fi, vi)));
        pos += snprintf(&buf[pos], bufsize - pos, ", ");
        pos += limitp_disp(&buf[pos], bufsize - pos,
                range_high(frag_range(fi, vi)));
        pos += snprintf(&buf[pos], bufsize - pos, "]");
        if (vi < nv)
            pos += snprintf(&buf[pos], bufsize - pos, "; ");
    }
    return pos;
}

void frag_dump(fid_t fi) {
    char buf[frag_dumpsize()];
    frag_disp(buf, sizeof(buf), fi);
    fprintf(stderr, "%s\n", buf);
}

void done_frags(void) {
    free(frags);
    done_limits();
}

void init_frags(void) {
    init_limits();
}

int denorm_addmul(fid_t fi, int *c, int q, uint vmax, int dir) {
    int cmax = c[vmax];
    dassert(cmax != 0);
    limit_t *lmax = ((dir < 0) ^ (cmax < 0))
        ? range_low(frag_range(fi, vmax))
        : range_high(frag_range(fi, vmax));
    cmax *= limitp_num(lmax);
    int qmax = limitp_den(lmax);
    if (qmax > 1) {
        int g = igcd(cmax, qmax);
        cmax /= g;
        qmax /= g;
    }
    for (uint vi = 0; vi < vmax; ++vi)
        c[vi] = c[vi] * qmax + lc_get(limitp_lc(lmax), vi) * cmax;
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
    dassert(vmax > 0);

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
    char vln[limit_size()];
    limit_t *ln = (limit_t *)&vln[0];
    limitp_set_norm(ln, cs, -cs[vmax], vmax - 1);

    fid_t fn = frag_dup(fi);
#ifdef DEBUG
    if (debug_split) {
        char buf1[limit_dumpsize()], buf2[limit_dumpsize()],
                buf3[limit_dumpsize()];
        limitp_disp(buf1, sizeof(buf1), range_low(frag_range(fi, vmax)));
        limitp_disp(buf2, sizeof(buf2), range_high(frag_range(fi, vmax)));
        limitp_disp(buf3, sizeof(buf3), ln);
        fprintf(stderr, "split %c[%s, %s] at %s to %u, %u\n",
                'a' + vmax - 1, buf1, buf2, buf3, frag_id(fi), frag_id(fn));
    }
#endif
    dassert(limitp_cmp(ln, range_low(frag_range(fi, vmax)), vmax - 1) != 0);
    dassert(limitp_cmp(ln, range_high(frag_range(fi, vmax)), vmax - 1) != 0);
    range_low_set(frag_range(fi, vmax), ln);
    range_high_set(frag_range(fn, vmax), ln);
    return fn;
}

fid_t base_fi;

void split_one(fid_t fi, int *c, uint vmax, uint pi, uint pj) {
    if (need_diag) {
        diag("split %u +%u", base_fi, nfrags);
        need_diag = 0;
    }

    find_range(fi, c, 1, vmax);
    int plow = FRC[0].num;
    int qlow = FRC[0].den;
    int phigh = FRC[1].num;
    int qhigh = FRC[1].den;
    dassert(plow * qhigh < qlow * phigh);
    if (phigh <= 0)
        return write_frag(fi, 1 << pj);
    if (plow >= 0)
        return write_frag(fi, 1 << pi);
    /* pivot straddles zero, so find a split */
    fid_t fj = find_split(fi, c, 1, vmax);
    split_one(fi, c, vmax, pi, pj);
    split_one(fj, c, vmax, pi, pj);
}

static inline uint fls(uint x) {
    uint c = __builtin_clz(x);
    return 8 * sizeof(x) - c - 1;
}

uint split_all_for(uint pi, uint pj) {
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

    uint count = 0;
    base_fi = 0;
    while (read_frag(new_frag())) {
        fid_t fi = nfrags - 1;
        ++base_fi;
        if (need_diag) {
            diag("split %u", base_fi);
            need_diag = 0;
        }
#ifdef DEBUG
        if (debug_split) {
            char buf[frag_dumpsize()];
            frag_disp(buf, sizeof(buf), fi);
            fprintf(stderr, "try split %s\n", buf);
        }
#endif
        if ((frag_ps(fi) & ps) == ps) {
            split_one(fi, &c[0], vmax, pi, pj);
            count += nfrags;
        } else {
#ifdef DEBUG
            if (debug_split)
                fprintf(stderr, ".. does not match\n");
#endif
            write_frag(fi, 0);
            count += 1;
        }
        reset_frags();
        if (need_log) {
            resolve_checkpoint();
            need_log = 0;
        }
    }
    return count;
}
