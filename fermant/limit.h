#ifndef LIMIT_H
#define LIMIT_H

#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "int.h"
#include "lincom.h"
#include "num.h"

typedef signed short limitz_t;
/* notionally:
    typedef struct {
        lincom_t lc;
        limitz_t num;
        limitz_t den;
    } limit_t;
*/
typedef void limit_t;
typedef enum {
    FRLOW = 0,  /* range for find_range() */
    FRHIGH,
    FSMIN,      /* range for find_split() */
    FSMAX,
    SPLIT,      /* target for split_one() */
    LIM0,       /* initial integration limits _0^1 */
    LIM1,
    ELIM_MAX
} limitid_t;

extern void *limits;
extern uint nlimits;
extern uint sizelimits;

__inline void *add_p(void *vp, uint off) {
    return (void *)(((char *)vp) + off);
}

__inline uint limit_size(void) {
    return lc_size() + 2 * sizeof(limitz_t);
}

__inline void resize_limits(uint extra) {
    if (nlimits + extra <= sizelimits)
        return;
    uint newsize = 3 * sizelimits / 2;
    while (nlimits + extra > newsize)
        newsize += 100;
    limits = realloc(limits, newsize * limit_size());
    sizelimits = newsize;
}

__inline limitid_t new_limit(void) {
    resize_limits(1);
    return (limitid_t)nlimits++;
}

__inline limit_t *limit_p(limitid_t l) {
    return (limit_t *)add_p(limits, l * limit_size());
}

__inline lincom_t *limit_lc(limitid_t l) {
    return (lincom_t *)add_p(limit_p(l), 0);
}

__inline limitz_t limit_num(limitid_t l) {
    limitz_t *lzp = (limitz_t *)add_p(limit_p(l), lc_size());
    return lzp[0];
}

__inline limitz_t limit_den(limitid_t l) {
    limitz_t *lzp = (limitz_t *)add_p(limit_p(l), lc_size());
    return lzp[1];
}

__inline limitz_t limit_num_set(limitid_t l, limitz_t v) {
    limitz_t *lzp = (limitz_t *)add_p(limit_p(l), lc_size());
    lzp[0] = v;
}

__inline limitz_t limit_den_set(limitid_t l, limitz_t v) {
    limitz_t *lzp = (limitz_t *)add_p(limit_p(l), lc_size());
    lzp[1] = v;
}

__inline void limit_set(limitid_t l, lincom_t *lc, int num, int den) {
    lc_copy(limit_lc(l), lc);
    limit_num_set(l, num);
    limit_den_set(l, den);
}

__inline int limit_cmp(limitid_t la, limitid_t lb, uint vmax) {
    int cmp = limit_num(la) * limit_den(lb) - limit_num(lb) * limit_den(la);
    return cmp ? cmp : lc_cmp(limit_lc(la), limit_lc(lb), vmax);
}

__inline void limit_dup(limitid_t ld, limitid_t ls) {
    if (ld != ls)
        limit_set(ld, limit_lc(ls), limit_num(ls), limit_den(ls));
}

__inline void limit_norm(limitid_t l, uint vmax) {
    int p = limit_num(l) * lc_norm(limit_lc(l), vmax);
    int q = limit_den(l);
    if (q < 0) {
        p = -p;
        q = -q;
    }
    int g = igcd(p, q);
    p /= g;
    q /= g;
    limit_num_set(l, p);
    limit_den_set(l, q);
}

__inline void limit_set_norm(limitid_t l, int *c, int q, uint vmax) {
    limit_dup(l, LIM0);
    int p = (q < 0) ? -1 : 1;
    q = abs(q);
    while (vmax && c[vmax] == 0)
        --vmax;
    int g = c[0];
    for (uint vi = 1; abs(g) != 1 && vi <= vmax; ++vi)
        g = igcd(g, c[vi]);
    for (uint vi = 0; vi <= vmax; ++vi)
        lc_set(limit_lc(l), vi, c[vi] / g);
    p *= g;
    g = igcd(p, q);
    limit_num_set(l, p / g);
    limit_den_set(l, q / g);
}

__inline void limit_addmul(limitid_t ld, limitid_t ls, int mult, uint vmax) {
    int p = mult * limit_num(ls) * limit_den(ld);
    int q = limit_den(ls) * limit_num(ld);
    if (q < 0) {
        p = -p;
        q = -q;
    }
    if (q > 1) {
        uint g = ugcd(abs(p), q);
        p /= g;
        q /= g;
    }
    for (uint vi = 0; vi <= vmax; ++vi)
        lc_set(limit_lc(ld), vi,
            q * lc_get(limit_lc(ld), vi)
            + p * lc_get(limit_lc(ls), vi)
        );
    limit_den_set(ld, q * limit_den(ld));
    limit_norm(ld, vmax);
}

extern uint limit_dumpsize(void);
extern uint limitp_disp(char *buf, uint bufsize, limit_t *lp);
extern uint limit_disp(char *buf, uint bufsize, limitid_t l);
extern void limit_dump(limitid_t l);

#endif /* LIMIT_H */
