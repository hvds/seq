#ifndef LIMIT_H
#define LIMIT_H

#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "int.h"
#include "lincom.h"
#include "num.h"

typedef lincomz_t limitz_t;
/* notionally:
    typedef struct {
        lincom_t lc;
        limitz_t num;
        limitz_t den;
    } limit_t;
*/
typedef void limit_t;

extern limit_t *LIM0;   /* /* initial integration limits _0^1 */
extern limit_t *LIM1;

extern void init_limits(void);
extern void done_limits(void);

__inline void *add_p(void *vp, size_t off) {
    return (void *)(((char *)vp) + off);
}

__inline uint limit_size(void) {
    return lc_size() + 2 * sizeof(limitz_t);
}

__inline lincom_t *limitp_lc(limit_t *lp) {
    return (lincom_t *)add_p(lp, 0);
}

__inline limitz_t limitp_num(limit_t *lp) {
    limitz_t *lzp = (limitz_t *)add_p(lp, lc_size());
    return lzp[0];
}

__inline limitz_t limitp_den(limit_t *lp) {
    limitz_t *lzp = (limitz_t *)add_p(lp, lc_size());
    return lzp[1];
}

__inline limitz_t limitp_num_set(limit_t *lp, limitz_t v) {
    limitz_t *lzp = (limitz_t *)add_p(lp, lc_size());
    lzp[0] = v;
}

__inline limitz_t limitp_den_set(limit_t *lp, limitz_t v) {
    limitz_t *lzp = (limitz_t *)add_p(lp, lc_size());
    lzp[1] = v;
}

__inline bool limitp_is_zero(limit_t *lp, uint vmax) {
    if (limitp_num(lp) == 0)
        return 1;
    return lc_is_zero(limitp_lc(lp), vmax);
}

__inline int limitp_cmp(limit_t *la, limit_t *lb, uint vmax) {
    int cmp = limitp_num(la) * limitp_den(lb) - limitp_num(lb) * limitp_den(la);
    return cmp ? cmp : lc_cmp(limitp_lc(la), limitp_lc(lb), vmax);
}

__inline void limitp_dup(limit_t *ld, limit_t *ls) {
    if (ld != ls)
        memcpy(ld, ls, limit_size());
}

__inline void limitp_set_norm(limit_t *lp, int *c, int q, uint vmax) {
    limitp_dup(lp, LIM0);
    int p = (q < 0) ? -1 : 1;
    q = abs(q);
    while (vmax && c[vmax] == 0)
        --vmax;
    int g = c[0];
    for (uint vi = 1; abs(g) != 1 && vi <= vmax; ++vi)
        g = igcd(g, c[vi]);
    for (uint vi = 0; vi <= vmax; ++vi)
        lc_set(limitp_lc(lp), vi, c[vi] / g);
    p *= g;
    g = igcd(p, q);
    limitp_num_set(lp, p / g);
    limitp_den_set(lp, q / g);
}

extern uint limit_dumpsize(void);
extern uint limitp_disp(char *buf, uint bufsize, limit_t *lp);
extern void limitp_dump(limit_t *lp);

#endif /* LIMIT_H */
