#ifndef LINCOM_H
#define LINCOM_H

#include <string.h>

#include "types.h"
#include "int.h"
#include "num.h"

typedef signed short lincomz_t;
/* notionally:
    typedef struct {
        lincomz_t value[nv + 1];
    } lincom_t;
*/
typedef void lincom_t;

__inline uint lc_size(void) {
    return (nv + 1) * sizeof(lincomz_t);
}

#define LINCOM_ALLOC(name) \
    lincomz_t z ## name[lc_size()]; \
    lincom_t *name = &z ## name[0]

__inline lincomz_t lc_get(lincom_t *lc, uint var) {
    return ((lincomz_t *)lc)[var];
}

__inline void lc_set(lincom_t *lc, uint var, lincomz_t v) {
    ((lincomz_t *)lc)[var] = v;
}

__inline int lc_is_const(lincom_t *lc, uint vmax) {
    for (uint vi = 1; vi <= vmax; ++vi)
        if (lc_get(lc, vi))
            return 0;
    return 1;
}

__inline int lc_cmp(lincom_t *lca, lincom_t *lcb, uint vmax) {
    return memcmp(lca, lcb, (vmax + 1) * sizeof(lincomz_t));
}

__inline void lc_copy(lincom_t *lcdst, lincom_t *lcsrc) {
    memcpy(lcdst, lcsrc, lc_size());
}

__inline uint lc_maxvar(lincom_t *lc) {
    for (uint i = nv; i > 0; --i)
        if (lc_get(lc, i))
            return i;
    return 0;
}

/* divide through by gcd, and return it; set sign such that most significant
 * variable is positive */
__inline int lc_norm(lincom_t *lc, uint vmax) {
    int g = lc_get(lc, vmax);
    for (uint vi = vmax; vi > 0; ) {
        --vi;
        g = (g == 0) ? igcd(g, lc_get(lc, vi)) : igcd(lc_get(lc, vi), g);
        if (g == 1)
            return 1;
    }
    if (g != 1 && g != 0)
        for (uint vi = 0; vi <= vmax; ++vi)
            lc_set(lc, vi, lc_get(lc, vi) / g);
    return g;
}

extern uint lc_dumpsize(void);
extern uint lc_disp(char *buf, uint bufsize, lincom_t *lc, uint vi);
extern void lc_dump(lincom_t *lc);

#endif  /* LINCOM_H */
