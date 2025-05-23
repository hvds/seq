#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "lincom.h"

extern inline uint lc_size(void);
extern inline lincomz_t lc_get(lincom_t *lc, uint var);
extern inline void lc_set(lincom_t *lc, uint var, lincomz_t v);
extern inline bool lc_is_const(lincom_t *lc, uint vmax);
extern inline bool lc_is_zero(lincom_t *lc, uint vmax);
extern inline int lc_cmp(lincom_t *lca, lincom_t *lcb, uint vmax);
extern inline void lc_copy(lincom_t *lcdst, lincom_t *lcsrc);
extern inline uint lc_maxvar(lincom_t *lc);
extern inline int lc_norm(lincom_t *lc, uint vmax);

uint lc_dumpsize(void) {
    assert(sizeof(lincom_t) <= 2);
    return (1 + 5 + 1) * (nv + 1) + 2;
}

uint lc_disp(char *buf, uint bufsize, lincom_t *lc, uint vmax) {
    uint pos = 0;
    for (uint vi = 0; vi <= vmax; ++vi) {
        int c = lc_get(lc, vi);
        if (c == 0)
            continue;
        if (c < 0 || pos)
            pos += snprintf(buf + pos, bufsize - pos, (c < 0) ? "-" : "+");
        c = abs(c);
        if (c != 1 || vi == 0)
            pos += snprintf(buf + pos, bufsize - pos, "%d", c);
        if (vi)
            pos += snprintf(buf + pos, bufsize - pos, "%c", 'a' - 1 + vi);
    }
    return pos;
}

void lc_dump(lincom_t *lc) {
    char buf[lc_dumpsize()];
    lc_disp(buf, sizeof(buf), lc, nv);
    fprintf(stderr, "%s\n", buf);
}
