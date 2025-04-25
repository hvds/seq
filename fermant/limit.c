#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define LIMIT_INLINE
#include "limit.h"
#undef LIMIT_INLINE

limit_t *LIM0;
limit_t *LIM1;

extern inline void *add_p(void *vp, uint off);
extern inline uint limit_size(void);
extern inline lincom_t *limitp_lc(limit_t *lp);
extern inline limitz_t limitp_num(limit_t *lp);
extern inline limitz_t limitp_den(limit_t *lp);
extern inline limitz_t limitp_num_set(limit_t *lp, limitz_t v);
extern inline limitz_t limitp_den_set(limit_t *lp, limitz_t v);
extern inline int limitp_cmp(limit_t *la, limit_t *lb, uint vmax);
extern inline void limitp_dup(limit_t *ld, limit_t *ls);
extern inline void limitp_set_norm(limit_t *lp, int *c, int q, uint vmax);

void init_limits(void) {
    LIM0 = calloc(1, limit_size());
    LIM1 = calloc(1, limit_size());
    limitp_den_set(LIM0, 1);
    lc_set(limitp_lc(LIM1), 0, 1);
    limitp_num_set(LIM1, 1);
    limitp_den_set(LIM1, 1);
}

void done_limits(void) {
    free(LIM0);
    free(LIM1);
}

uint limit_dumpsize(void) {
    /* "-123(lincom)/456\n" */
    assert(sizeof(limitz_t) <= 2);
    return 7 + lc_dumpsize() + 9;
}

uint limitp_disp(char *buf, uint bufsize, limit_t *lp) {
    lincom_t *lc = (lincom_t *)add_p(lp, 0);
    limitz_t *lzp = (limitz_t *)add_p(lp, lc_size());
    int p = lzp[0];
    int q = lzp[1];
    int terms = 0;
    for (uint i = 0; i <= nv; ++i)
        if (lc_get(lc, i))
            ++terms;

    uint pos = 0;
    if (terms == 0) {
        pos += snprintf(&buf[pos], bufsize - pos, "0");
    } else if (terms == 1) {
        p *= lc_get(lc, lc_maxvar(lc));
        pos += snprintf(&buf[pos], bufsize - pos, "%d", p);
        if (q != 1)
            pos += snprintf(&buf[pos], bufsize - pos, "/%d", q);
    } else if (p == 1 && q == 1)
        pos += lc_disp(&buf[pos], bufsize - pos, lc, nv);
    else {
        if (p == -1)
            pos += snprintf(&buf[pos], bufsize - pos, "-");
        else if (p != 1)
            pos += snprintf(&buf[pos], bufsize - pos, "%d", p);
        pos += snprintf(&buf[pos], bufsize - pos, "(");
        pos += lc_disp(&buf[pos], bufsize - pos, lc, nv);
        pos += snprintf(&buf[pos], bufsize - pos, ")");
        if (q != 1)
            pos += snprintf(&buf[pos], bufsize - pos, "/%d", q);
    }
    return pos;
}

void limitp_dump(limit_t *lp) {
    char buf[limit_dumpsize()];
    limitp_disp(buf, sizeof(buf), lp);
    fprintf(stderr, "%s\n", buf);
}
