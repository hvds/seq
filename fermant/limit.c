#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define LIMIT_INLINE
#include "limit.h"
#undef LIMIT_INLINE

void *limits = NULL;
uint nlimits = (uint)ELIM_MAX;
uint sizelimits = 0;

extern inline void *add_p(void *vp, uint off);
extern inline uint limit_size(void);
extern inline void resize_limits(uint extra);
extern inline limitid_t new_limit(void);
extern inline limit_t *limit_p(limitid_t l);
extern inline lincom_t *limit_lc(limitid_t l);
extern inline limitz_t limit_num(limitid_t l);
extern inline limitz_t limit_den(limitid_t l);
extern inline limitz_t limit_num_set(limitid_t l, limitz_t v);
extern inline limitz_t limit_den_set(limitid_t l, limitz_t v);
extern inline void limit_set(limitid_t l, lincom_t *lc, int num, int den);
extern inline int limit_cmp(limitid_t la, limitid_t lb, uint vmax);
extern inline void limit_dup(limitid_t ld, limitid_t ls);
extern inline void limit_norm(limitid_t l, uint vmax);
extern inline void limit_set_norm(limitid_t l, int *c, int q, uint vmax);
extern inline void limit_addmul(limitid_t ld, limitid_t ls, int mult, uint vmax);

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

uint limit_disp(char *buf, uint bufsize, limitid_t l) {
    return limitp_disp(buf, bufsize, limit_p(l));
}

void limit_dump(limitid_t l) {
    char buf[limit_dumpsize()];
    limit_disp(buf, sizeof(buf), l);
    fprintf(stderr, "%s\n", buf);
}
