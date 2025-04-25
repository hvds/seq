#ifndef FRAG_H
#define FRAG_H

#include <assert.h>

#include "types.h"
#include "path.h"
#include "limit.h"

/* notionally:
    typedef struct {
        limit_t low;
        limit_t high;
    } range_t;
*/
typedef void range_t;

typedef struct {
    int num;
    int den;
} crange_t;

typedef uint fid_t;
typedef struct {
    pathset_t ps;   /* paths reachable */
    char r[0];      /* range_t r[nv] */
} frag_t;

extern crange_t FRC[2];
extern frag_t **frags;
extern uint nfrags;
extern uint sizefrags;
extern uint next_fragid;

extern void init_frags(void);
extern void done_frags(void);
extern void split_all_for(uint pi, uint pj);
extern uint frag_dumpsize(void);
extern uint frag_disp(char *buf, uint bufsize, fid_t l);
extern void frag_dump(fid_t l);

__inline uint range_size(void) {
    return limit_size() * 2;
}

__inline limit_t *range_low(range_t *rp) {
    return (limit_t *)add_p(rp, 0);
}

__inline limit_t *range_high(range_t *rp) {
    return (limit_t *)add_p(rp, limit_size());
}

__inline void range_low_set(range_t *rp, limit_t *lp) {
    limitp_dup(range_low(rp), lp);
}

__inline void range_high_set(range_t *rp, limit_t *lp) {
    limitp_dup(range_high(rp), lp);
}

__inline uint frag_size(void) {
    return sizeof(frag_t) + nv * range_size();
}

__inline frag_t *frag_p(fid_t fi) {
    return frags[fi];
}

__inline void resize_frags(uint extra) {
    if (nfrags + extra <= sizefrags)
        return;
    uint newsize = 3 * sizefrags / 2;
    while (nfrags + extra > newsize)
        newsize += 100;
    frags = realloc(frags, newsize * sizeof(frag_t *));
    sizefrags = newsize;
}

__inline fid_t new_frag(void) {
    fid_t fi = next_fragid++;
    frag_t *fp = calloc(1, frag_size());
    resize_frags(1);
    frags[nfrags++] = fp;
    return fi;
}

__inline void free_frag(frag_t *fp) {
    free(fp);
}

__inline fid_t frag_dup(fid_t fi) {
    fid_t fn = new_frag();
    memcpy(frag_p(fn), frag_p(fi), frag_size());
    return fn;
}

__inline pathset_t frag_ps(fid_t fi) {
    return frag_p(fi)->ps;
}

__inline void frag_ps_set(fid_t fi, pathset_t ps) {
    frag_p(fi)->ps = ps;
}

__inline range_t *frag_range(fid_t fi, uint vi) {
    assert(vi > 0);
    return (range_t *)add_p(frag_p(fi),
            sizeof(frag_t) + (vi - 1) * range_size());
}

#endif /* FRAG_H */
