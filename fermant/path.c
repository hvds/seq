#include <stdlib.h>
#include <strings.h>

#include "int.h"
#include "path.h"
#include "num.h"
#include "frag.h"
#include "source.h"
#include "../inc/diag.h"

path_t *paths = NULL;
uint npaths = 0;
uint sizepaths = 0;
resolve_t *resolves = NULL;
uint nresolve = 0;
uint sizeresolve = 0;

extern inline path_t path_p(uint pi);
extern inline path_t include_line(path_t cur, uint var);
extern inline uint path_first(path_t p);
extern inline path_t path_xor(uint pi, uint pj);
extern inline uint path_count(path_t p);
extern inline path_t pathset_first(pathset_t ps);
extern inline uint pathset_count(pathset_t ps);

/* Lookup table (used only for the initialization of paths) assigning
 * a variable (0 .. nv-1) to each line segment.
 * The first (na-1)nb entries are for eastward line segments; the remaining
 * na(nb-1) are northward ones.
 */
uint *transvar;

/* Returns the (zero-based) index of the variable representing the line
 * segment that travels in direction (0 = east, 1 = north) from the node
 * at (x, y).
 */
uint translate_var(uint x, uint y, uint direction) {
    uint offset = (direction == 0)
        ? (na - 1) * y + x
        : ((na - 1) * nb) + na * y + x;
    return transvar[offset];
}

/* Initializes the transvar lookup table according to the specified
 * strategy.
 */
void init_variable_allocation(int strategy) {
    transvar = calloc(nv, sizeof(uint));

    /* TODO: add more strategies - plain, nearest-sw, nearest-se */
    /* strategy boustrophedon: easting left to right, northing right to left
     *  *k*l*
     *  j i h
     *  *f*g*
     *  e d c
     *  *a*b*
     */
    uint ix = 0;
    for (uint y = 0; y < nb; ++y)
        for (uint x = 0; x < na - 1; ++x)
            transvar[ix++] = y * (2 * na - 1) + x;
    for (uint y = 0; y < nb - 1; ++y)
        for (uint x = 0; x < na; ++x)
            transvar[ix++] = (y + 1) * (2 * na - 1) - (x + 1);
    return;
}

/* Sort a+b < a+c < a < b+c */
int path_comparator(const void *a, const void *b) {
    path_t pa = *(path_t *)a;
    path_t pb = *(path_t *)b;
    path_t diff = pa ^ pb;
    if (diff == 0)
        return 0;
    path_t mindiff = 1 << path_first(diff);
    return (pa & mindiff) ? -1 : 1;
}

int resolve_comparator(const void *a, const void *b) {
    resolve_t *ra = (resolve_t *)a;
    resolve_t *rb = (resolve_t *)b;
    path_t pa = path_xor(ra->pi, ra->pj);
    path_t pb = path_xor(rb->pi, rb->pj);
    uint ca = path_count(pa);
    uint cb = path_count(pb);
    if (ca != cb)
        return (int)ca - (int)cb;
    return -path_comparator(&pa, &pb);
}

/* Writes a zero-terminated text representation of the path with the specified
 * index into the provided buffer. Returns the length written.
 */
uint render_path(char *buf, uint buflen, uint pi) {
    uint off = 0;
    if (pi > npaths)
        fail("Invalid path index %d, max known is %d", pi, npaths);
    path_t p = paths[pi];
    while (p && off < buflen) {
        uint index = path_first(p);
        p &= ~(1 << index);
        if (off < buflen)
            buf[off++] = 'a' + index;
        if (off < buflen)
            buf[off++] = p ? '+' : 0;
    }
    if (off == 0 && off < buflen)
        buf[off++] = 0;
    return off;
}

/* Recursive helper to generate paths: generate additional paths given
 * we are at node (x, y) and took the path cur to reach this point.
 * Sets paths[] and npaths.
 */
void gen_paths_r(path_t cur, uint x, uint y) {
    if (x >= na - 1) {
        if (y >= nb - 1) {
            if (npaths >= sizepaths)
                fail("panic: index %d exceeds expected path count %d",
                        npaths, sizepaths);
            paths[npaths++] = cur;
            return;
        }
        goto try_north;
    }
    /* try east */
    gen_paths_r(include_line(cur, translate_var(x, y, 0)), x + 1, y);
    if (y < nb - 1) {
      try_north:
        gen_paths_r(include_line(cur, translate_var(x, y, 1)), x, y + 1);
    }
}

void done_paths(void) {
    free(paths);
    free(transvar);
    free(resolves);
}

/* Generate all paths, using the variable allocation defined by the
 * specified strategy.
 */
void init_paths(int strategy) {
    /* path_t must fit nv bits */
    if (sizeof(path_t) * 8 < nv)
        fail("panic: nv %d > sizeof(path_t) %d", nv, sizeof(path_t) * 8);

    /* there will be comb(a+b-2, a-1) distinct paths */
    sizepaths = icomb(na + nb - 2, na - 1);
    if (sizeof(pathset_t) * 8 < sizepaths)
        fail("panic: npaths %d > sizeof(pathset_t) %d",
                sizepaths, sizeof(pathset_t) * 8);
    paths = calloc(sizepaths, sizeof(path_t));

    init_variable_allocation(strategy);

    gen_paths_r((path_t)0, 0, 0);
    qsort(paths, npaths, sizeof(path_t), &path_comparator);

    char buf[nv * 2];
    for (uint pi = 0; pi < npaths; ++pi) {
        render_path(buf, sizeof(buf), pi);
        report("000 path %d: %s (%0x)\n", pi, buf, paths[pi]);
    }

    sizeresolve = npaths * (npaths - 1) / 2;
    resolves = calloc(sizeresolve, sizeof(resolve_t));
    uint ri = 0;
    for (uint pi = 0; pi < npaths; ++pi) {
        for (uint pj = pi + 1; pj < npaths; ++pj) {
            if (ri >= sizeresolve)
                fail("panic: resolves overflow");
            resolves[ri].pi = pi;
            resolves[ri].pj = pj;
            ++ri;
        }
    }
    nresolve = ri;
    qsort(resolves, nresolve, sizeof(resolve_t), &resolve_comparator);

    char buf1[2 * nv], buf2[2 * nv];
    for (uint ri = 0; ri < nresolve; ++ri) {
        resolve_t *r = &resolves[ri];
        render_path(buf1, sizeof(buf1), r->pi);
        render_path(buf2, sizeof(buf2), r->pj);
        report("000 resolve %d: [%s] [%s]\n", ri, buf1, buf2);
    }
}

pathset_t all_paths(void) {
    return (pathset_t)((1 << npaths) - 1);
}

uint split_all(void) {
    uint count;
    for (uint ri = 0; ri < nresolve; ++ri) {
        resolve_t *rp = &resolves[ri];
        if (!resolve_open(ri))
            continue;
        count = split_all_for(rp->pi, rp->pj);
        resolve_close(ri);
        report("000 split_all_for(%u, %u) gives %u frags\n",
                rp->pi, rp->pj, count);
    }
    return count;
}

