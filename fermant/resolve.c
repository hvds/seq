#include <stdlib.h>

#include "int.h"
#include "resolve.h"
#include "path.h"
#include "frag.h"
#include "num.h"

resolve_t *resolves = NULL;
uint nresolve = 0;
uint sizeresolve = 0;

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

void done_resolve(void) {
    free(resolves);
}

void init_resolve(void) {
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
        report("resolve %d: [%s] [%s]\n", ri, buf1, buf2);
    }
}

void split_all(void) {
    for (uint ri = 0; ri < nresolve; ++ri) {
        resolve_t *rp = &resolves[ri];
        split_all_for(rp->pi, rp->pj);
        report("split_all_for(%u, %u) gives %u frags\n",
                rp->pi, rp->pj, nfrags);
    }
    report("after split_all have %u frags\n", nfrags);
}
