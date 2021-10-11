#include <stdlib.h>
#include <string.h>

#include "sym.h"

extern int n;

/* We're not handling order-4 symmetries for now (rot90/rot270) */
/*
 *  sym_t sym_order[SYM_ORDER] = { xY, Xy, XY, yx, yX, Yx, YX };
 */
sym_t sym_order[SYM_ORDER] = { xY, Xy, XY, yx, YX };

/*
    Return true if the dimensions of a grid are transposed by this symmetry.
*/
static bool is_transpose(sym_t s) {
    return (s >= yx) ? 1 : 0;
}

/*
    Apply this symmetry to the supplied location.
*/
loc_t sym_transloc(sym_t s, span_t span, loc_t l) {
#define INVx(z) (span.min.x + span.max.x - (z))
#define INVy(z) (span.min.y + span.max.y - (z))
#define SWAPx(z) ((z) - span.min.x + span.min.y)
#define SWAPy(z) ((z) - span.min.y + span.min.x)
#define INVSWAPx(z) INVy(SWAPx(z))
#define INVSWAPy(z) INVx(SWAPy(z))
    switch (s) {
        case xy:
            return (loc_t){ l.x, l.y };
        case xY:
            return (loc_t){ l.x, INVy(l.y) };
        case Xy:
            return (loc_t){ INVx(l.x), l.y };
        case XY:
            return (loc_t){ INVx(l.x), INVy(l.y) };
        case yx:
            return (loc_t){ SWAPy(l.y), SWAPx(l.x) };
        case yX:
            return (loc_t){ SWAPy(l.y), INVSWAPx(l.x) };
        case Yx:
            return (loc_t){ INVSWAPy(l.y), SWAPx(l.x) };
        case YX:
            return (loc_t){ INVSWAPy(l.y), INVSWAPx(l.x) };
    }
}

/*
    Return the set of symmetries shown by this arrangement of points
*/
sym_t sym_check(loclist_t *ll, span_t span, int size) {
    sym_t s = 0;
    loc_t copy[n];

    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        int remain;

        if (is_transpose(si)
            && span.max.x - span.min.x != span.max.y - span.min.y
        )
            continue;

        memcpy(copy, ll->list, sizeof(loc_t) * size);
        remain = size;
        for (int j = 0; j < remain; ++j) {
            loc_t this = copy[j];
            loc_t that = sym_transloc(si, span, this);

            /* We only deal with order-2 symmetries here, so each point
             * is either its own pair, or has a partner that we match up
             */
            if (loc_eq(this, that))
                continue;
            for (int k = j + 1; k < remain; ++k) {
                if (loc_eq(copy[k], that)) {
                    copy[k] = copy[--remain];
                    goto this_ok;
                }
            }
            goto this_not_ok;
          this_ok:
            ;
        }
        s |= si;
      this_not_ok:
        ;
    }

    return s;
}

/*
    Return TRUE if this point is canonical under the known symmetries
*/
bool sym_best(sym_t s, span_t span, loc_t p1) {
    if (s == 0)
        return 1;
    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        if (s & si) {
            loc_t p2 = sym_transloc(si, span, p1);
            if (loc_lt(p2, p1))
                return 0;
        }
    }
    return 1;
}

/*
    Return TRUE if this pair of points is canonical under the known
    symmetries.
    We first order them so that p1 < p2, then (assuming p1' < p2') require
    that p1' < p1 || (p1' == p1 && p2' <= p2).
*/
bool sym_best2(sym_t s, span_t span, loc_t p1, loc_t p2) {
    loc_t t;

    if (s == 0)
        return 1;

    if (loc_lt(p2, p1)) {
        t = p1;
        p1 = p2;
        p2 = t;
    }
    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        if (s & si) {
            loc_t p3 = sym_transloc(si, span, p1);
            loc_t p4 = sym_transloc(si, span, p2);
            if (loc_lt(p3, p4)) {
                if (loc_lt(p3, p1) || (loc_eq(p3, p1) && loc_lt(p4, p2)))
                    return 0;
            } else {
                if (loc_lt(p4, p1) || (loc_eq(p4, p1) && loc_lt(p3, p2)))
                    return 0;
            }
        }
    }
    return 1;
}

/*
    Return TRUE if the edge formed by this pair of points lies on
    a known axis of symmetry.
*/
bool sym_axis(sym_t s, span_t span, loc_t p1, loc_t p2) {
    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        if ((s & si) && loc_eq(p1, sym_transloc(si, span, p1))
                && loc_eq(p2, sym_transloc(si, span, p2)))
            return 1;
    }
    return 0;
}
