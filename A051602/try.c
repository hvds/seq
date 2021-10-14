#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/times.h>
#include <unistd.h>

#include "loc.h"
#include "sym.h"

long clock_tick;

/* A structure of context used at each level of recursion, which
 * encapsulates information about the arrangement up to this point,
 * and information about what extensions have already been tried.
 * It does not list the points themselves, those are stored globally
 * in point[].
 * At all times the 
 */
typedef struct {
    int squares;    /* Number of squares in this arrangement */
    int power;      /* Original unit square is now 2^power */
    span_t span;    /* Min/max extent of this arrangement */
    int try3[3];    /* The next triple to try (list of 3 point indices) */
    int try2[3];    /* The next pair to try, plus direction */
    loclist_t *seen;/* List of points tried, that we don't need to try again */
    pairlist_t *pairs;  /* List of pairs tried */
} cx_t;

int n;              /* We're trying to find A051602(n) */
int verbose = 0;    /* report every maximum (1) or iteration (2) */
int quiet = 0;      /* skip reporting for arrangements with < n points */
int best;           /* Greatest number of squares seen in any arrangement */
loc2plist_t *point; /* List of points in the current arrangement */
cx_t *context;      /* List of context objects */
loc2p_t minspan;    /* { x, y } size of smallest maximal solution */
loc2p_t maxspan;    /* { x, y } size of greatest maximal solution */
unsigned long visit;        /* Count of iterations */
unsigned long lim_visit = 0;/* Stop after this many iterations */
unsigned long last_new;     /* Iteration at which last new result found */

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

/* Show a result-so-far, consisting of i points */
void report(int i) {
    cx_t *cx = &context[i];
    int power = cx->power;
    loc_t span = loc_diff(cx->span.min, cx->span.max);

    printf("(%lu) p%d %dx%d ", visit, power,
            (span.x >> 1) + 1, (span.y >> 1) + 1);
    printf("%d:", cx->squares);
    for (int j = 0; j < i; ++j) {
        loc_t p = list2p_get(point, j, cx->power);
        printf(" %d:%d", p.x - cx->span.min.x, p.y - cx->span.min.y);
    }
    printf("\n");
}

void init(void) {
    loclist_t *first_seen = new_loclist(10);
    pairlist_t *first_pairs = new_pairlist(10);
    context = (cx_t *)malloc((n + 1) * sizeof(cx_t));

    /* used to calculate timings */
    clock_tick = sysconf(_SC_CLK_TCK);

    /* for less-verbose output don't buffer, even if stdout is not a tty */
    if (verbose != 2)
        setvbuf(stdout, (char*)NULL, _IONBF, 0);

    /* Start off with 4 points making a unit square */
    point = new_loc2plist(n + 1);
    list2p_set(point, 0, (loc2p_t){ 0, 0, 1 });
    list2p_set(point, 1, (loc2p_t){ 0, 2, 1 });
    list2p_set(point, 2, (loc2p_t){ 2, 0, 1 });
    list2p_set(point, 3, (loc2p_t){ 2, 2, 1 });

    context[4] = (cx_t){
        1,                      /* squares */
        1,                      /* power */
        { { 0, 0 }, { 2, 2 } }, /* span */
        /* We know that no three of our first four points will form an
         * empty triple, so we'll start looking for triples only from
         * the next point.
         */
        { 4, 1, 0 },    /* Next triple to look for is point[4] + [1] + [0] */
        { 1, 0, 0 },    /* Next pair to look for is point[1] + [0], in the
                         * first direction (of 2) */
        first_seen,     /* Seed with empty lists */
        first_pairs
    };

    best = 1;
    minspan = (loc2p_t){ 2, 2, 1 };
    maxspan = (loc2p_t){ 2, 2, 1 };
    visit = 0UL;
    last_new = 0UL;
}

void finish(void) {
    free_loc2plist(point);
    free_loclist(context[4].seen);
    free_pairlist(context[4].pairs);
    free(context);
}

/* We have a double recursion: try_next() finds the next extension to try,
 * and calls try_with() which applies that extension and calculates the
 * effects before calling try_next() again.
 */
void try_next(int depth, int new);

/* Attempt to find a given 2-point extension; on failure, returns FALSE;
 * on success, returns TRUE with the 2 new points appended to point[].
 * For each pair of candidate points, we can extend in 2 directions;
 * considering the candidate as an edge of a square, we can make the
 * second edge by rotating it either 90 or 270 degrees.
 */
int try_test2(int points, int try2[3], int power) {
    loc_t pi = list2p_get(point, try2[0], power);
    loc_t pj = list2p_get(point, try2[1], power);
    int dir = try2[2];
    loc_t diff;
    loc2p_t pk, pl;

    diff = loc_diff(pi, pj);
    if (dir == 1) {
        pk = (loc2p_t){ loc_rot90(pi, diff), power };
        pl = (loc2p_t){ loc_rot90(pj, diff), power };
    } else if (dir == 0) {
        pk = (loc2p_t){ loc_rot270(pi, diff), power };
        pl = (loc2p_t){ loc_rot270(pj, diff), power };
    } else {
        pk = (loc2p_t){ loc_diag1(pi, pj), power };
        pl = (loc2p_t){ loc_diag2(pi, pj), power };
    }
    /* This is a valid extension only if the two new points are not
     * already present in the arrangement.
     */
    if (list2p_exists_lim(point, points, pk)
        || list2p_exists_lim(point, points, pl)
    )
        return 0;
    list2p_set(point, points, pk);
    list2p_set(point, points + 1, pl);
    return 1;
}

/* Return TRUE if the new pair just found by try_test2 would be
 * a duplicate of a pair we've already tried.
 */
int duplicate_pair(int points, sym_t sym, int try2[3], int power) {
    int ii = try2[0], ij = try2[1], ik = points, il = points + 1;
    loc_t pi = list2p_get(point, ii, power);
    loc_t pj = list2p_get(point, ij, power);
    loc_t pk = list2p_get(point, ik, power);
    loc_t pl = list2p_get(point, il, power);
    loc2p_t pm = (loc2p_t){ pk.x * 2 - pi.x, pk.y * 2 - pi.y, power };
    loc2p_t pn = (loc2p_t){ pl.x * 2 - pj.x, pl.y * 2 - pj.y, power };
    int im, in;

    /* if the new pair lies on an axis of symmetry, duplication has already
     * been catered for
     */
    if (sym_axis(sym, context[points].span, pk, pl))
        return 0;

    /* if the opposite pair is not already in the arrangement, there's
     * no duplication
     */
    if ((im = list2p_find_lim(point, points, pm)) < 0)
        return 0;
    if ((in = list2p_find_lim(point, points, pn)) < 0)
        return 0;

    /* if the opposite pair comes earlier in the order than the original
     * pair, it's a duplicate, else it's not
     */
    return ((im > in)
        ? (im < ii || (im == ii && in < ij))
        : (in < ii || (in == ii && im < ij))
    );
}

/* Attempt to find a given 1-point extension; on failure, returns FALSE;
 * on success, returns TRUE with the new points appended to point[].
 * We need to determine a way of joining the 3 input points as two edges
 * so that we can extrapolate a fourth point making a square; we do this
 * by attempting to rotate one edge by 90 or 270 degrees to make it
 * coincide with the third point.
 */
int try_test3(int points, int try3[3], int power) {
    loc_t pi = list2p_get(point, try3[0], power);
    loc_t pj = list2p_get(point, try3[1], power);
    loc_t pk = list2p_get(point, try3[2], power);
    loc_t pl, diff;

    diff = loc_diff(pi, pj);
    if (loc_eq(pk, loc_rot90(pi, diff))) {
        pl = loc_rot90(pj, diff);
    } else if (loc_eq(pk, loc_rot90(pj, diff))) {
        pl = loc_rot90(pi, diff);
    } else if (loc_eq(pk, loc_rot270(pi, diff))) {
        pl = loc_rot270(pj, diff);
    } else if (loc_eq(pk, loc_rot270(pj, diff))) {
        pl = loc_rot270(pi, diff);
    } else {
        diff = loc_diff(pi, pk);
        if (loc_eq(pj, loc_rot90(pk, diff))) {
            pl = loc_rot90(pi, diff);
        } else if (loc_eq(pj, loc_rot270(pk, diff))) {
            pl = loc_rot270(pi, diff);
        } else {
            return 0;
        }
    }
    /* This is a valid extension only if the new point is not already
     * present in the arrangement.
     */
    if (list2p_exists_lim(point, points, (loc2p_t){ pl, power }))
        return 0;
    list2p_set(point, points, (loc2p_t){ pl, power });
    return 1;
}

/* Count the new squares formed by adding point[i] to the previous
 * arrangement. To avoid double-counting, we look at the two edges
 * point[i]-point[j] and point[i]-point[k], and count only when
 * j < k.
 */
int find_squares(int i, int power) {
    loc_t pi = list2p_get(point, i, power);
    int count = 0;

    for (int j = 0; j < i; ++j) {
        loc_t pj = list2p_get(point, j, power);
        loc_t diff = loc_diff(pi, pj);
        loc_t pk = loc_rot90(pi, diff);
        loc_t pl = loc_rot90(pj, diff);
        if (list2p_find_lim(point, i, (loc2p_t){ pk, power }) > j
           && list2p_exists_lim(point, i, (loc2p_t){ pl, power })
        )
            ++count;

        pk = loc_rot270(pi, diff);
        pl = loc_rot270(pj, diff);
        if (list2p_find_lim(point, i, (loc2p_t){ pk, power }) > j
           && list2p_exists_lim(point, i, (loc2p_t){ pl, power })
        )
            ++count;
    }
    return count;
}

/* Extend the existing arrangement of 'points' points with 'new'
 * additional points (which have already been placed in point[]).
 * We update context[points+new].squares and .span, and check
 * whether we're reached a new maximum, or the same maximum with
 * a new grid size.
 * This finishes with a recursive call to try_next(); on return,
 * all arrangements with this extension will have been checked.
 */
void try_with(int points, int new) {
    cx_t *ocx = &context[points];
    cx_t *ncx = &context[points + new];
    int power = ocx->power, raising = 0;

    if (is_loc2p_odd(list2p_get2p(point, points), power)) {
        raising = 1;
        ++power;
    }

    /* The caller will already have filled in most of the new context,
     * but some parts may have been overwritten by subsequent work.
     */
    ncx->power = power;
    ncx->span = ocx->span;
    if (raising) {
        ncx->span.min = raise2p(ncx->span.min);
        ncx->span.max = raise2p(ncx->span.max);
    }

    for (int i = points; i < points + new; ++i) {
        loc_t pi = list2p_get(point, i, power);

        /* Update the count of squares */
        ncx->squares += find_squares(i, power);

        /* Track the 4 limits of the arrangement */
        if (ncx->span.min.x > pi.x) ncx->span.min.x = pi.x;
        if (ncx->span.min.y > pi.y) ncx->span.min.y = pi.y;
        if (ncx->span.max.x < pi.x) ncx->span.max.x = pi.x;
        if (ncx->span.max.y < pi.y) ncx->span.max.y = pi.y;
    }

    if ((quiet == 0 || points + new == n) && ncx->squares >= best) {
        loc_t span = loc_diff(ncx->span.min, ncx->span.max);
        int newspan = 0;

        /* minspan/maxspan are canonicalized to call the larger dimension 'x' */
        if (span.x < span.y) {
            int tmp = span.x;
            span.x = span.y;
            span.y = tmp;
        }
        if (ncx->squares > best) {
            /* New record: reset minspan/maxspan, and always report this */
            minspan = (loc2p_t){ span, power };
            maxspan = (loc2p_t){ span, power };
            newspan = 1;
            last_new = visit;
        } else {
            /* Match of existing record, report it only if new gridsize */
            int spower = (minspan.power > maxspan.power) ? minspan.power : maxspan.power;
            if (spower < power) {
                minspan = lift2p(minspan, power);
                maxspan = lift2p(maxspan, power);
            } else if (power < spower) {
                span = convert2p((loc2p_t){ span, power }, spower);
            }

            if (minspan.p.x > span.x) minspan.p.x = span.x, newspan = 1;
            if (minspan.p.y > span.y) minspan.p.y = span.y, newspan = 1;
            if (maxspan.p.x < span.x) maxspan.p.x = span.x, newspan = 1;
            if (maxspan.p.y < span.y) maxspan.p.y = span.y, newspan = 1;
        }
        best = ncx->squares;
        if (newspan || verbose == 1) {
            if (verbose == 2) {
                /* distinguish from the per-iteration report */
                printf("* ");
            }
            report(points + new);
        }
    }

    try_next(points + new, new);
}

/* Mark a point as seen, either while applying it (present=true) or after
 * returning from applying it (present=false).
 * If present, we remove it from the seen list; if not, we add it.
 * Similarly, any pair including this point is removed from the pairs list,
 * and if present we promote the other point of the pair to the seen list.
 */
void seen_point(loclist_t *seen, pairlist_t *pairs, loc_t p, int power, int present) {
    if (present)
        list_remove(seen, p);
    else
        list_append(seen, p);

    int used = pairs->used;
    for (int i = used - 1; i >= 0; --i) {
        pair_t pair = pair_get(pairs, i);
        if (loc_eq(p, pair.p[0])) {
            if (present)
                list_append(seen, pair.p[1]);
            if (used--)
                pair_set(pairs, i, pair_get(pairs, used));
        } else if (loc_eq(p, pair.p[1])) {
            if (present)
                list_append(seen, pair.p[0]);
            if (used--)
                pair_set(pairs, i, pair_get(pairs, used));
        }
    }
    pairs->used = used;
}

/* Return TRUE if this pair of points is canonical under the known
 * symmetries.
 */
bool sym_best2(int points, sym_t s, span_t span, int power, int i1, int i2) {
    loc_t p1 = list2p_get(point, i1, power);
    loc_t p2 = list2p_get(point, i2, power);

    if (s == 0)
        return 1;

    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        if (s & si) {
            loc_t p3 = sym_transloc(si, span, p1);
            loc_t p4 = sym_transloc(si, span, p2);
            int i3 = list2p_find_lim(point, points, (loc2p_t){ p3, power });
            int i4 = list2p_find_lim(point, points, (loc2p_t){ p4, power });
            if (i3 < i4) {
                if (i4 < i1 || (i4 == i1 && i3 < i2))
                    return 0;
            } else {
                if (i3 < i1 || (i3 == i1 && i4 < i2))
                    return 0;
            }
        }
    }
    return 1;
}

/* Return TRUE if this triplet of points is canonical under the known
 * symmetries: ie if its highest-indexed point is less than the highest
 * indexed point of a given transform, or it is equal and the next-highest
 * is less, etc. The input indices are given in descending order.
 */
bool sym_best3(
    int points, sym_t s, span_t span, int power, int i1, int i2, int i3
) {
    loc_t p1 = list2p_get(point, i1, power);
    loc_t p2 = list2p_get(point, i2, power);
    loc_t p3 = list2p_get(point, i3, power);

    if (s == 0)
        return 1;

    for (int i = 0; i < SYM_ORDER; ++i) {
        sym_t si = sym_order[i];
        if (s & si) {
            loc_t p4 = sym_transloc(si, span, p1);
            loc_t p5 = sym_transloc(si, span, p2);
            loc_t p6 = sym_transloc(si, span, p3);
            int i4 = list2p_find_lim(point, points, (loc2p_t){ p4, power });
            int i5 = list2p_find_lim(point, points, (loc2p_t){ p5, power });
            int i6 = list2p_find_lim(point, points, (loc2p_t){ p6, power });
            if (i4 > i5 && i4 > i6) {
                if (i4 < i1) return 0;
                if (i4 == i1) {
                    if (i5 > i6) {
                        if (i5 < i2 || (i5 == i2 && i6 < i3))
                            return 0;
                    } else {
                        if (i6 < i2 || (i6 == i2 && i5 < i3))
                            return 0;
                    }
                }
            } else if (i5 > i4 && i5 > i6) {
                if (i5 < i1) return 0;
                if (i5 == i1) {
                    if (i4 > i6) {
                        if (i4 < i2 || (i4 == i2 && i6 < i3))
                            return 0;
                    } else {
                        if (i6 < i2 || (i6 == i2 && i4 < i3))
                            return 0;
                    }
                }
            } else {
                if (i6 < i1) return 0;
                if (i6 == i1) {
                    if (i4 > i5) {
                        if (i4 < i2 || (i4 == i2 && i5 < i3))
                            return 0;
                    } else {
                        if (i5 < i2 || (i5 == i2 && i4 < i3))
                            return 0;
                    }
                }
            }
        }
    }
    return 1;
}

/* Try all possible extensions of the existing 'points'-point arrangement. */
void try_next(int points, int new) {
    cx_t *cx = &context[points];
    cx_t *cx1 = &context[points + 1];
    cx_t *cx2 = &context[points + 2];
    loclist_t *seen = dup_loclist(cx->seen);
    pairlist_t *pairs = dup_pairlist(cx->pairs);
    int power = cx->power;
    int raising = (power > context[points - new].power) ? 1 : 0;
    sym_t sym = sym_check(point, cx->span, points, power);

    if (lim_visit && visit >= lim_visit)
        return;
    if (verbose == 2 && (quiet == 0 || points == n))
        report(points);
    ++visit;

    if (raising) {
        for (int i = 0; i < seen->used; ++i)
            list_set(seen, i, raise2p(list_get(seen, i)));
        for (int i = 0; i < pairs->used; ++i) {
            pair_t p = pair_get(pairs, i);
            pair_set(pairs, i, (pair_t){ raise2p(p.p[0]), raise2p(p.p[1]) });
        }
    }

    for (int i = points - new; i < points; ++i)
        seen_point(seen, pairs, list2p_get(point, i, power), power, 1);

    if (points + 1 <= n) {
        /* Try to extend 3 points into a square. */
        memcpy(cx1, cx, sizeof(cx_t));
        cx1->seen = seen;
        cx1->pairs = pairs;

        /* Try all triples we haven't already tried */
        while (cx1->try3[0] < points) {
            if (
                /* If this triple forms a square with a missing 4th point */
                try_test3(points, cx1->try3, power)
                /* .. and that point isn't on the list to be suppressed */
                && !list_exists(seen, list2p_get(point, points, power))
                /* .. and it's canonical under symmetries of this arrangement */
                && sym_best3(points, sym, cx->span, power,
                        cx1->try3[0], cx1->try3[1], cx1->try3[2])
            ) {
                /* then apply this extension (and recurse) */
                try_with(points, 1);
                /* suppress it from further extensions of this arrangement */
                seen_point(seen, pairs,
                        list2p_get(point, points, power), power, 0);
                /* restore this, it will have been overwritten */
                cx1->squares = cx->squares;
            }
            /* advance to the next untried triple */
            ++cx1->try3[2];
            if (cx1->try3[2] == cx1->try3[1]) {
                cx1->try3[2] = 0;
                ++cx1->try3[1];
                if (cx1->try3[1] == cx1->try3[0]) {
                    cx1->try3[1] = 1;
                    ++cx1->try3[0];
                }
            }
        }
    }

    if (points + 2 <= n) {
        /* Try to extend 2 points into a square, two ways. Any point
         * already tried on its own in the previous section (at this
         * level or higher up the call chain) need not be tried again.
         */
        memcpy(cx2, cx, sizeof(cx_t));
        /* propagate progress made */
        memcpy(&cx2->try3, &cx1->try3, sizeof(cx2->try3));
        cx2->seen = seen;
        cx2->pairs = pairs;

        /* Try all pairs and directions we haven't already tried */
        while (cx2->try2[0] < points) {
            if (
                /* If this pair is canonical for symmetry of this arrangement */
                sym_best2(points, sym, cx->span, power,
                        cx2->try2[0], cx2->try2[1])
                /* .. and its direction is canonical */
                && !(cx2->try2[2] == 1 && sym_axis(sym, cx->span,
                        list2p_get(point, cx2->try2[0], power),
                        list2p_get(point, cx2->try2[1], power)))
                /* .. and it forms a square with two missing points */
                && try_test2(points, cx2->try2, power)
                /* .. and the points are not on a list to be suppressed */
                && !list_exists(seen, list2p_get(point, points, power))
                && !list_exists(seen, list2p_get(point, points + 1, power))
                && !pair_exists(pairs, (pair_t){
                    list2p_get(point, points, power),
                    list2p_get(point, points + 1, power)
                })
            ) {
                /* then apply this extension (and recurse) */
                try_with(points, 2);
                /* suppress the pair from further extensions */
                pair_append(pairs, (pair_t){
                    list2p_get(point, points, power),
                    list2p_get(point, points + 1, power)
                });
                /* restore this, it will have been overwritten */
                cx2->squares = cx->squares;
            }
            /* advance to the next untried pair/dir */
            ++cx2->try2[2];
            if (cx2->try2[2] == 3) {
                cx2->try2[2] = 0;
                ++cx2->try2[1];
                if (cx2->try2[1] == cx2->try2[0]) {
                    cx2->try2[1] = 0;
                    ++cx2->try2[0];
                }
            }
        }
    }
    free_loclist(seen);
    free_pairlist(pairs);
    return;
}

int main(int argc, char** argv) {
    int arg = 1;

    while (arg < argc && argv[arg][0] == '-') {
        char *s = argv[arg++];
        if (strcmp("--", s) == 0)
            break;
        else if (strcmp("-m", s) == 0)
            verbose = 1;
        else if (strcmp("-v", s) == 0)
            verbose = 2;
        else if (strcmp("-q", s) == 0)
            quiet = 1;
        else if (strcmp("-i", s) == 0)
            lim_visit = atol(argv[arg++]);
        else {
            fprintf(stderr, "Unknown option '%s'\n", s);
            exit(1);
        }
    }
    if (arg + 1 != argc) {
        fprintf(stderr, "Usage: try [-v] [-i maxiter] <n>\n");
        return 1;
    }

    n = atoi(argv[arg]);
    if (n < 0) {
        fprintf(stderr, "Need n >= 0\n");
        return 1;
    }
    if (n < 4) {
        printf("%d 0 (0) %.2fs\n", n, timing());
        return 0;
    }
    if (n == 4) {
        printf("4 1 2x2 2x2 (0) %.2fs\n", timing());
        return 0;
    }

    init();
    try_next(4, 0);

    /* report final results */
    printf("%d %d %lu %dx%d %dx%d (%lu) %.2fs\n",
        n, best, last_new,
        (minspan.p.x >> 1) + 1, (minspan.p.y >> 1) + 1,
        (maxspan.p.x >> 1) + 1, (maxspan.p.y >> 1) + 1,
        visit, timing()
    );
    finish();
    return 0;
}
