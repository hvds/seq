#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/times.h>
#include <unistd.h>

#include "loc.h"
#include "sym.h"

#define MAXN 64

long clock_tick;

/* A structure of context used at each level of recursion, which
 * encapsulates information about the arrangement up to this point,
 * and information about what extensions have already been tried.
 * It does not list the points themselves, those are stored globally
 * in point[].
 */
typedef struct {
    int squares;    /* Number of squares in this arrangement */
    loc_t span[2];  /* { xmin, ymin }, { xmax, ymax } of this arrangement */
    int try3[3];    /* The next triple to try (list of 3 point indices) */
    int try2[3];    /* The next pair to try, plus direction */
    loclist_t *seen;/* List of points tried, that we don't need to try again */
} cx_t;

int n;              /* We're trying to find A051602(n) */
int best;           /* Greatest number of squares seen in any arrangement */
loclist_t *point;   /* List of points in the current arrangement */
cx_t context[MAXN]; /* List of context objects */
loc_t minspan;      /* { x, y } size of smallest maximal solution */
loc_t maxspan;      /* { x, y } size of greatest maximal solution */
unsigned long visit;/* Count of iterations */

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

/* Show a result-so-far, consisting of i points */
void report(int i) {
    cx_t *cx = &context[i];
    loc_t span = loc_diff(cx->span[0], cx->span[1]);

    printf("%dx%d ", span.x + 1, span.y + 1);
    printf("%d:", cx->squares);
    for (int j = 0; j < i; ++j) {
        loc_t p = list_get(point, j);
        printf(" %d:%d", p.x, p.y);
    }
    printf("\n");
}

void init(void) {
    loclist_t *first_seen = new_loclist(10);

    /* Start off with 4 points making a unit square */
    point = new_loclist(MAXN);
    list_set(point, 0, (loc_t){ 0, 0 });
    list_set(point, 1, (loc_t){ 0, 1 });
    list_set(point, 2, (loc_t){ 1, 0 });
    list_set(point, 3, (loc_t){ 1, 1 });

    context[4] = (cx_t){
        1,                      /* squares */
        { { 0, 0 }, { 1, 1 } }, /* span */
        /* We know that no three of our first four points will form an
         * empty triple, so we'll start looking for triples only from
         * the next point.
         */
        { 4, 1, 0 },    /* Next triple to look for is point[4] + [1] + [0] */
        { 1, 0, 0 },    /* Next pair to look for is point[1] + [0], in the
                         * first direction (of 2) */
        first_seen      /* Seed with an empty list */
    };

    best = 1;
    minspan = (loc_t){ 1, 1 };
    maxspan = (loc_t){ 1, 1 };
    visit = 0UL;
}

void finish(void) {
    free_loclist(point);
    free_loclist(context[4].seen);
}

/* We have a double recursion: try_next() finds the next extension to try,
 * and calls try_with() which applies that extension and calculates the
 * effects before calling try_next() again.
 */
void try_next(int depth);

/* Attempt to find a given 2-point extension; on failure, returns FALSE;
 * on success, returns TRUE with the 2 new points appended to point[].
 * For each pair of candidate points, we can extend in 2 directions;
 * considering the candidate as an edge of a square, we can make the
 * second edge by rotating it either 90 or 270 degrees.
 */
int try_test2(int points, int try2[3]) {
    loc_t pi = list_get(point, try2[0]);
    loc_t pj = list_get(point, try2[1]);
    int dir = try2[2];
    loc_t pk, pl, diff;

    diff = loc_diff(pi, pj);
    if (dir) {
        pk = loc_rot90(pi, diff);
        pl = loc_rot90(pj, diff);
    } else {
        pk = loc_rot270(pi, diff);
        pl = loc_rot270(pj, diff);
    }
    /* This is a valid extension only if the two new points are not
     * already present in the arrangement.
     */
    if (list_exists_lim(point, points, pk)
        || list_exists_lim(point, points, pl)
    )
        return 0;
    list_set(point, points, pk);
    list_set(point, points + 1, pl);
    return 1;
}

/* Attempt to find a given 1-point extension; on failure, returns FALSE;
 * on success, returns TRUE with the new points appended to point[].
 * We need to determine a way of joining the 3 input points as two edges
 * so that we can extrapolate a fourth point making a square; we do this
 * by attempting to rotate one edge by 90 or 270 degrees to make it
 * coincide with the third point.
 */
int try_test3(int points, int try3[3]) {
    loc_t pi = list_get(point, try3[0]);
    loc_t pj = list_get(point, try3[1]);
    loc_t pk = list_get(point, try3[2]);
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
    if (list_exists_lim(point, points, pl))
        return 0;
    list_set(point, points, pl);
    return 1;
}

/* Count the new squares formed by adding point[i] to the previous
 * arrangement. To avoid double-counting, we look at the two edges
 * point[i]-point[j] and point[i]-point[k], and count only when
 * j < k.
 */
int find_squares(int i) {
    loc_t pi = list_get(point, i);
    int count = 0;

    for (int j = 0; j < i; ++j) {
        loc_t pj = list_get(point, j);
        loc_t diff = loc_diff(pi, pj);
        loc_t pk = loc_rot90(pi, diff);
        loc_t pl = loc_rot90(pj, diff);
        if (list_find_lim(point, i, pk) > j
           && list_exists_lim(point, i, pl)
        )
            ++count;

        pk = loc_rot270(pi, diff);
        pl = loc_rot270(pj, diff);
        if (list_find_lim(point, i, pk) > j
           && list_exists_lim(point, i, pl)
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

    /* The caller will already have filled in most of the new context,
     * but some parts may have been overwritten by subsequent work.
     */
    ncx->span[0] = ocx->span[0];
    ncx->span[1] = ocx->span[1];

    for (int i = points; i < points + new; ++i) {
        loc_t pi = list_get(point, i);

        /* Update the count of squares */
        ncx->squares += find_squares(i);

        /* Track the 4 limits of the arrangement */
        if (ncx->span[0].x > pi.x) ncx->span[0].x = pi.x;
        if (ncx->span[0].y > pi.y) ncx->span[0].y = pi.y;
        if (ncx->span[1].x < pi.x) ncx->span[1].x = pi.x;
        if (ncx->span[1].y < pi.y) ncx->span[1].y = pi.y;
    }

    if (ncx->squares >= best) {
        loc_t span = loc_diff(ncx->span[0], ncx->span[1]);

        /* minspan/maxspan are canonicalized to call the larger dimension 'x' */
        if (span.x < span.y) {
            int tmp = span.x;
            span.x = span.y;
            span.y = tmp;
        }
        if (ncx->squares > best) {
            /* New record: reset minspan/maxspan, and always report this */
            minspan = span;
            maxspan = span;
            report(points + new);
        } else {
            /* Match of existing record, report it only if new gridsize */
            int newspan = 0;
            if (minspan.x > span.x) minspan.x = span.x, newspan = 1;
            if (minspan.y > span.y) minspan.y = span.y, newspan = 1;
            if (maxspan.x < span.x) maxspan.x = span.x, newspan = 1;
            if (maxspan.y < span.y) maxspan.y = span.y, newspan = 1;
            if (newspan)
                report(points + new);
        }
        best = ncx->squares;
    }
    try_next(points + new);
}

/* Try all possible extensions of the existing 'points'-point arrangement. */
void try_next(int points) {
    cx_t *cx = &context[points];
    cx_t *cx1 = &context[points + 1];
    cx_t *cx2 = &context[points + 2];
    loclist_t *seen = dup_loclist(cx->seen);
    /* Used for symmetry checks: centre point is span/2 */
    loc_t span = loc_sum(cx->span[0], cx->span[1]);
    sym_t sym = sym_check(point, span, points);

    ++visit;

    if (points + 1 <= n) {
        /* Try to extend 3 points into a square. */
        memcpy(cx1, cx, sizeof(cx_t));
        cx1->seen = seen;

        /* Try all triples we haven't already tried */
        while (cx1->try3[0] < points) {
            if (
                /* If this triple forms a square with a missing 4th point */
                try_test3(points, cx1->try3)
                /* .. and that point isn't on the list to be suppressed */
                && !list_exists(seen, list_get(point, points))
                /* .. and it's canonical under symmetries of this arrangement */
                && sym_best(sym, span, list_get(point, points))
            ) {
                /* then apply this extension (and recurse) */
                try_with(points, 1);
                /* suppress it from further extensions of this arrangement */
                list_append(seen, list_get(point, points));
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

        /* Try all pairs and directions we haven't already tried */
        while (cx2->try2[0] < points) {
            if (
                /* If this pair/dir forms a square with two missing points */
                try_test2(points, cx2->try2)
                /* .. and neither new point is on the list to be suppressed */
                && !list_exists(seen, list_get(point, points))
                && !list_exists(seen, list_get(point, points + 1))
                /* .. and it's canonical under symmetries of this arrangement */
                && sym_best2(sym, span, list_get(point, points),
                        list_get(point, points + 1))
            ) {
                /* then apply this extension (and recurse) */
                try_with(points, 2);
                /* restore this, it will have been overwritten */
                cx2->squares = cx->squares;
            }
            /* advance to the next untried pair/dir */
            if (cx2->try2[2] == 0) {
                cx2->try2[2] = 1;
            } else {
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
    return;
}

int main(int argc, char** argv) {
    /* used to calculate timings */
    clock_tick = sysconf(_SC_CLK_TCK);
    /* don't buffer output, even if stdout is not a tty */
    setvbuf(stdout, (char*)NULL, _IONBF, 0);

    if (argc != 2) {
        fprintf(stderr, "Usage: try <n>\n");
        return 1;
    }
    n = atoi(argv[1]);
    if (n < 0 || n > MAXN) {
        fprintf(stderr, "Need 0 <= n <= %d\n", MAXN);
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
    try_next(4);

    /* report final results */
    printf("%d %d %dx%d %dx%d (%lu) %.2fs\n",
        n, best, minspan.x + 1, minspan.y + 1,
        maxspan.x + 1, maxspan.y + 1, visit, timing()
    );
    finish();
    return 0;
}
