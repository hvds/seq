#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <sys/times.h>
#include <unistd.h>

#include "loc.h"
#include "grid.h"

#define MAXN 64

long clock_tick;

typedef struct {
    int squares;
    int best;
    loc_t span[2];
    int try3[3];
    int try2[3];
} cx_t;

int n;
loc_t point[MAXN];
grid_t *grid[MAXN];
cx_t context[MAXN];
loc_t minspan;
loc_t maxspan;
unsigned long visit;

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

void report(int i) {
    cx_t *cx = &context[i];
    loc_t span = loc_diff(cx->span[0], cx->span[1]);

    printf("[ %d %d ] ", span.x + 1, span.y + 1);
    printf("%d:", cx->squares);
    for (int j = 0; j < i; ++j)
        printf(" %d:%d", point[j].x, point[j].y);
    printf("\n");
}

void init(void) {
    clock_tick = sysconf(_SC_CLK_TCK);

    /* FIXME: try2[] needs to start [1, 0, 0] but use symmetries */
    cx_t start = { 1, 1, { { 0, 0 }, { 1, 1 } }, { 4, 1, 0 }, { 3, 2, 0 } };
    point[0].x = 0; point[0].y = 0;
    point[1].x = 0; point[1].y = 1;
    point[2].x = 1; point[2].y = 0;
    point[3].x = 1; point[3].y = 1;
    grid[4] = new_grid(0, 1, 0, 1);
    for (int i = 0; i < 4; ++i)
        set_grid(grid[4], point[i].x, point[i].y, i + 1);
    memcpy(&context[4], &start, sizeof(cx_t));
    minspan.x = 1; minspan.y = 1;
    maxspan.x = 1; maxspan.y = 1;
    visit = 0UL;
}

void try_next(int depth);

int try_test2(int points, int try2[3]) {
    loc_t pi = point[try2[0]];
    loc_t pj = point[try2[1]];
    int dir = try2[2];
    loc_t *pk = &point[points];
    loc_t *pl = &point[points + 1];
    loc_t diff;

    diff = loc_diff(pi, pj);
    if (dir) {
        *pk = loc_rot90(pi, diff);
        *pl = loc_rot90(pj, diff);
    } else {
        *pk = loc_rot270(pi, diff);
        *pl = loc_rot270(pj, diff);
    }
    if (get_grid(grid[points], pk->x, pk->y) || get_grid(grid[points], pl->x, pl->y))
        return 0;
    return 1;
}

int try_test3(int points, int try3[3]) {
    loc_t pi = point[try3[0]];
    loc_t pj = point[try3[1]];
    loc_t pk = point[try3[2]];
    loc_t *pl = &point[points];
    loc_t diff;

    diff = loc_diff(pi, pj);
    if (loc_eq(pk, loc_rot90(pi, diff))) {
        *pl = loc_rot90(pj, diff);
    } else if (loc_eq(pk, loc_rot90(pj, diff))) {
        *pl = loc_rot90(pi, diff);
    } else if (loc_eq(pk, loc_rot270(pi, diff))) {
        *pl = loc_rot270(pj, diff);
    } else if (loc_eq(pk, loc_rot270(pj, diff))) {
        *pl = loc_rot270(pi, diff);
    } else {
        diff = loc_diff(pi, pk);
        if (loc_eq(pj, loc_rot90(pk, diff))) {
            *pl = loc_rot90(pi, diff);
        } else if (loc_eq(pj, loc_rot270(pk, diff))) {
            *pl = loc_rot270(pi, diff);
        } else {
            return 0;
        }
    }
    if (get_grid(grid[points], pl->x, pl->y))
        return 0;
    return 1;
}

int find_squares(grid_t *g, int i) {
    loc_t pi = point[i];
    int count = 0;

    for (int j = 0; j < i; ++j) {
        loc_t pj = point[j];
        loc_t diff = loc_diff(pi, pj);
        loc_t pk = loc_rot90(pi, diff);
        loc_t pl = loc_rot90(pj, diff);
        if (get_grid(g, pk.x, pk.y) > j
            && get_grid(g, pl.x, pl.y)
        )
            ++count;

        pk = loc_rot270(pi, diff);
        pl = loc_rot270(pj, diff);
        if (get_grid(g, pk.x, pk.y) > j
            && get_grid(g, pl.x, pl.y)
        )
            ++count;
    }
    return count;
}

void try_with(int points, int new) {
    cx_t *ocx = &context[points];
    cx_t *ncx = &context[points + new];
    grid_t *og = grid[points];
    grid_t *ng = grid[points + new];

    if (ng)
        free_grid(ng);
    ng = grid[points + new] = dup_grid(og);

    ncx->span[0] = ocx->span[0];
    ncx->span[1] = ocx->span[1];
    for (int i = points; i < points + new; ++i) {
        loc_t pi = point[i];
        set_grid(ng, pi.x, pi.y, i + 1);
        ncx->squares += find_squares(ng, i);
        if (ncx->span[0].x > pi.x) ncx->span[0].x = pi.x;
        if (ncx->span[0].y > pi.y) ncx->span[0].y = pi.y;
        if (ncx->span[1].x < pi.x) ncx->span[1].x = pi.x;
        if (ncx->span[1].y < pi.y) ncx->span[1].y = pi.y;
    }
    if (ncx->squares >= ncx->best) {
        loc_t span = loc_diff(ncx->span[0], ncx->span[1]);
        if (span.x < span.y) {
            int tmp = span.x;
            span.x = span.y;
            span.y = tmp;
        }
        if (ncx->squares > ncx->best) {
            minspan = span;
            maxspan = span;
            report(points + new);
        } else {
            int newspan = 0;
            if (minspan.x > span.x) minspan.x = span.x, newspan = 1;
            if (minspan.y > span.y) minspan.y = span.y, newspan = 1;
            if (maxspan.x < span.x) maxspan.x = span.x, newspan = 1;
            if (maxspan.y < span.y) maxspan.y = span.y, newspan = 1;
            if (newspan)
                report(points + new);
        }
        ocx->best = ncx->best = ncx->squares;
    }
    try_next(points + new);
}

void try_next(int points) {
    cx_t *cx = &context[points];
    cx_t *cx1 = &context[points + 1];
    cx_t *cx2 = &context[points + 2];
    grid_t *og = grid[points];
    grid_t *seen = new_grid(og->xmin, og->xmax, og->ymin, og->ymax);

    ++visit;

    /* FIXME: we must try3() first, and can suppress try2() attempts for
     * any point seen during that phase */
    if (points + 2 <= n) {
        /* Try to extend 2 points into a square, two ways. Any of the points
         * tried here need not be re-tried in the following section.
         */
        memcpy(cx2, cx, sizeof(cx_t));
        while (cx2->try2[0] < points) {
            if (try_test2(points, cx2->try2)) {
                try_with(points, 2);
                set_grid(seen, point[points].x, point[points].y, 1);
                set_grid(seen, point[points + 1].x, point[points + 1].y, 1);
                cx2->squares = cx->squares;
            }
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

        if (cx->best < cx2->best) {
            cx->best = cx2->best;
        }
    }

    if (points + 1 <= n) {
        /* Try to extend 3 points into a square. */
        memcpy(cx1, cx, sizeof(cx_t));
        if (points + 2 <= n) {
            /* propagate progress made */
            memcpy(&cx1->try2, &cx2->try2, sizeof(cx1->try2));
        }
        while (cx1->try3[0] < points) {
            if (try_test3(points, cx1->try3)
                && !get_grid(seen, point[points].x, point[points].y)
            ) {
                try_with(points, 1);
                set_grid(seen, point[points].x, point[points].y, 1);
                cx1->squares = cx->squares;
            }
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

        if (cx->best < cx1->best) {
            cx->best = cx1->best;
        }
    }
    return;
}

int main(int argc, char** argv) {
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
        printf("%d 0 (0)\n", n);
        return 0;
    }
    if (n == 4) {
        printf("4 1 [ 2 2 ] [ 2 2 ] (0)\n");
        return 0;
    }

    init();
    try_next(4);
    printf("%d %d [ %d %d ] [ %d %d ] (%lu) %.2fs\n",
        n, context[4].best, minspan.x + 1, minspan.y + 1,
        maxspan.x + 1, maxspan.y + 1, visit, timing()
    );
    return 0;
}
