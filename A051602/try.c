#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <sys/times.h>
#include <unistd.h>

#include "loc.h"

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
loclist_t *point;
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
    for (int j = 0; j < i; ++j) {
        loc_t p = list_get(point, j);
        printf(" %d:%d", p.x, p.y);
    }
    printf("\n");
}

void init(void) {
    clock_tick = sysconf(_SC_CLK_TCK);
    setvbuf(stdout, (char*)NULL, _IONBF, 0);

    point = new_loclist(MAXN);
    list_set(point, 0, (loc_t){ 0, 0 });
    list_set(point, 1, (loc_t){ 0, 1 });
    list_set(point, 2, (loc_t){ 1, 0 });
    list_set(point, 3, (loc_t){ 1, 1 });

    context[4] = (cx_t){
        1, 1,
        { { 0, 0 }, { 1, 1 } },
        { 4, 1, 0 },
        /* FIXME: try2[] needs to start [1, 0, 0] but use symmetries */
        { 3, 2, 0 }
    };

    minspan = (loc_t){ 1, 1 };
    maxspan = (loc_t){ 1, 1 };
    visit = 0UL;
}

void finish(void) {
    free_loclist(point);
}

void try_next(int depth);

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
    if (list_exists_lim(point, points, pk)
        || list_exists_lim(point, points, pl)
    )
        return 0;
    list_set(point, points, pk);
    list_set(point, points + 1, pl);
    return 1;
}

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
    if (list_exists_lim(point, points, pl))
        return 0;
    list_set(point, points, pl);
    return 1;
}

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

void try_with(int points, int new) {
    cx_t *ocx = &context[points];
    cx_t *ncx = &context[points + new];

    ncx->span[0] = ocx->span[0];
    ncx->span[1] = ocx->span[1];
    for (int i = points; i < points + new; ++i) {
        loc_t pi = list_get(point, i);
        ncx->squares += find_squares(i);
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
    loclist_t *seen = new_loclist(10); /* will grow as needed */

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
                /* FIXME: do we need to deduplicate here? */
                list_append(seen, list_get(point, points));
                list_append(seen, list_get(point, points + 1));
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
                && !list_exists(seen, list_get(point, points))
            ) {
                try_with(points, 1);
                list_append(seen, list_get(point, points));
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
    free_loclist(seen);
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
    finish();
    return 0;
}
