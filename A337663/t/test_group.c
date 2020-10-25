#include <stdlib.h>
#include <stdarg.h>

#include "../sym.h"
#include "../group.h"
#include "common.h"

typedef struct pgroup_s {
    int x;
    int y;
    int syms;
    char* str;
} pgroup_t;

int all_sym = (1 << xY) | (1 << Xy) | (1 << XY)
        | (1 << yx) | (1 << yX) | (1 << Yx) | (1 << YX);

group_t *_parse_group(pgroup_t pg) {
    int *vals = parse_vals(pg.x, pg.y, pg.str);
    group_t *g = new_group(pg.x, pg.y, pg.syms, vals);
    ref_group(g);
    return g;
}

void vis_group(group_t *got, group_t *expect, char *legend, va_list argp) {
    int bad = 0, broken = 0;
    if (got->x != expect->x) {
        broken = 1;
        printf("# got x=%d, expected %d\n", got->x, expect->x);
    }
    if (got->y != expect->y) {
        broken = 1;
        printf("# got y=%d, expected %d\n", got->y, expect->y);
    }
    if (got->sym != expect->sym) {
        bad = 1;
        printf("# got sym=%d, expected %d\n", got->sym, expect->sym);
    }
    if (!broken)
        for (int i = 0; i < got->x; ++i)
            for (int j = 0; j < got->y; ++j)
                if (got->vals[i * got->y + j] != expect->vals[i * got->y + j]) {
                    bad = 1;
                    printf("# got value(%d, %d)=%d, expected %d\n",
                        i, j, got->vals[i * got->y + j],
                        expect->vals[i * got->y + j]
                    );
                }
    if (bad || broken)
        vfail(legend, argp);
    else
        vok(legend, argp);
}
void is_group(group_t *got, group_t *expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_group(got, expect, legend, argp);
    va_end(argp);
}

int _cmpgroup(const void *a, const void *b) {
    group_t *ga = *(group_t **)a, *gb = *(group_t **)b;
    int result;

    result = ga->x - gb->x;
    if (result == 0)
        result = ga->y - gb->y;
    for (int i = 0; i < ga->x; ++i)
        for (int j = 0; j < ga->y; ++j)
            if (result == 0)
                result = ga->vals[i * ga->y + j] - gb->vals[i * ga->y + j];
    return result;
}

void vis_grouplist(grouplist_t *got, grouplist_t *expect, char *legend, va_list argp) {
    qsort(&got->g[0], got->count, sizeof(group_t *), _cmpgroup);
    qsort(&expect->g[0], expect->count, sizeof(group_t *), _cmpgroup);

    int bad = 0, gi = 0, ei = 0;
    while (gi < got->count || ei < expect->count) {
        int diff = (gi < got->count)
            ? (ei < got->count)
                ? _cmpgroup(&got->g[gi], &expect->g[ei])
                : -1
            : 1;
        if (diff < 0) {
            printf("# Got unexpected group ");
            print_group(got->g[gi]);
            printf("\n");
            bad = 1;
            ++gi;
        } else if (diff > 0) {
            printf("# Missed expected group ");
            print_group(expect->g[ei]);
            printf("\n");
            bad = 1;
            ++ei;
        } else {
            if (got->g[gi]->sym != expect->g[ei]->sym) {
                printf("# (TODO) Got sym=%d, expected %d in ",
                        got->g[gi]->sym, expect->g[ei]->sym);
                print_group(expect->g[ei]);
                printf("\n");
                /* TODO */
                /* bad = 1; */
            }
            ++gi, ++ei;
        }
    }
    if (bad)
        vfail(legend, argp);
    else
        vok(legend, argp);
}
void is_grouplist(grouplist_t *got, grouplist_t *expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_grouplist(got, expect, legend, argp);
    va_end(argp);
}

void test_place(void) {
    group_t *g = _parse_group((pgroup_t){ 1, 1, all_sym, "2" });

    group_t *got = group_place(g, (loc_t){ -1, 0 }, 3);
    ref_group(got);
    group_t *expect = _parse_group((pgroup_t){ 2, 1, 1 << xY, "3; 2" });
    is_group(got, expect, "group_place above");
    unref_group(got);
    unref_group(expect);

    got = group_place(g, (loc_t){ 0, 1 }, 3);
    ref_group(got);
    expect = _parse_group((pgroup_t){ 1, 2, 1 << Xy, "2 3" });
    is_group(got, expect, "group_place right");
    unref_group(got);
    unref_group(expect);

    unref_group(g);
}

void test_place_with(void) {
    group_t *g = _parse_group((pgroup_t){ 1, 1, all_sym, "2" });

    grouplist_t *got = group_place_with(g, (loc_t){ -1, -1 }, 3, 2);
    pgroup_t pexpect[6] = {
        (pgroup_t){ 2, 3, 0, "1 3 0; 1 0 2" },
        (pgroup_t){ 3, 3, 0, "1 0 0; 0 3 0; 1 0 2" },
        (pgroup_t){ 3, 3, yx, "0 0 1; 0 3 0; 1 0 2" },
        (pgroup_t){ 3, 3, 0, "1 0 0; 1 3 0; 0 0 2" },
        (pgroup_t){ 3, 3, yx, "0 1 0; 1 3 0; 0 0 2" },
        (pgroup_t){ 3, 3, 0, "0 0 1; 1 3 0; 0 0 2" },
    };
    grouplist_t *expect = new_grouplist(6);
    for (int i = 0; i < 6; ++i)
        expect->g[i] = _parse_group(pexpect[i]);
    is_grouplist(got, expect, "group_place_with above left with 2");
    free_grouplist(got);
    free_grouplist(expect);

    unref_group(g);
}

/*

--g1--  --g2--  place  with 1s  reflect
   *        *
 2      1 3       4       2       YX
        1 .

         xy  (nothing, groups clash)

             1 1 . .   1 . 1 .   . 1 1 .
. * . .  ->  . 4 . .   . 4 . .   . 4 . .  reflects to yX
2 . 3 1  xY  2 . 3 1   2 . 3 1   2 . 3 1
. . . 1      . . . 1   . . . 1   . . . 1

1 . .        1 . . .   1 . . .   1 . . .
1 3 .    ->  1 3 . 1   1 3 . 1   1 3 . .  reflects to Yx
. . *    Xy  . . 4 1   . . 4 .   . . 4 1
. 2 .        . 2 . .   . 2 . 1   . 2 . 1

. . . 1      . . . 1
. . 3 1  ->  1 . 3 1  reflects to YX
. * . .  XY  . 4 . .
2 . . .      2 . 1 .

         yx  (nothing, groups clash)

1 1 .        1 1 . .   1 1 . .   1 1 . .
. 3 .    ->  . 3 . 1   . 3 . 1   . 3 . .  reflects to xY
. . *    yX  . . 4 1   . . 4 .   . . 4 1
. 2 .        . 2 . .   . 2 . 1   . 2 . 1

             1 1 . .   1 . 1 .   . 1 1 .
. * . .  ->  . 4 . .   . 4 . .   . 4 . .  reflects to Xy
2 . 3 .  Yx  2 . 3 .   2 . 3 .   2 . 3 .
. . 1 1      . . 1 1   . . 1 1   . . 1 1

. . 1 1      . . 1 1
. . 3 .  ->  1 . 3 .  reflects to XY
. * . .  yx  . 4 . .
2 . . .      2 . 1 .

*/
void test_coalesce(void) {
    group_t *g1 = _parse_group((pgroup_t){ 1, 1, all_sym, "2" });
    group_t *g2 = _parse_group((pgroup_t){ 2, 2, 0, "1 3; 1 0" });

    grouplist_t *got = coalesce_group(
        g1, (loc_t){ -1, 1 }, g2, (loc_t){ -1, 2 }, 4, 2
    );
    pgroup_t pexpect[7] = {
        /* xY */
        (pgroup_t){ 4, 4, 0, "1 1 0 0; 0 4 0 0; 2 0 3 1; 0 0 0 1" },
        (pgroup_t){ 4, 4, 0, "1 0 1 0; 0 4 0 0; 2 0 3 1; 0 0 0 1" },
        (pgroup_t){ 4, 4, 0, "0 1 1 0; 0 4 0 0; 2 0 3 1; 0 0 0 1" },
        /* Xy */
        (pgroup_t){ 4, 4, 0, "1 0 0 0; 1 3 0 1; 0 0 4 1; 0 2 0 0" },
        (pgroup_t){ 4, 4, 0, "1 0 0 0; 1 3 0 1; 0 0 4 0; 0 2 0 1" },
        (pgroup_t){ 4, 4, 0, "1 0 0 0; 1 3 0 0; 0 0 4 1; 0 2 0 1" },
        /* YX */
        (pgroup_t){ 4, 4, 0, "0 0 1 1; 1 0 3 0; 0 4 0 0; 2 0 1 0" },
    };
    grouplist_t *expect = new_grouplist(7);
    for (int i = 0; i < 7; ++i)
        expect->g[i] = _parse_group(pexpect[i]);
    is_grouplist(got, expect, "group_coalesce with 2");
    free_grouplist(got);
    free_grouplist(expect);

    unref_group(g1);
    unref_group(g2);
}

int main(void) {
    init_test();
    init_sym();
    init_group();

    test_place();
    test_place_with();
    test_coalesce();

    finish_group();
    finish_sym();
    done_testing();
}
