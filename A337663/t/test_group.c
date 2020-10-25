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

void _is_group(group_t *got, group_t *expect, char *legend, va_list argp) {
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
    _is_group(got, expect, legend, argp);
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

void _is_grouplist(grouplist_t *got, grouplist_t *expect, char *legend, va_list argp) {
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
    _is_grouplist(got, expect, legend, argp);
    va_end(argp);
}

void test_place(void) {
    group_t *g = _parse_group((pgroup_t){ 1, 1, all_sym, "2" });
    group_t *got = group_place(g, (loc_t){ -1, -1 }, 3);
    ref_group(got);
    group_t *expect = _parse_group((pgroup_t){ 2, 2, 1 << yx, "3 0; 0 2" });
    is_group(got, expect, "group_place above left");
    unref_group(g);
    unref_group(got);
    unref_group(expect);
}

int main(void) {
    init_test();
    init_sym();
    init_group();

    test_place();

    finish_group();
    finish_sym();
    done_testing();
}
