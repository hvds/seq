#include <stdlib.h>

#include "../sym.h"
#include "../group.h"
#include "common.h"

typedef struct pgroup_s {
    int x;
    int y;
    int syms;
    char* str;
} pgroup_t;

int all_sym = (1 << xY) | (1 << Xy) | (1 << XY) | (1 << yx) | (1 << yX) | (1 << Yx) | (1 << YX);

group_t *_parse_group(pgroup_t pg) {
    int *vals = parse_vals(pg.x, pg.y, pg.str);
    group_t *g = new_group(pg.x, pg.y, pg.syms, vals);
    ref_group(g);
    return g;
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
