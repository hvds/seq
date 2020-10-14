#include <stdlib.h>

#include "../sym.h"
#include "../group.h"
#include "common.h"
#include "test_sym.h"

void test_basics(void) {
    sym_t not_tr[4] = { xy, xY, Xy, XY };
    sym_t is_tr[4] = { yx, yX, Yx, YX };

    is_int(MAXSYM + 1, 8, "there are 8 syms");
    ok("all_bits(0xff) finds them all (skipped)");
    for (int i = 0; i < 4; ++i)
        is_bool(is_transpose(not_tr[i]), false,
                "%d is_transpose", not_tr[i]);
    for (int i = 0; i < 4; ++i)
        is_bool(is_transpose(is_tr[i]), true,
                "%d is_transpose", is_tr[i]);
}

void test_sym(void) {
    for (int i = 0; i < data_sym_count; ++i) {
        test_sym_t data = data_sym[i];
        sym_t s = data.s;
        int x = data.x, y = data.y;
        int *vals = parse_vals(x, y, data.vals);

        is_bool(sym_check(s, x, y, vals), true,
                "%d sym_check %s", s, data.vals);

        int *got = sym_transform(s, x, y, vals);
        is_grid(x, y, got, vals,
                "%d sym_transform %s", s, data.vals);

        free(vals);
        free(got);
    }
}

void test_asym(void) {
    for (int i = 0; i < data_asym_count; ++i) {
        test_asym_t data = data_asym[i];
        sym_t s = data.s;
        int x = data.x, y = data.y;
        int *from = parse_vals(x, y, data.from);
        int *to = is_transpose(s)
            ? parse_vals(y, x, data.to)
            : parse_vals(x, y, data.to);

        is_bool(sym_check(s, x, y, from), false,
                "%d sym_check %s", s, data.from);

        int *got = sym_transform(s, x, y, from);
        is_grid(x, y, got, to,
                "%d sym_transform %s", s, data.from);

        free(from);
        free(to);
        free(got);
    }
}

void test_loc(void) {
    for (int i = 0; i < data_loc_count; ++i) {
        test_loc_t data = data_loc[i];
        sym_t s = data.s;
        int x = data.x, y = data.y;
        int xfrom = data.xfrom, yfrom = data.yfrom;
        int xto = data.xto, yto = data.yto;
        loc_t from = (loc_t){ xfrom, yfrom };
        bool same = (xfrom == xto && yfrom == yto);

        loc_t got = sym_transloc(s, x, y, from);
        is_loc(got, (loc_t){ xto, yto },
                "%d sym_transloc %d.%d", s, xfrom, yfrom);

        is_bool(sym_checkloc(s, x, y, from), same,
                "%d sym_checkloc %d.%d", s, xfrom, yfrom);
    }
}

int main(void) {
    init_test();
    test_basics();
    test_sym();
    test_asym();
    test_loc();
    done_testing();
}
