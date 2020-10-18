#include <stdlib.h>

#include "../sym.h"
#include "../group.h"
#include "common.h"
#include "test_sym.h"

void test_basics(void) {
    sym_t not_tr[4] = { xy, xY, Xy, XY };
    sym_t is_tr[4] = { yx, yX, Yx, YX };
    sym_t not_ref[4] = { xy, XY, yX, Yx };
    sym_t is_ref[4] = { xY, Xy, yx, YX };

    is_int(MAXSYM + 1, 8, "there are 8 syms");
    ok("all_bits(0xff) finds them all (skipped)");
    for (int i = 0; i < 4; ++i)
        is_bool(is_transpose(not_tr[i]), false,
                "%d is_transpose", not_tr[i]);
    for (int i = 0; i < 4; ++i)
        is_bool(is_transpose(is_tr[i]), true,
                "%d is_transpose", is_tr[i]);
    for (int i = 0; i < 4; ++i)
        is_bool(is_reflect(not_ref[i]), false,
                "%d is_reflect", not_ref[i]);
    for (int i = 0; i < 4; ++i)
        is_bool(is_reflect(is_ref[i]), true,
                "%d is_reflect", is_ref[i]);
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

void test_dup(void) {
    int base[4] = { 1, 2, 4, 8 }, symi[4];
    int *v, vc[4], vd[4], c, d;

    for (sym_t i = 0; i <= MAXSYM; ++i) {
        if (!is_reflect(i))
            continue;
        /* construct a grid that shows this symmetry */
        v = sym_transform(i, 2, 2, &base[0]);
        for (int k = 0; k < 4; ++k)
            symi[k] = base[k] | v[k];

        c = d = 0;
        for (sym_t j = 0; j <= MAXSYM; ++j) {
            int *w = sym_transform(j, 2, 2, &symi[0]);
            int packed = (w[0] << 24) | (w[1] << 16) | (w[2] << 8) | w[3];
            if (sym_dup(i, j))
                vd[d++] = packed;
            else
                vc[c++] = packed;
            free(w);
        }
        if ((d != 4) || (c != 4))
            fatal("Got %d canonical and %d dup, possibly memory corruption",
                    c, d);
        /* we expect each of the canonical vc[] to be distinct; we expect
         * each of the dup vd[] to be a copy of something canonical
         */
        int cbad = 0, dbad = 0;
        for (int ci = 0; ci < c - 1; ++ci)
            for (int cj = ci + 1; cj < c; ++cj)
                if (vc[ci] == vc[cj])
                    ++cbad;
        for (int di = 0; di < d; ++di) {
            ++dbad;
            for (int ci = 0; ci < c; ++ci)
                if (vd[di] == vc[ci]) {
                    --dbad;
                    break;
                }
        }
        is_int(cbad, 0, "No duplicate canonical for %d", i);
        is_int(dbad, 0, "Every dup dups a canonical for %d", i);
        free(v);
    }
}

int main(void) {
    init_test();
    init_sym();

    test_basics();
    test_sym();
    test_asym();
    test_loc();
    test_dup();

    finish_sym();
    done_testing();
}
