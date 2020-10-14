#include "../sym.h"

typedef struct test_sym_s {
    sym_t s;
    int x;
    int y;
    char *vals;
} test_sym_t;

#define data_sym_count 20
test_sym_t data_sym[data_sym_count] = {
    (test_sym_t){ xy, 2, 2, "1 2; 3 4" },
    (test_sym_t){ xy, 1, 4, "1 3 5 7" },
    (test_sym_t){ xy, 4, 1, "1; 3; 5; 7" },
    (test_sym_t){ xY, 4, 1, "1; 2; 3; 4" },
    (test_sym_t){ xY, 3, 2, "1 1; 2 2; 3 3" },
    (test_sym_t){ xY, 2, 3, "1 3 1; 2 4 2" },
    (test_sym_t){ Xy, 1, 4, "1 2 3 4" },
    (test_sym_t){ Xy, 2, 3, "1 2 3; 1 2 3" },
    (test_sym_t){ Xy, 3, 2, "1 2; 3 4; 1 2" },
    (test_sym_t){ XY, 2, 3, "1 1 0; 0 1 1" },
    (test_sym_t){ XY, 3, 3, "1 1 0; 0 0 0; 0 1 1" },
    (test_sym_t){ XY, 3, 2, "1 0; 1 1; 0 1" },
    (test_sym_t){ yx, 2, 2, "1 2; 2 3" },
    (test_sym_t){ yx, 3, 3, "1 3 5; 3 9 13; 5 13 21" },
    (test_sym_t){ yX, 2, 2, "1 1; 1 1" },
    (test_sym_t){ yX, 3, 3, "1 0 1; 0 2 0; 1 0 1" },
    (test_sym_t){ Yx, 2, 2, "1 1; 1 1" },
    (test_sym_t){ Yx, 3, 3, "1 0 1; 0 2 0; 1 0 1" },
    (test_sym_t){ YX, 2, 2, "0 1; 1 0" },
    (test_sym_t){ YX, 3, 3, "0 1 1; 3 0 1; 5 3 0" }
};

typedef struct test_asym_s {
    sym_t s;
    int x;
    int y;
    char *from;
    char *to;
} test_asym_t;

#define data_asym_count 19
test_asym_t data_asym[data_asym_count] = {
    (test_asym_t){ xY, 2, 2, "1 2; 3 4", "2 1; 4 3" },
    (test_asym_t){ xY, 1, 4, "1 2 3 4", "4 3 2 1" },
    (test_asym_t){ Xy, 2, 2, "1 2; 3 4", "3 4; 1 2" },
    (test_asym_t){ Xy, 4, 1, "1; 2; 3; 4", "4; 3; 2; 1" },
    (test_asym_t){ XY, 2, 2, "1 2; 3 4", "4 3; 2 1" },
    (test_asym_t){ XY, 4, 1, "1; 2; 3; 4", "4; 3; 2; 1" },
    (test_asym_t){ XY, 1, 4, "1 2 3 4", "4 3 2 1" },
    (test_asym_t){ yx, 2, 2, "1 2; 3 4", "1 3; 2 4" },
    (test_asym_t){ yx, 4, 1, "1; 2; 3; 4", "1 2 3 4" },
    (test_asym_t){ yx, 1, 4, "1 2 3 4", "1; 2; 3; 4" },
    (test_asym_t){ yX, 2, 2, "1 2; 3 4", "3 1; 4 2" },
    (test_asym_t){ yX, 4, 1, "1; 2; 3; 4", "4 3 2 1" },
    (test_asym_t){ yX, 1, 4, "1 2 3 4", "1; 2; 3; 4" },
    (test_asym_t){ Yx, 2, 2, "1 2; 3 4", "2 4; 1 3" },
    (test_asym_t){ Yx, 4, 1, "1; 2; 3; 4", "1 2 3 4" },
    (test_asym_t){ Yx, 1, 4, "1 2 3 4", "4; 3; 2; 1" },
    (test_asym_t){ YX, 2, 2, "1 2; 3 4", "4 2; 3 1" },
    (test_asym_t){ YX, 4, 1, "1; 2; 3; 4", "4 3 2 1" },
    (test_asym_t){ YX, 1, 4, "1 2 3 4", "4; 3; 2; 1" }
};

typedef struct test_loc_s {
    sym_t s;
    int x;
    int y;
    int xfrom;
    int yfrom;
    int xto;
    int yto;
} test_loc_t;

#define data_loc_count 21
test_loc_t data_loc[data_loc_count] = {
    (test_loc_t){ xy, 3, 3, 0, 0, 0, 0 },
    (test_loc_t){ xy, 9, 1, 2, 2, 2, 2 },
    (test_loc_t){ xy, 9, 1, -2, -2, -2, -2 },
    (test_loc_t){ xY, 3, 3, 0, 0, 0, 2 },
    (test_loc_t){ xY, 9, 1, 2, 2, 2, -2 },
    (test_loc_t){ xY, 9, 1, -2, -2, -2, 2 },
    (test_loc_t){ Xy, 3, 3, 0, 0, 2, 0 },
    (test_loc_t){ Xy, 9, 1, 2, 2, 6, 2 },
    (test_loc_t){ Xy, 9, 1, -2, -2, 10, -2 },
    (test_loc_t){ XY, 3, 3, 0, 0, 2, 2 },
    (test_loc_t){ XY, 9, 1, 2, 2, 6, -2 },
    (test_loc_t){ XY, 9, 1, -2, -2, 10, 2 },
    (test_loc_t){ yX, 3, 3, 0, 0, 0, 2 },
    (test_loc_t){ yX, 9, 1, 2, 2, 2, 6 },
    (test_loc_t){ yX, 9, 1, -2, -2, -2, 10 },
    (test_loc_t){ Yx, 3, 3, 0, 0, 2, 0 },
    (test_loc_t){ Yx, 9, 1, 2, 2, -2, 2 },
    (test_loc_t){ Yx, 9, 1, -2, -2, 2, -2 },
    (test_loc_t){ YX, 3, 3, 0, 0, 2, 2 },
    (test_loc_t){ YX, 9, 1, 2, 2, -2, 6 },
    (test_loc_t){ YX, 9, 1, -2, -2, 2, 10 }
};
