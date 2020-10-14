#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "common.h"

int t = 0;
int failed = 0;

void init_test(void) {
    ;
}

void done_testing(void) {
    printf("1..%d\n", t);
    if (failed)
        printf("# Looks like you failed %d tests of %d.\n", failed, t);
}

void _ok(char *legend, va_list argp) {
    printf("ok %d - ", ++t);
    vprintf(legend, argp);
    printf("\n");
}

void ok(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _ok(legend, argp);
    va_end(argp);
}

void fail(char *legend, va_list argp) {
    printf("not ok %d - ", ++t);
    vprintf(legend, argp);
    printf("\n");
    ++failed;
}

void print_grid(int x, int y, int *vals) {
    for (int i = 0; i < x; ++i) {
        if (i) printf("; ");
        for (int j = 0; j < y; ++j) {
            if (j) printf(" ");
            printf("%d", vals[i * y + j]);
        }
    }
}

char *_expect(char *str, char *expect) {
    size_t len = strlen(expect);
    if (strncmp(str, expect, len)) {
        fprintf(stderr, "Expected '%s', not '%s'\n", expect, str);
        exit(1);
    }
    return str + len;
}

int *parse_vals(int x, int y, char *str) {
    int *vals = malloc(x * y * sizeof(int));
    char *next;

    for (int i = 0; i < x; ++i) {
        if (i)
            str = _expect(str, "; ");
        for (int j = 0; j < y; ++j) {
            if (j)
                str = _expect(str, " ");
            vals[i * y + j] = strtod(str, &next);
            if (str == next) {
                fprintf(stderr, "Expected number at '%s'", str);
                exit(1);
            }
            str = next;
        }
    }
    if (*str) {
        fprintf(stderr, "Expected end of string at '%s'", str);
        exit(1);
    }
    return vals;
}

void is_bool(bool got, bool expect, char *legend, ...) {
    va_list argp;

    va_start(argp, legend);
    if (got == expect) {
        _ok(legend, argp);
    } else {
        fail(legend, argp);
        printf("#          got: %s\n#   expected: %s\n",
                got ? "true" : "false", expect ? "true" : "false");
    }
    va_end(argp);
}

void is_int(int got, int expect, char *legend, ...) {
    va_list argp;

    va_start(argp, legend);
    if (got == expect) {
        _ok(legend, argp);
    } else {
        fail(legend, argp);
        printf("#          got: %d\n#   expected: %d\n", got, expect);
    }
    va_end(argp);
}

void is_loc(loc_t got, loc_t expect, char *legend, ...) {
    va_list argp;

    va_start(argp, legend);
    if (got.x == expect.x && got.y == expect.y) {
        _ok(legend, argp);
    } else {
        fail(legend, argp);
        printf("#          got: { %d, %d }\n#   expected: { %d, %d }\n",
                got.x, got.y, expect.x, expect.y);
    }
    va_end(argp);
}

void is_grid(int x, int y, int *got, int *expect, char *legend, ...) {
    int fail_count = 0;
    va_list argp;

    va_start(argp, legend);
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            if (got[i * y + j] != expect[i * y + j])
                ++fail_count;
    if (fail_count == 0) {
        _ok(legend, argp);
    } else {
        fail(legend, argp);
        printf("#          got: "); print_grid(x, y, got); printf("\n");
        printf("#     expected: "); print_grid(x, y, expect); printf("\n");
    }
    va_end(argp);
}

