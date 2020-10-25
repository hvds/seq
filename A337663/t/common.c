#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "common.h"

int t = 0;
int failed = 0;

void init_test(void) {
    return;
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

void _fail(char *legend, va_list argp) {
    printf("not ok %d - ", ++t);
    vprintf(legend, argp);
    printf("\n");
    ++failed;
}
void fail(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _fail(legend, argp);
    va_end(argp);
}

void fatal(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vfprintf(stderr, legend, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(1);
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
    if (strncmp(str, expect, len))
        fatal("Expected '%s', not '%s'\n", expect, str);
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
            if (str == next)
                fatal("Expected number at '%s'", str);
            str = next;
        }
    }
    if (*str)
        fatal("Expected end of string at '%s'", str);
    return vals;
}

void _is_bool(bool got, bool expect, char *legend, va_list argp) {
    if (got == expect) {
        _ok(legend, argp);
    } else {
        _fail(legend, argp);
        printf("#          got: %s\n#   expected: %s\n",
                got ? "true" : "false", expect ? "true" : "false");
    }
}
void is_bool(bool got, bool expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _is_bool(got, expect, legend, argp);
    va_end(argp);
}

void _is_int(int got, int expect, char *legend, va_list argp) {
    if (got == expect) {
        _ok(legend, argp);
    } else {
        _fail(legend, argp);
        printf("#          got: %d\n#   expected: %d\n", got, expect);
    }
}
void is_int(int got, int expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _is_int(got, expect, legend, argp);
    va_end(argp);
}

void _is_loc(loc_t got, loc_t expect, char *legend, va_list argp) {
    if (got.x == expect.x && got.y == expect.y) {
        _ok(legend, argp);
    } else {
        _fail(legend, argp);
        printf("#          got: { %d, %d }\n#   expected: { %d, %d }\n",
                got.x, got.y, expect.x, expect.y);
    }
}
void is_loc(loc_t got, loc_t expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _is_loc(got, expect, legend, argp);
    va_end(argp);
}

void _is_grid(int x, int y, int *got, int *expect, char *legend, va_list argp) {
    int fail_count = 0;

    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            if (got[i * y + j] != expect[i * y + j])
                ++fail_count;
    if (fail_count == 0) {
        _ok(legend, argp);
    } else {
        _fail(legend, argp);
        printf("#          got: "); print_grid(x, y, got); printf("\n");
        printf("#     expected: "); print_grid(x, y, expect); printf("\n");
    }
}
void is_grid(int x, int y, int *got, int *expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_grid(x, y, got, expect, legend, argp);
    va_end(argp);
}

void _is_group(group_t *got, group_t *expect, char *legend, va_list argp) {
    int xyfail = failed;
    _is_int(got->x, expect->x, legend, argp);
    _is_int(got->y, expect->y, legend, argp);
    xyfail -= failed;
    _is_int(got->sym, expect->sym, legend, argp);
    if (xyfail)
        fail("skip grid, x/y don't match");
    else
        _is_grid(got->x, got->y, got->vals, expect->vals, legend, argp);
}
void is_group(group_t *got, group_t *expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    _is_group(got, expect, legend, argp);
    va_end(argp);
}
