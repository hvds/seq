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
        printf("# Looks like you failed %d test%s of %d.\n",
                failed, failed == 1 ? "" : "s", t);
}

void vok(char *legend, va_list argp) {
    printf("ok %d - ", ++t);
    vprintf(legend, argp);
    printf("\n");
}
void ok(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vok(legend, argp);
    va_end(argp);
}

void vfail(char *legend, va_list argp) {
    printf("not ok %d - ", ++t);
    vprintf(legend, argp);
    printf("\n");
    ++failed;
}
void fail(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vfail(legend, argp);
    va_end(argp);
}

void vfatal(char *legend, va_list argp) {
    vfprintf(stderr, legend, argp);
    fprintf(stderr, "\n");
    exit(1);
}
void fatal(char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vfatal(legend, argp);
    va_end(argp);
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

void vis_bool(bool got, bool expect, char *legend, va_list argp) {
    if (got == expect) {
        vok(legend, argp);
    } else {
        vfail(legend, argp);
        printf("#          got: %s\n#   expected: %s\n",
                got ? "true" : "false", expect ? "true" : "false");
    }
}
void is_bool(bool got, bool expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_bool(got, expect, legend, argp);
    va_end(argp);
}

void vis_int(int got, int expect, char *legend, va_list argp) {
    if (got == expect) {
        vok(legend, argp);
    } else {
        vfail(legend, argp);
        printf("#          got: %d\n#   expected: %d\n", got, expect);
    }
}
void is_int(int got, int expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_int(got, expect, legend, argp);
    va_end(argp);
}

void vis_loc(loc_t got, loc_t expect, char *legend, va_list argp) {
    if (got.x == expect.x && got.y == expect.y) {
        vok(legend, argp);
    } else {
        vfail(legend, argp);
        printf("#          got: { %d, %d }\n#   expected: { %d, %d }\n",
                got.x, got.y, expect.x, expect.y);
    }
}
void is_loc(loc_t got, loc_t expect, char *legend, ...) {
    va_list argp;
    va_start(argp, legend);
    vis_loc(got, expect, legend, argp);
    va_end(argp);
}

void vis_grid(int x, int y, int *got, int *expect, char *legend, va_list argp) {
    int fail_count = 0;

    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            if (got[i * y + j] != expect[i * y + j])
                ++fail_count;
    if (fail_count == 0) {
        vok(legend, argp);
    } else {
        vfail(legend, argp);
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
