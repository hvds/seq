#include <sys/times.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

int clk_tck;
int gtime;

int a, b, ab;
unsigned long long count;

#define MAX_A 8
#define MAX_B 8
#define SIZE_A 1
typedef struct vec_s {
    unsigned char c[SIZE_A * MAX_B];
} vec_t;

vec_t edge[4];
vec_t nb[MAX_A * MAX_B];
vec_t set, unset, surface;

inline void vec_set(vec_t *v, int i, int j) {
    v->c[j] |= (1 << i);
}
inline void vec_unset(vec_t *v, int i, int j) {
    v->c[j] &= ~(1 << i);
}
inline int vec_test(vec_t *v, int i, int j) {
    if (v->c[j] & (1 << i))
        return 1;
    return 0;
}
inline void vec_copy(vec_t *dest, vec_t *src) {
    memcpy((void *)dest, (void *)src, sizeof(vec_t));
}
inline void vec_or(vec_t *dest, vec_t *src) {
    int i, j;
    for (i = 0; i < SIZE_A * MAX_B; ++i) {
        dest->c[i] |= src->c[i];
    }
}
inline vec_t *neighbour(int i, int j) {
    return &nb[j * a + i];
}

double difftime(clock_t t0, clock_t t1) {
    return ((double)t1 - t0) / clk_tck;
}

int curtime(void) {
    struct tms t;
    times(&t);
    return (int) t.tms_utime;
}

int line_count = 0;
void clear_line(void) {
    while (line_count > 0) {
        --line_count;
        printf("\x08 \x08");
    }
}
void unclear_line(int count) {
    line_count += count;
}
void report(void) {
    int output = 0, i, j;
    clear_line();
    output += printf("  (%.2fs) %llu -", difftime(gtime, curtime()), count);
    for (i = 0; i < a; ++i) {
        output += printf(" ");
        for (j = 0; j < b; ++j) {
            output += printf("%c",
                vec_test(&set, i, j) ? '1'
                : vec_test(&unset, i, j) ? '0' : '.');
        }
    }
    fflush(stdout);
    unclear_line(output);
}

void init(int a, int b) {
    int i, j;

    clk_tck = sysconf(_SC_CLK_TCK);
    gtime = curtime();

    ab = a * b;
    for (i = 0; i < a; ++i) {
        vec_set(&edge[0], i, 0);
        vec_set(&edge[1], i, b - 1);
    }
    for (i = 0; i < b; ++i) {
        vec_set(&edge[2], 0, i);
        vec_set(&edge[3], a - 1, i);
    }
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            vec_t *n = neighbour(i, j);
            if (i > 0)
                vec_set(n, i - 1, j);
            if (i < a - 1)
                vec_set(n, i + 1, j);
            if (j > 0)
                vec_set(n, i, j - 1);
            if (j < b - 1)
                vec_set(n, i, j + 1);
        }
    }
}

void check_solution(void) {
    int i, j;
    for (i = 0; i < 4; ++i) {
        int edge_seen = 0;
        for (j = 0; j < SIZE_A * MAX_B; ++j) {
            if ((set.c[j] & edge[i].c[j]) != 0) {
                edge_seen = 1;
                break;
            }
        }
        if (!edge_seen)
            return;
    }
    /* all edges seen, it's a solution */
    ++count;
    if ((count & 0x3ffffff) == 0) {
        report();
    }
}

void recurse(int i, int j) {
    vec_t old_surface, old_unset;
    int x, y;

    vec_set(&set, i, j);
    vec_copy(&old_unset, &unset);
    vec_copy(&old_surface, &surface);
    vec_or(&surface, neighbour(i, j));
    for (x = 0; x < a; ++x) {
        for (y = 0; y < b; ++y) {
            if (vec_test(&surface, x, y)
                && !vec_test(&set, x, y)
                && !vec_test(&unset, x, y)
            ) {
                recurse(x, y);
                vec_set(&unset, x, y);
            }
        }
    }
    check_solution();

    vec_copy(&surface, &old_surface);
    vec_copy(&unset, &old_unset);
    vec_unset(&set, i, j);
}

int main(int argc, char** argv) {
    int i;

    if (argc != 3) {
        fprintf(stderr, "usage: %s <a> <b>\n", argv[0]);
        return 1;
    }
    a = atoi(argv[1]);
    b = atoi(argv[2]);

    init(a, b);

    for (i = 0; i < a; ++i) {
        recurse(i, 0);
        vec_set(&unset, i, 0);
    }
    clear_line();
    printf("(%.2fs) f(%d, %d) = %lld\n", difftime(gtime, curtime()), a, b, count);
    return 0;
}

