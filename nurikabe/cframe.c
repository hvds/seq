#include <sys/times.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

int clk_tck;
int gtime;

int a, b;
unsigned long long count;

#define MAX_A 8
#define MAX_B 8
typedef unsigned long long vec_t;

vec_t edge[4];
vec_t nb[MAX_A * MAX_B];
vec_t fullvec;

vec_t set, unset, surface;

#define v1 ((vec_t)1)
#define voff(x, y) ((y) * MAX_A + (x))
#define vbit(x, y) (v1 << voff(x, y))

inline vec_t *neighbour(int i, int j) {
    return &nb[voff(i, j)];
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
void keep_line(void) {
    line_count = 0;
    printf("\n");
}
void report(void) {
    int output = 0, i, j;
    clear_line();
    output += printf("  (%.2fs) %llu -", difftime(gtime, curtime()), count);
    for (i = 0; i < a; ++i) {
        output += printf(" ");
        for (j = 0; j < b; ++j) {
            output += printf("%c",
                (set & vbit(i, j)) ? '1'
                : (unset * vbit(i, j)) ? '0' : '.');
        }
    }
    fflush(stdout);
    unclear_line(output);
}

void init(int a, int b) {
    int i, j;

    clk_tck = sysconf(_SC_CLK_TCK);
    gtime = curtime();

    for (i = 0; i < a; ++i) {
        edge[0] |= vbit(i, 0);
        edge[1] |= vbit(i, b - 1);
    }
    for (i = 0; i < b; ++i) {
        edge[2] |= vbit(0, i);
        edge[3] |= vbit(a - 1, i);
    }
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            vec_t *n = neighbour(i, j);
            if (i > 0)
                *n |= vbit(i - 1, j);
            if (i < a - 1)
                *n |= vbit(i + 1, j);
            if (j > 0)
                *n |= vbit(i, j - 1);
            if (j < b - 1)
                *n |= vbit(i, j + 1);
            fullvec |= vbit(i, j);
        }
    }
}

void record_solution(unsigned long long found) {
    count += found;

#ifdef DEBUG
    report();
    keep_line();
#else
    if ((count & 0x3ffffff) < found)
        report();
#endif
}

void recurse(int i, int j) {
    vec_t old_surface, old_unset, mask;
    int x, y;

    old_unset = unset;
    old_surface = surface;

    set |= vbit(i, j);
    surface |= *neighbour(i, j);
    mask = surface & ~(set | unset);

    /* Choice of starting position guarantees we touch edge[0] */
    if ((set & edge[1]) && (set & edge[2]) && (set & edge[3])) {
        /* all edges seen, it's a solution */

        if ((surface | set | unset) == fullvec) {
            /* all remaining permutations will be solutions: count the
             * bits in the mask representing the remaining surface, and
             * report 2 to that power.
             */
            vec_t var = mask;
            int bits = 0;
            while (var) {
                var &= var - 1;
                ++bits;
            }
            record_solution(1 << bits);
            goto done;
        }

        record_solution(1);
    }

    for (x = 0; x < a; ++x) {
        for (y = 0; y < b; ++y) {
            if (mask & vbit(x, y)) {
                recurse(x, y);
                unset |= vbit(x, y);
                mask &= ~vbit(x, y);
                if (!mask)
                    goto done;
            }
        }
    }

  done:
    surface = old_surface;
    unset = old_unset;
    set &= ~vbit(i, j);
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
        unset |= vbit(i, 0);
    }
    clear_line();
    printf("(%.2fs) f(%d, %d) = %lld\n", difftime(gtime, curtime()), a, b, count);
    return 0;
}

