#include <sys/times.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

int clk_tck;
int gtime;

int a, b;
typedef unsigned long long count_t;
count_t count;

#define MAX_A 8
#define MAX_B 8
#define MAX_AB (MAX_A * MAX_B)
typedef unsigned long long vec_t;

vec_t edge[4];
vec_t nb[MAX_AB];
vec_t fullvec;

typedef struct state_s {
    vec_t set;
    vec_t unset;
    vec_t surface;
    vec_t mask;
} state_t;
state_t state[MAX_AB];

#define v1 ((vec_t)1)
#define voff(x, y) ((y) * a + (x))
#define vbitx(x) (v1 << (x))
#define vbit(x, y) vbitx(voff(x, y))

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
void report(int index) {
    int output = 0, i, j;
    clear_line();
    output += printf("  (%.2fs) %llu -", difftime(gtime, curtime()), count);
    for (i = 0; i < a; ++i) {
        output += printf(" ");
        for (j = 0; j < b; ++j) {
            output += printf("%c",
                (state[index].set & vbit(i, j)) ? '1'
                : (state[index].unset * vbit(i, j)) ? '0' : '.');
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

void record_solution(count_t found, int index) {
    count += found;

#ifdef DEBUG
    report(index);
    keep_line();
#else
    if ((count & 0x3ffffff) < found)
        report(index);
#endif
}

inline int bit_count(vec_t v) {
    int i = 0;
    while (v) {
        v &= v - 1;
        ++i;
    }
    return i;
}

void recurse(int index, int i) {
    int x;
    state_t *s = &state[index];

    *s = state[index - 1];
    s->set |= vbitx(i);
    s->surface |= nb[i];
    s->mask = s->surface & ~(s->set | s->unset);

    /* Choice of starting position guarantees we touch edge[0] */
    if ((s->set & edge[1]) && (s->set & edge[2]) && (s->set & edge[3])) {
        /* all edges seen, it's a solution */

        if ((s->surface | s->set | s->unset) == fullvec) {
            /* all remaining permutations will be solutions: count the
             * bits in the mask representing the remaining surface, and
             * report 2 to that power.
             */
            record_solution(1 << bit_count(s->mask), index);
            return;
        }

        record_solution(1, index);
    }

    for (x = 0; x < a * b; ++x) {
        if (s->mask & vbitx(x)) {
            recurse(index + 1, x);
            s->unset |= vbitx(x);
            s->mask &= ~vbitx(x);
            if (!s->mask)
                return;
        }
    }
}

int main(int argc, char** argv) {
    int i;
    state_t *s = &state[0];

    if (argc != 3) {
        fprintf(stderr, "usage: %s <a> <b>\n", argv[0]);
        return 1;
    }
    a = atoi(argv[1]);
    b = atoi(argv[2]);

    init(a, b);

    for (i = 0; i < a; ++i) {
        recurse(1, voff(i, 0));
        s->unset |= vbit(i, 0);
    }
    clear_line();
    printf("(%.2fs) f(%d, %d) = %llu\n",
            difftime(gtime, curtime()), a, b, count);
    return 0;
}

