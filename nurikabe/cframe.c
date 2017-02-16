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
    int off;
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
void report(state_t *s) {
    int output = 0, i, j;

    clear_line();
    output += printf("  (%.2fs) %llu -", difftime(gtime, curtime()), count);
    for (j = 0; j < b; ++j) {
        output += printf(" ");
        for (i = 0; i < a; ++i) {
            output += printf("%c",
                (s->set & vbit(i, j)) ? '1'
                : (s->unset * vbit(i, j)) ? '0' : '.');
        }
    }
    fflush(stdout);
    unclear_line(output);
}

void init(void) {
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

void record_solution(count_t found, state_t *s) {
    count += found;

#ifdef DEBUG
    report(s);
    keep_line();
#else
    if ((count & 0x3ffffff) < found)
        report(s);
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

void calc(void) {
    int si;

    state[0].set = (vec_t)0;
    state[0].unset = (vec_t)0;
    state[0].surface = (vec_t)0;
    state[0].mask = edge[0];
    state[0].off = -1;
    si = 0;

    while (si >= 0) {
        state_t *scur = &state[si], *snext;
        int off = scur->off;
        vec_t mask = scur->mask;

        if (!mask) {
            /* nothing left to try, so derecurse */
            --si;
            continue;
        }

        /* is this the first attempt at this level? */
        if (off >= 0) {
            /* no - make sure we don't try again what we've just tried */
            scur->unset |= vbitx(off);
        }

        /* find the next thing to try at this level */
        while ((mask & vbitx(++off)) == 0)
            ;

        /* found something to try, so update this level and set up the next */
        ++si;
        snext = &state[si];
        snext->set = scur->set | vbitx(off);
        snext->unset = scur->unset;
        snext->surface = scur->surface | nb[off];
        snext->mask = snext->surface & ~(snext->set | snext->unset);
        snext->off = -1;

        scur->off = off;
        scur->mask &= ~vbitx(off); 
        scur->unset |= vbitx(off);  /* must be _after_ snext->unset init */

        if ((snext->set & edge[1])
            && (snext->set & edge[2])
            && (snext->set & edge[3])
        ) {
            /* all edges seen, it's a solution */

            if ((snext->surface | snext->set | snext->unset) == fullvec) {
                /* all remaining permutations will be solutions: count the
                 * bits in the mask representing the remaining surface, and
                 * report 2 to that power.
                 */
                record_solution(1 << bit_count(snext->mask), snext);
                --si;
                continue;
            }
            record_solution(1, snext);
        }
        /* now recurse with "next" as the new "cur" */
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

    init();
    calc();

    clear_line();
    printf("(%.2fs) f(%d, %d) = %llu\n",
            difftime(gtime, curtime()), a, b, count);
    return 0;
}

