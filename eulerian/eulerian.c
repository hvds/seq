#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "diag.h"

typedef unsigned int uint;
typedef unsigned long ulong;

ulong count;
ulong count_try;
uint n;
uint halfn;
uint v;
uint e;

typedef enum {
    e_in = 0,
    e_out
} e_direction;

uint *vectors = NULL;

typedef struct s_vertex {
    uint allocated;
    uint in;
    uint out;
} t_vertex;
t_vertex *vertices = NULL;

/* edge from left to right, left < right */
typedef struct s_edge {
    uint left;
    uint right;
    uint allocated;
} t_edge;
t_edge *edges = NULL;

typedef struct s_stack {
    uint edge;
    uint direction;
    uint more;
} t_stack;
t_stack *stack;
uint si;

double diag_delay = 1, diagt = 1;
double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

static inline void diag_plain(void) {
    double t1 = utime();
    char buf[2 * e + 1];
    uint j = 0;
    for (uint i = 0; i < si; ++i) {
        if (stack[i].more)
            buf[j++] = '*';
        buf[j++] = stack[i].direction ? '1' : '0';
    }
    buf[j] = 0;
    diag(buf);
    diagt = t1 + diag_delay;
    if (t1 > 60) {
        keep_diag();
        exit(1);
    }
}

static inline uint set_edge(uint ei, e_direction d, uint more) {
    t_stack *sp = &stack[si];
    t_edge *ep = &edges[ei];
    uint vector = ep->left ^ ep->right;
    t_vertex *lp = &vertices[ep->left];
    t_vertex *rp = &vertices[ep->right];

    /* allocate */
    sp->edge = ei;
    sp->direction = d;
    sp->more = more;
    ++si;

    ep->allocated = 1;
    lp->allocated |= vector;
    rp->allocated |= vector;
    if (d == e_in) {
        ++lp->in;
        ++rp->out;
        if (lp->in > halfn || rp->out > halfn)
            return 0;
    } else {
        ++lp->out;
        ++rp->in;
        if (lp->out > halfn || rp->in > halfn)
            return 0;
    }
    return 1;
}

static inline void unset_edge(uint ei, e_direction d) {
    t_edge *ep = &edges[ei];
    uint vector = ep->left ^ ep->right;
    t_vertex *lp = &vertices[ep->left];
    t_vertex *rp = &vertices[ep->right];

    ep->allocated = 0;
    lp->allocated &= ~vector;
    rp->allocated &= ~vector;
    if (d == e_in) {
        --lp->in;
        --rp->out;
    } else {
        --lp->out;
        --rp->in;
    }
}

void recurse(void) {
    while (1) {
        ++count_try;
        if (utime() >= diagt)
            diag_plain();
        /* go deeper */
        {
            /* find first unallocated edge */
            uint ei = si;

            /* if all edges were allocated, we have a solution */
            if (ei >= e) {
                ++count;
                goto destack;
            }

            t_edge *ep = &edges[ei];
            t_vertex *lp = &vertices[ep->left];
            t_vertex *rp = &vertices[ep->right];
            uint more = 0;
            e_direction d;

            if (lp->in == halfn) {
                /* forced out */
                if (rp->in == halfn)
                    goto destack; /* but that's not possible */
                d = e_out;
            } else if (lp->out == halfn) {
                /* forced in */
                if (rp->out == halfn)
                    goto destack; /* but that's not possible */
                d = e_in;
            } else {
                /* both are possible */
                more = 1;
                d = e_in;
            }
            if (set_edge(ei, d, more))
                continue;   /* recurse deeper */
            goto destack;   /* or fail */
        }
      destack:
        {
            if (si == 0)
                break;
            --si;
            t_stack *sp = &stack[si];
            unset_edge(sp->edge, sp->direction);
            if (sp->more && set_edge(sp->edge, e_out, 0))
                continue;   /* recurse deeper */
            goto destack;
        }
    }
}

void run(void) {
    count = 0UL;
    count_try = 0UL;
    si = 0;
    for (uint i = 0; i < n; ++i)
        set_edge(i, (i < halfn) ? e_in : e_out, 0);
    uint mult = 1;
    for (uint i = n; i > 0; --i) {
        if (i > halfn)
            mult *= i;
        else
            mult /= i;
    }
    recurse();
    keep_diag();
    printf("For n=%u found %lu eulerian graphs (recurse %lu)\n",
            n, count * mult, count_try);
}

void init(void) {
    /* number of vertices is 2^n; number of edges is n 2^{n-1} */
    v = 1 << n;
    e = n * (1 << (n - 1));
    halfn = n / 2;

    /* for vertices i, j, there is an edge between them if i ^ j in
     * this list of vectors */
    vectors = (uint *)malloc(n * sizeof(uint));
    for (uint i = 0; i < n; ++i)
        vectors[i] = 1 << i;

    /* for each vertex we track the count of +ve edges allocated,
     * and the set of vectors allocated */
    vertices = (t_vertex *)calloc(v, sizeof(t_vertex));

    edges = (t_edge *)malloc(e * sizeof(t_edge));
    uint ec = 0;
    for (uint i = 0; i < v; ++i) {
        for (uint vj = 0; vj < n; ++vj) {
            uint j = i ^ vectors[vj];
            if (j < i)
                continue;
            edges[ec].left = i;
            edges[ec].right = j;
            edges[ec].allocated = 0;
            ++ec;
        }
    }

    stack = (t_stack *)malloc((e + 1) * sizeof(t_stack));
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    n = strtoul(argv[1], NULL, 10);
    if (n & 1) {
        fprintf(stderr, "<n> must be even, not %u\n", n);
        return 1;
    }
    if (n < 2) {
        fprintf(stderr, "<n> must be at least 2, not %u\n", n);
        return 1;
    }
    if (n > 8)
        fprintf (stderr, "n > 8 may cause problems, use at own risk\n");

    init_diag();
    init();
    run();
}
