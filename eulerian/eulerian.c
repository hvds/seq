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

void diag_plain(void) {
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

void try_set(uint ei, uint d);
void try_next(void) {
    /* find first unallocated edge */
    uint ei = si;

    /* if all edges were allocated, we have a solution */
    if (ei >= e) {
        ++count;
        return;
    }

    t_edge *ep = &edges[ei];
    t_vertex *lp = &vertices[ep->left];
    t_vertex *rp = &vertices[ep->right];

    if (lp->in == halfn) {
        /* forced out */
        if (rp->in == halfn)
            return; /* but that's not possible */
        stack[si].more = 0;
        try_set(ei, e_out);
        return;
    }
    if (lp->out == halfn) {
        /* forced in */
        if (rp->out == halfn)
            return; /* but that's not possible */
        stack[si].more = 0;
        try_set(ei, e_in);
        return;
    }
    /* both are possible */
    stack[si].more = 1;
    try_set(ei, e_in);
    stack[si].more = 0;
    try_set(ei, e_out);
}

void try_set(uint ei, uint d) {
    ++count_try;
    if (utime() >= diagt)
        diag_plain();
    uint tnsi = si;
    {
        t_stack *sp = &stack[si];
        t_edge *ep = &edges[ei];
        uint vector = ep->left ^ ep->right;
        t_vertex *lp = &vertices[ep->left];
        t_vertex *rp = &vertices[ep->right];

        /* allocate */
        sp->edge = ei;
        sp->direction = d;
        ++si;

        ep->allocated = 1;
        lp->allocated |= vector;
        rp->allocated |= vector;
        if (d == e_in) {
            ++lp->in;
            ++rp->out;
            if (lp->in > halfn || rp->out > halfn)
                goto destack;
        } else {
            ++lp->out;
            ++rp->in;
            if (lp->out > halfn || rp->in > halfn)
                goto destack;
        }
    }

    /* find the next edge to try, and recurse until done */
    try_next();

  destack:
    while (si > tnsi) {
        --si;
        t_stack *sp = &stack[si];
        t_edge *ep = &edges[sp->edge];
        uint vector = ep->left ^ ep->right;
        t_vertex *lp = &vertices[ep->left];
        t_vertex *rp = &vertices[ep->right];

        ep->allocated = 0;
        lp->allocated &= ~vector;
        rp->allocated &= ~vector;
        if (sp->direction == e_in) {
            --lp->in;
            --rp->out;
        } else {
            --lp->out;
            --rp->in;
        }
    }
}

void run(void) {
    count = 0UL;
    count_try = 0UL;
    si = 0;
    try_set(0, 1);
    keep_diag();
    printf("For n=%u found %lu eulerian graphs (recurse %lu)\n",
            n, count * 2, count_try);
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
