#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
char *available = NULL;
uint *choice = NULL;
uint ci;

typedef struct s_vertex {
    uint in;
    uint out;
    uint *edges;
} t_vertex;
t_vertex *vertices = NULL;

/* edge from left to right, left < right */
typedef struct s_edge {
    uint left;
    uint right;
} t_edge;
t_edge *edges = NULL;

typedef struct s_alloc {
    uint edge;
    e_direction direction;
} t_alloc;
t_alloc *stack;
uint si;
t_alloc *pend = NULL;
uint pi;

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
    /* overkill worst case for d=6, will fail for d > 6 */
    char buf[192 * 4];
    uint j = 0;
    for (uint i = 0; i < ci; ++i)
        j += sprintf(&buf[j], "%u ", choice[i]);
    buf[j] = 0;
    diag(buf);
    diagt = t1 + diag_delay;
    if (t1 > 60) {
        keep_diag();
        exit(1);
    }
}

static inline void add_pend(uint v, e_direction d) {
    t_vertex *vp = &vertices[v];
    uint spare = halfn - ((d == e_in) ? vp->in : vp->out);
    for (uint i = 0; i < n; ++i) {
        uint ei = vp->edges[i];
        if (available[ei / 8] & (1 << (ei % 8))) {
            t_alloc *pp = &pend[pi++];
            pp->edge = ei;
            pp->direction = (edges[ei].left == v) ? d : 1 - d;
            --spare;
            if (spare == 0)
                return;
        }
    }
    fprintf(stderr, "out of edges\n");
    exit(1);
}

static inline uint set_edge(uint ei, e_direction d) {
    t_alloc *sp = &stack[si];
    t_edge *ep = &edges[ei];
    uint vector = ep->left ^ ep->right;
    t_vertex *lp = &vertices[ep->left];
    t_vertex *rp = &vertices[ep->right];

    /* allocate */
    sp->edge = ei;
    sp->direction = d;
    ++si;

    available[ei / 8] &= ~(1 << (ei % 8));
    if (d == e_in) {
        ++lp->in;
        ++rp->out;
        if (lp->in > halfn || rp->out > halfn)
            return 0;
        if (lp->in == halfn && lp->out < halfn)
            add_pend(ep->left, e_out);
        if (rp->out == halfn && rp->in < halfn)
            add_pend(ep->right, e_in);
    } else {
        ++lp->out;
        ++rp->in;
        if (lp->out > halfn || rp->in > halfn)
            return 0;
        if (lp->out == halfn && lp->in < halfn)
            add_pend(ep->left, e_in);
        if (rp->in == halfn && rp->out < halfn)
            add_pend(ep->right, e_out);
    }
    return 1;
}

static inline void unset_edge(uint ei, e_direction d) {
    t_edge *ep = &edges[ei];
    uint vector = ep->left ^ ep->right;
    t_vertex *lp = &vertices[ep->left];
    t_vertex *rp = &vertices[ep->right];

    available[ei / 8] |= 1 << (ei % 8);
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

        /* make forced moves first */
        while (pi) {
            t_alloc *pp = &pend[--pi];
            uint ei = pp->edge;
            if (!(available[ei / 8] & (1 << (ei % 8))))
                continue;
            if (!set_edge(ei, pp->direction))
                goto destack;
        }

        /* go deeper */
        {
            /* find first unallocated edge */
            uint ei;
            for (uint i = 0; i < e / 8; ++i)
                if (available[i])
                    for (uint j = 0; j < 8; ++j)
                        if (available[i] & (1 << j)) {
                            ei = i * 8 + j;
                            goto got_ei;
                        }

            /* if all edges were allocated, we have a solution */
            ++count;
            goto destack;

          got_ei:
            ;
            t_edge *ep = &edges[ei];
            t_vertex *lp = &vertices[ep->left];
            t_vertex *rp = &vertices[ep->right];
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
                choice[ci++] = si;
                d = e_in;
            }
            if (set_edge(ei, d))
                continue;   /* recurse deeper */
            goto destack;   /* or fail */
        }

      destack:
        {
            pi = 0; /* remaining forced moves no longer relevant */
            if (ci == 0)
                break;
            uint target = choice[--ci];
            while (si > target) {
                t_alloc *sp = &stack[--si];
                unset_edge(sp->edge, sp->direction);
            }
            if (set_edge(stack[si].edge, e_out))
                continue;
            goto destack;
        }
    }
}

void run(void) {
    count = 0UL;
    count_try = 0UL;
    si = 0;
    ci = 0;
    pi = 0;
    for (uint i = 0; i < n; ++i)
        set_edge(i, (i < halfn) ? e_in : e_out);
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
    for (uint i = 0; i < v; ++i)
        vertices[i].edges = (uint *)malloc(n * sizeof(uint));

    edges = (t_edge *)malloc(e * sizeof(t_edge));
    uint ec = 0;
    for (uint i = 0; i < v; ++i) {
        for (uint vj = 0; vj < n; ++vj) {
            uint j = i ^ vectors[vj];
            if (j < i)
                continue;
            edges[ec].left = i;
            edges[ec].right = j;
            vertices[i].edges[vertices[i].in++] = ec;
            if (vertices[i].in == n)
                vertices[i].in = 0;
            vertices[j].edges[vertices[j].in++] = ec;
            if (vertices[j].in == n)
                vertices[j].in = 0;
            ++ec;
        }
    }

    available = (char*)malloc((e + 7) / 8);
    memset(available, 0xff, (e + 7) / 8);
    for (uint i = e; i & 7; ++i)
        available[i / 8] &= ~(1 << (i % 8));

    stack = (t_alloc *)malloc((e + 1) * sizeof(t_alloc));
    choice = (uint *)malloc(e * sizeof(uint));
    pend = (t_alloc *)malloc(e * sizeof(t_alloc));
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
    if (n > 6) {
        fprintf(stderr, "For n > 6, must at least increase diag buf\n");
        return 1;
    }

    init_diag();
    init();
    run();
}
