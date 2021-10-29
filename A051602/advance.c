#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/times.h>
#include <unistd.h>

#include "loc.h"
#include "sym.h"
#include "grid.h"

long clock_tick;

typedef unsigned char uchar;

typedef struct {
    size_t size;
    size_t used;
    grid_t *arr;
} collection_t;

typedef struct {
    int n;          /* number of points */
    int squares;    /* best count of squares found */
    /* structures found with (squares - diff) .. (squares) squares */
    collection_t *coll;
} found_t;

int n;              /* We're trying to find A051602(i) up to n */
/* For A051602(i) we'll try to extend structures found for A051602(i-1)
 * and A051602(i-2) that are within <diff> of the max found.
 */
int diff;
int verbose = 0;    /* report every maximum (1) or iteration (2) */
int quiet = 0;      /* skip reporting for arrangements with < n points */
int best;           /* Greatest number of squares seen in any arrangement */
found_t found[3];  /* Best structures found */

/* Reusable structure to hold symmetries for in_coll() */
struct {
    size_t size;
    uchar *grids;
} in_coll_sym;

unsigned long visit = 0;        /* Count of iterations */

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

/* Show a result-so-far, consisting of i points */
void report(int i) {
    found_t *f = &found[2];
    printf("%d %d:", f->n, f->squares);
    for (int i = diff; i >= 0; --i)
        printf(" %ld", f->coll[i].used);
    printf(" [%.2fs]\n", timing());
}

void report_grid(int points, int squares, grid_t *g) {
    loc_t p;

    printf("(%lu) p%d %dx%d %d:", visit, points, g->span.x, g->span.y, squares);
    for (p.x = 0; p.x < g->span.x; ++p.x)
        for (p.y = 0; p.y < g->span.y; ++p.y)
            if (is_grid(g, p))
                printf(" %d:%d", p.x, p.y);
    printf("\n");
}

void clear_coll(collection_t *c) {
    for (int i = 0; i < c->used; ++i)
        free(c->arr[i].grid);
    c->used = 0;
}

static inline uchar *sym_grid(int i) {
    return &in_coll_sym.grids[in_coll_sym.size * i];
}

static inline void sym_swap_x(loc_t span, uchar *gs, uchar *gd) {
    size_t bytes_per_col = (span.y + 7) >> 3;
    for (int i = 0; i < span.x; ++i)
        memcpy(gd + bytes_per_col * (span.x - 1 - i),
                gs + bytes_per_col * i, bytes_per_col);
}

static inline void sym_swap_y(loc_t span, uchar *gs, uchar *gd) {
    /* FIXME: we can make this faster, eg with 8-bit reversing lookup */
    /* FIXME: would prefer to use grid.h methods */
    size_t bytes_per_col = (span.y + 7) >> 3;
    memset(gd, 0, span.x * bytes_per_col);
    for (int x = 0; x < span.x; ++x)
        for (int y = 0; y < span.y; ++y) {
            uchar *gsx = gs + x * bytes_per_col + (y >> 3);
            uchar *gdx = gd + x * bytes_per_col + ((span.y - 1 - y) >> 3);
            if (*gsx & (1 << (y & 7)))
                *gdx |= 1 << ((span.y - 1 - y) & 7);
        }
}

static inline void sym_trans(loc_t span, uchar *gs, uchar *gd) {
    size_t sbytes_per_col = (span.y + 7) >> 3;
    size_t dbytes_per_col = (span.x + 7) >> 3;
    memset(gd, 0, span.y * dbytes_per_col);
    for (int x = 0; x < span.x; ++x)
        for (int y = 0; y < span.y; ++y) {
            uchar *gsx = gs + x * sbytes_per_col + (y >> 3);
            if (*gsx & (1 << (y & 7))) {
                uchar *gdx = gd + y * dbytes_per_col + (x >> 3);
                *gdx |= (1 << (x & 7));
            }
        }
}

int in_coll(collection_t *c, grid_t *g) {
    grid_t gsym = (grid_t){ g->span, NULL };
    grid_t gtrans = (grid_t){ (loc_t){ gsym.span.y, gsym.span.x }, NULL };
    size_t ssym = grid_size(&gsym);
    size_t strans = grid_size(&gtrans);
    size_t size = (ssym > strans) ? ssym : strans;
    if (size > in_coll_sym.size) {
        free(in_coll_sym.grids);
        in_coll_sym.grids = (uchar *)malloc(size * 8);
        in_coll_sym.size = size;
    }

    /* 0: null symmetry */
    memcpy(sym_grid(0), g->grid, ssym);
    /* 1: x -> -x */
    sym_swap_x(gsym.span, sym_grid(0), sym_grid(1));
    /* 2: y -> -y */
    sym_swap_y(gsym.span, sym_grid(0), sym_grid(2));
    /* 3: x -> -x, y -> -y */
    sym_swap_x(gsym.span, sym_grid(2), sym_grid(3));
    /* 4: x -> y, y -> x */
    sym_trans(gsym.span, sym_grid(0), sym_grid(4));
    /* 5: x -> -y, y -> x */
    sym_swap_x(gtrans.span, sym_grid(4), sym_grid(5));
    /* 6: x -> y, y -> -x */
    sym_swap_y(gtrans.span, sym_grid(4), sym_grid(6));
    /* 7: x -> -y, y -> -x */
    sym_swap_x(gtrans.span, sym_grid(6), sym_grid(7));

    /* FIXME: use a hash */
    for (int i = 0; i < c->used; ++i) {
        grid_t *g2 = &c->arr[i];
        if (loc_eq(g2->span, gsym.span))
            for (int j = 0; j <= 3; ++j) {
                gsym.grid = sym_grid(j);
                if (same_grid(g2, &gsym))
                    return 1;
            }
        if (loc_eq(g2->span, gtrans.span))
            for (int j = 4; j <= 7; ++j) {
                gtrans.grid = sym_grid(j);
                if (same_grid(g2, &gtrans))
                    return 1;
            }
    }
    return 0;
}

int save_coll(collection_t *c, grid_t *g) {
    if (in_coll(c, g))
        return 0;
    if (c->used >= c->size) {
        size_t new_size = c->size ? c->size * 3 / 2 : 20;
        c->arr = (grid_t *)realloc(c->arr, new_size * sizeof(grid_t));
        memset(&c->arr[c->size], 0, new_size - c->size);
        c->size = new_size;
    }
    c->arr[c->used++] = *g;
    free(g);
    return 1;
}

void rotate_found(void) {
    found_t f;

    f = found[0];
    memmove(&found[0], &found[1], sizeof(found_t) * 2);
    found[2] = f;
    found[2].n += 3;
    found[2].squares = found[1].squares;
    for (int i = 0; i <= diff; ++i)
        clear_coll(&found[2].coll[i]);
}

void advance_found(int squares) {
    found_t *f = &found[2];
    int off = squares - f->squares;
    if (off <= diff) {
        for (int i = 0; i <= diff - off; ++i) {
            collection_t t = f->coll[i];
            f->coll[i] = f->coll[i + off];
            f->coll[i + off] = t;
        }
        for (int i = diff + 1 - off; i <= diff; ++i)
            clear_coll(&f->coll[i]);
    } else {
        for (int i = 0; i <= diff; ++i)
            clear_coll(&f->coll[i]);
    }
    f->squares = squares;
}

int save_found(grid_t *g, int squares) {
    found_t *f = &found[2];
    if (squares < f->squares - diff)
        return 0;
    if (squares > f->squares)
        advance_found(squares);

    collection_t *c = &f->coll[squares - (f->squares - diff)];
    if (save_coll(c, g)) {
        /* g has been freed */
        if (verbose && (verbose > 1 || c->used == 1)
            && (!quiet || f->n == n)
        )
            report_grid(f->n, squares, &c->arr[c->used - 1]);
        ++visit;
        return 1;
    } else
        return 0;
}

void init(void) {
    /* used to calculate timings */
    clock_tick = sysconf(_SC_CLK_TCK);

    /* for less-verbose output don't buffer, even if stdout is not a tty */
    if (verbose != 2)
        setvbuf(stdout, (char*)NULL, _IONBF, 0);

    in_coll_sym.size = 0;
    in_coll_sym.grids = (uchar *)NULL;

    for (int i = 0; i < 3; ++i) {
        found[i].n = i + 2;
        found[i].squares = (i == 2) ? 1 : 0;
        found[i].coll = (collection_t *)calloc(diff + 1, sizeof(collection_t));
    }
    grid_t *g0 = new_grid((loc_t){ 2, 2 });
    set_grid(g0, (loc_t){ 0, 0 });
    set_grid(g0, (loc_t){ 0, 1 });
    set_grid(g0, (loc_t){ 1, 0 });
    set_grid(g0, (loc_t){ 1, 1 });
    save_found(g0, 1);
}

void finish(void) {
    for (int i = 0; i < 3; ++i) {
        collection_t *c = found[i].coll;
        for (int j = 0; j <= diff; ++j) {
            clear_coll(&c[j]);
            free(c[j].arr);
        }
        free(c);
    }
    free(in_coll_sym.grids);
}

int find_squares(grid_t *g, loc_t p0) {
    int squares = 0;
    loc_t p1, p2, p3;

    for (p1.x = 0; p1.x < g->span.x; ++p1.x)
        for (p1.y = 0; p1.y < g->span.y; ++p1.y) {
            if (!is_grid(g, p1))
                continue;
            p2 = loc_rot90(p0, loc_diff(p1, p0));
            if (loc_lt(p1, p2) && is_grid(g, p2)) {
                p3 = loc_rot270(p1, loc_diff(p0, p1));
                if (is_grid(g, p3))
                    ++squares;
            }
            p2 = loc_rot270(p0, loc_diff(p1, p0));
            if (loc_lt(p1, p2) && is_grid(g, p2)) {
                p3 = loc_rot90(p1, loc_diff(p0, p1));
                if (is_grid(g, p3))
                    ++squares;
            }
        }
    return squares;
}

void try1(grid_t *gs, loc_t p, int ss) {
    loc_t off = { p.x < 0 ? -p.x : 0, p.y < 0 ? -p.y : 0 };
    loc_t span = {
        off.x ? gs->span.x + off.x : (p.x >= gs->span.x) ? p.x + 1 : gs->span.x,
        off.y ? gs->span.y + off.y : (p.y >= gs->span.y) ? p.y + 1 : gs->span.y
    };
    grid_t *gd = new_grid(span);
    loc_t q;

    for (q.x = 0; q.x < gs->span.x; ++q.x)
        for (q.y = 0; q.y < gs->span.y; ++q.y)
            if (is_grid(gs, q))
                set_grid(gd, loc_sum(q, off));

    p = loc_sum(p, off);

    int snew = ss + find_squares(gd, p);
    set_grid(gd, p);
    if (!save_found(gd, snew))
        free_grid(gd);
}

void try2(grid_t *gs, loc_t p, loc_t q, int ss) {
    loc_t off = { 0, 0 };
    if (p.x < -off.x) off.x = -p.x;
    if (q.x < -off.x) off.x = -q.x;
    if (p.y < -off.y) off.y = -p.y;
    if (q.y < -off.y) off.y = -q.y;

    loc_t span = gs->span;
    if (p.x >= span.x) span.x = p.x + 1;
    if (q.x >= span.x) span.x = q.x + 1;
    if (p.y >= span.y) span.y = p.y + 1;
    if (q.y >= span.y) span.y = q.y + 1;
    span = loc_sum(span, off);

    grid_t *gd = new_grid(span);
    loc_t r;
    for (r.x = 0; r.x < gs->span.x; ++r.x)
        for (r.y = 0; r.y < gs->span.y; ++r.y)
            if (is_grid(gs, r))
                set_grid(gd, loc_sum(r, off));
    p = loc_sum(p, off);
    q = loc_sum(q, off);

    int snew = ss + find_squares(gd, p);
    set_grid(gd, p);
    snew += find_squares(gd, q);
    set_grid(gd, q);
    if (!save_found(gd, snew))
        free_grid(gd);
}

void try_all1(grid_t *g, int ss) {
    loc_t p0, p1, p2, p3;
    int b2, b3;

    for (p0.x = 0; p0.x < g->span.x; ++p0.x) {
        for (p0.y = 0; p0.y < g->span.y; ++p0.y) {
            if (!is_grid(g, p0))
                continue;
            for (p1.x = p0.x; p1.x < g->span.x; ++p1.x) {
                int p1ybase = (p1.x == p0.x) ? p0.y + 1 : 0;
                for (p1.y = p1ybase; p1.y < g->span.y; ++p1.y) {
                    if (!is_grid(g, p1))
                        continue;
                    p2 = loc_rot90(p0, loc_diff(p1, p0));
                    b2 = is_grid(g, p2);
                    p3 = loc_rot270(p1, loc_diff(p0, p1));
                    b3 = is_grid(g, p3);
                    if (b2 && !b3) {
                        try1(g, p3, ss);
                    } else if (b3 && !b2) {
                        try1(g, p2, ss);
                    }
                    p2 = loc_rot270(p0, loc_diff(p1, p0));
                    b2 = is_grid(g, p2);
                    p3 = loc_rot90(p1, loc_diff(p0, p1));
                    b3 = is_grid(g, p3);
                    if (b2 && !b3) {
                        try1(g, p3, ss);
                    } else if (b3 && !b2) {
                        try1(g, p2, ss);
                    }
                }
            }
        }
    }
}

void try_all2(grid_t *g, int ss) {
    loc_t p0, p1, p2, p3;
    loc_t doffset;
    grid_t *dg = expand_grid(g, &doffset);
    loc_t dp0, dp1, dp2, dp3;
    int par0;

    for (p0.x = 0; p0.x < g->span.x; ++p0.x) {
        for (p0.y = 0; p0.y < g->span.y; ++p0.y) {
            if (!is_grid(g, p0))
                continue;
            dp0 = loc_sum(doffset, loc_expand(p0));
            par0 = loc_parity(p0);
            for (p1.x = p0.x; p1.x < g->span.x; ++p1.x) {
                int p1ybase = (p1.x == p0.x) ? p0.y + 1 : 0;
                for (p1.y = p1ybase; p1.y < g->span.y; ++p1.y) {
                    if (!is_grid(g, p1))
                        continue;

                    p2 = loc_rot90(p0, loc_diff(p1, p0));
                    p3 = loc_rot270(p1, loc_diff(p0, p1));
                    if (!is_grid(g, p2) && !is_grid(g, p3))
                        try2(g, p2, p3, ss);

                    p2 = loc_rot270(p0, loc_diff(p1, p0));
                    p3 = loc_rot90(p1, loc_diff(p0, p1));
                    if (!is_grid(g, p2) && !is_grid(g, p3))
                        try2(g, p2, p3, ss);

                    /* use doubled points only if needed */
                    if (loc_parity(p1) == par0) {
                        p2 = loc_diag1(p0, p1);
                        p3 = loc_diag2(p0, p1);
                        if (!is_grid(g, p2) && !is_grid(g, p3))
                            try2(g, p2, p3, ss);
                    } else {
                        dp1 = loc_sum(doffset, loc_expand(p1));
                        dp2 = loc_diag1(dp0, dp1);
                        dp3 = loc_diag2(dp0, dp1);
                        if (!is_grid(dg, dp2) && !is_grid(dg, dp3))
                            try2(dg, dp2, dp3, ss);
                    }
                }
            }
        }
    }
    free_grid(dg);
}

void try_coll1(int j) {
    found_t *fs = &found[1];
    collection_t *cs = &fs->coll[diff - j];
    int ss = fs->squares - j;

    for (int k = 0; k < cs->used; ++k)
        try_all1(&cs->arr[k], ss);
}

void try_coll2(int j) {
    found_t *fs = &found[0];
    collection_t *cs = &fs->coll[diff - j];
    int ss = fs->squares - j;

    for (int k = 0; k < cs->used; ++k)
        try_all2(&cs->arr[k], ss);
}

/* Report current state and advance to next state. */
void advance(int i) {
    report(i);
    if (i < n) {
        rotate_found();
        for (int j = 0; j <= diff; ++j) {
            try_coll1(j);
            try_coll2(j);
        }
    }
}

int main(int argc, char** argv) {
    int arg = 1;

    while (arg < argc && argv[arg][0] == '-') {
        char *s = argv[arg++];
        if (strcmp("--", s) == 0)
            break;
        else if (strcmp("-m", s) == 0)
            verbose = 1;
        else if (strcmp("-v", s) == 0)
            verbose = 2;
        else if (strcmp("-q", s) == 0)
            quiet = 1;
        else {
            fprintf(stderr, "Unknown option '%s'\n", s);
            exit(1);
        }
    }
    if (arg + 2 != argc) {
        fprintf(stderr, "Usage: try [-v | -m] [-q] <n> <diff>\n");
        return 1;
    }

    n = atoi(argv[arg++]);
    if (n < 4) {
        fprintf(stderr, "Need n >= 4\n");
        return 1;
    }
    diff = atoi(argv[arg++]);
    if (diff < 0) {
        fprintf(stderr, "Need diff >= 0\n");
        return 1;
    }

    init();
    for (int i = 4; i <= n; ++i) {
        advance(i);
    }

    finish();
    return 0;
}
