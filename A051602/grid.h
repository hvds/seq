#ifndef GRID_H
#define GRID_H

#include <stdlib.h>
#include "loc.h"

typedef unsigned char uchar;

typedef struct {
    loc_t span;
    uchar *grid;
} grid_t;

inline static size_t grid_bytespercol(grid_t *g) {
    return (g->span.y + 7) >> 3;
}

inline static size_t grid_byte(grid_t *g, loc_t p) {
    return p.x * grid_bytespercol(g) + (p.y >> 3);
}

inline static size_t grid_size(grid_t *g) {
    return grid_byte(g, (loc_t){ g->span.x, 0 });
}

inline static uchar grid_off(grid_t *g, loc_t p) {
    return 1 << (p.y & 7);
}

inline static grid_t *new_grid(loc_t span) {
    grid_t *g = (grid_t *)malloc(sizeof(grid_t));
    g->span = span;
    g->grid = (uchar *)calloc(grid_size(g), 1);
    return g;
}

inline static void free_grid(grid_t *g) {
    free(g->grid);
    free(g);
}

inline static void set_grid(grid_t *g, loc_t p) {
    g->grid[grid_byte(g, p)] |= grid_off(g, p);
}

inline static int is_grid(grid_t *g, loc_t p) {
    if (p.x < 0 || p.y < 0 || p.x >= g->span.x || p.y >= g->span.y)
        return 0;
    return (g->grid[grid_byte(g, p)] & grid_off(g, p)) ? 1 : 0;
}

/*
 * loc_expand() all grid points: (x,y) -> (y+x,y-x)
 * we size to worst case, then adjust for observed span
 */
inline static grid_t *expand_grid(grid_t *g, loc_t *offset) {
    loc_t xspan = { g->span.x + g->span.y - 1, g->span.x + g->span.y - 1 };
    loc_t minspan = xspan, maxspan = { 0, 0 };
    loc_t off = { 0, g->span.x - 1 };
    loc_t p;
    grid_t *xg = new_grid(xspan);
    size_t col = grid_bytespercol(xg);

    for (p.x = 0; p.x < g->span.x; ++p.x)
        for (p.y = 0; p.y < g->span.y; ++p.y)
            if (is_grid(g, p)) {
                loc_t xp = loc_sum(off, loc_expand(p));
                set_grid(xg, xp);
                if (minspan.x > xp.x) minspan.x = xp.x;
                if (minspan.y > xp.y) minspan.y = xp.y;
                if (maxspan.x <= xp.x) maxspan.x = xp.x + 1;
                if (maxspan.y <= xp.y) maxspan.y = xp.y + 1;
            }
    if (minspan.x > 0) {
        size_t shift = minspan.x * col;
        memmove(&xg->grid[0], &xg->grid[shift], grid_size(xg) - shift);
        maxspan.x -= minspan.x;
    }
    if (minspan.y > 0) {
        int ybytes = minspan.y >> 3;
        int ybits = minspan.y & 7;
        for (int x = 0; x < maxspan.x; ++x) {
            uchar *s = &xg->grid[x * col];
            if (ybytes) {
                memmove(s, s + ybytes, col - ybytes);
                memset(s + col - ybytes, 0, ybytes);
            }
            if (ybits) {
                uchar shift = 0;
                for (int i = col - 1; i >= 0; --i) {
                    uchar tmp = s[i];
                    s[i] = (tmp >> ybits) | shift;
                    shift = tmp << (8 - ybits);
                }
            }
        }
    }
    size_t newcol = (maxspan.y - minspan.y + 7) >> 3;
    if (newcol < col) {
        for (int x = 1; x < maxspan.x; ++x)
            memmove(&xg->grid[x * newcol], &xg->grid[x * col], newcol);
        memset(&xg->grid[maxspan.x * newcol], 0, maxspan.x * (col - newcol));
    }
    if (maxspan.x < xspan.x)
        memset(&xg->grid[maxspan.x * col], 0, col * (xspan.x - maxspan.x));

    xg->span = (loc_t){ maxspan.x, maxspan.y - minspan.y };
    *offset = loc_diff(minspan, off);
    return xg;
}

inline static int same_grid(grid_t *g1, grid_t *g2) {
    if (g1->span.x != g2->span.x)
        return 0;
    if (g1->span.y != g2->span.y)
        return 0;
    return memcmp(g1->grid, g2->grid, grid_size(g1)) ? 0 : 1;
}

#endif
