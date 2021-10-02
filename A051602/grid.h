#ifndef GRID_H
#define GRID_H

#include <assert.h>
#include <strings.h>

typedef struct {
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int *grid;
} grid_t;

grid_t *new_grid(int xmin, int xmax, int ymin, int ymax);
void free_grid(grid_t *g);
void copy_grid(grid_t *dest, grid_t *src);
grid_t *dup_grid(grid_t *g);
void resize_grid(grid_t *g, int xmin, int xmax, int ymin, int ymax);

static int grid_xspan(grid_t *g) {
    return g->xmax - g->xmin + 1;
}

static int grid_yspan(grid_t *g) {
    return g->ymax - g->ymin + 1;
}

static int grid_size(grid_t *g) {
    return grid_xspan(g) * grid_yspan(g);
}

static int grid_off(grid_t *g, int x, int y) {
    int off = (x - g->xmin) + (y - g->ymin) * grid_xspan(g);
    assert(off >= 0 && off < grid_size(g));
    return off;
}

static int get_grid(grid_t *g, int x, int y) {
    if (x < g->xmin || x > g->xmax || y < g->ymin || y > g->ymax)
        return 0;
    return g->grid[grid_off(g, x, y)];
}
    
static int set_grid(grid_t *g, int x, int y, int val) {
    if (x < g->xmin || x > g->xmax || y < g->ymin || y > g->ymax) {
        int xmin = (x < g->xmin) ? x : g->xmin;
        int xmax = (x > g->xmax) ? x : g->xmax;
        int ymin = (y < g->ymin) ? y : g->ymin;
        int ymax = (y > g->ymax) ? y : g->ymax;
        resize_grid(g, xmin, xmax, ymin, ymax);
    }
    g->grid[grid_off(g, x, y)] = val;
}

static void reset_grid(grid_t *g) {
    bzero(g->grid, sizeof(int) * grid_xspan(g) * grid_yspan(g));
}

#endif
