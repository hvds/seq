#include <stdlib.h>
#include <string.h>
#include "grid.h"

grid_t *new_grid(int xmin, int xmax, int ymin, int ymax) {
    grid_t *g = (grid_t *)calloc(sizeof(grid_t), 1);
    g->xmin = xmin;
    g->xmax = xmax;
    g->ymin = ymin;
    g->ymax = ymax;
    g->grid = (int *)calloc(sizeof(int), grid_xspan(g) * grid_yspan(g));
    return g;
}

void free_grid(grid_t *g) {
    free(g->grid);
    free(g);
}

static void _copy_grid(grid_t *gd, grid_t *gs) {
    for (int x = gs->xmin; x <= gs->xmax; ++x) {
        for (int y = gs->ymin; y <= gs->ymax; ++y) {
            set_grid(gd, x, y, get_grid(gs, x, y));
        }
    }
}

void copy_grid(grid_t *gd, grid_t *gs) {
    resize_grid(gd, gs->xmin, gs->xmax, gs->ymin, gs->ymax);
    _copy_grid(gd, gs);
}

grid_t *dup_grid(grid_t *g) {
    grid_t *g2 = (grid_t *)malloc(sizeof(grid_t));

    memcpy(g2, g, sizeof(grid_t));
    g2->grid = (int *)malloc(sizeof(int) * grid_xspan(g) * grid_yspan(g));
    memcpy(g2->grid, g->grid, sizeof(int) * grid_xspan(g) * grid_yspan(g));
    return g2;
}

void resize_grid(grid_t *g, int xmin, int xmax, int ymin, int ymax) {
    grid_t *new;
    int *tmp;

    if (xmin >= g->xmin && xmax <= g->xmax
        && ymin >= g->ymin && ymax <= g->ymax
    )
        return;
    if (xmin > g->xmin) xmin = g->xmin;
    if (xmax < g->xmax) xmax = g->xmax;
    if (ymin > g->ymin) ymin = g->ymin;
    if (ymax < g->ymax) ymax = g->ymax;
    new = new_grid(xmin, xmax, ymin, ymax);
    _copy_grid(new, g);
    tmp = g->grid;
    memcpy(g, new, sizeof(grid_t));
    new->grid = tmp;
    free_grid(new);
}
