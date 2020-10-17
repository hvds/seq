#include <stdlib.h>

#include "sym.h"
#include "group.h"

void init_sym(void) {
    return;
}

/*
    Return true if the dimensions of a grid are transposed by this symmetry.
*/
bool is_transpose(sym_t s) {
    return (s & 4) ? 1 : 0;
}

/*
    Return true if the specified grid is invariant under this symmetry.
*/
bool sym_check(sym_t s, int x, int y, int *v) {
    int xm = x - 1, ym = y - 1;

    if (x != y && is_transpose(s))
        return 0;

    switch (s) {
        case xy:
            return 1;
        case xY:
            for (int i = 0; i <= xm; ++i)
                for (int j = 0; j + j < ym; ++j)
                    if (v[i * y + j] != v[i * y + ym - j])
                        return 0;
            return 1;
        case Xy:
            for (int i = 0; i + i < xm; ++i)
                for (int j = 0; j <= ym; ++j)
                    if (v[i * y + j] != v[(xm - i) * y + j]) 
                        return 0;
            return 1;
        case XY:
            for (int i = 0; i <= xm; ++i)
                for (int j = 0; j <= ym; ++j)
                    if (v[i * y + j] != v[(xm - i) * y + ym - j])
                        return 0;
            return 1;
        case yx:
            for (int i = 0; i <= xm; ++i)
                for (int j = i; j <= ym; ++j)
                    if (v[i * y + j] != v[j * y + i])
                        return 0;
            return 1;
        case yX:
            for (int i = 0; i <= xm; ++i)
                for (int j = 0; j <= ym; ++j) 
                    if (v[i * y + j] != v[j * y + ym - i])
                        return 0;
            return 1;
        case Yx:
            for (int i = 0; i <= xm; ++i)
                for (int j = 0; j <= ym; ++j) 
                    if (v[i * y + j] != v[(xm - j) * y + i])
                        return 0;
            return 1;
        case YX:            for (int i = 0; i <= xm; ++i)
                for (int j = 0; j <= ym; ++j) 
                    if (v[i * y + j] != v[(xm - j) * y + ym - i])
                        return 0;
            return 1;
    }
}

/*
    Apply this symmetry to the supplied grid.

    Returns a newly malloced array of ints with the transformed grid.
    The dimensions of the grid will be swapped if is_transpose(s).
*/
int *sym_transform(sym_t s, int x, int y, int *vals) {
    int *v = malloc(sizeof(int *) * x * y);
    int xm = x - 1, ym = y - 1;

    switch (s) {
        case xy:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[i * y + j] = vals[i * y + j];
            break;
        case xY:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j)
                    v[i * y + ym - j] = vals[i * y + j];
            break;
        case Xy:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[(xm - i) * y + j] = vals[i * y + j];
            break;
        case XY:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[(xm - i) * y + ym - j] = vals[i * y + j];
            break;
        case yx:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[j * x + i] = vals[i * y + j];
            break;
        case yX:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j)
                    v[j * x + xm - i] = vals[i * y + j];
            break;
        case Yx:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[(ym - j) * x + i] = vals[i * y + j];
            break;
        case YX:
            for (int i = 0; i < x; ++i)
                for (int j = 0; j < y; ++j) 
                    v[(ym - j) * x + xm - i] = vals[i * y + j];
            break;
    }
    return v;
}

/*
    Return true if the specified location is invariant under this symmetry.
*/
bool sym_checkloc(sym_t s, int x, int y, loc_t l) {
    switch (s) {
        case xy:
            return 1;
        case xY:
            return (l.y << 1) == y - 1;
        case Xy:
            return (l.x << 1) == x - 1;
        case XY:
            return ((l.x << 1) == x - 1) && ((l.y << 1) == y - 1);
        case yx:
            return l.x == l.y;
        case yX:
            return ((l.x << 1) == x - 1) && ((l.y << 1) == y - 1);
        case Yx:
            return ((l.x << 1) == x - 1) && ((l.y << 1) == y - 1);
        case YX:
            return ((l.x + l.y) >> 1) == x + y - 2;
    }
}

/*
    Apply this symmetry to the supplied location.
*/
loc_t sym_transloc(sym_t s, int x, int y, loc_t l) {
    switch (s) {
        case xy:
            return (loc_t){ l.x, l.y };
        case xY:
            return (loc_t){ l.x, y - 1 - l.y };
        case Xy:
            return (loc_t){ x - 1 - l.x, l.y };
        case XY:
            return (loc_t){ x - 1 - l.x, y - 1 - l.y };
        case yx:
            return (loc_t){ l.y, l.x };
        case yX:
            return (loc_t){ l.y, x - 1 - l.x };
        case Yx:
            return (loc_t){ y - 1 - l.y, l.x };
        case YX:
            return (loc_t){ y - 1 - l.y, x - 1 - l.x };
    }
}
