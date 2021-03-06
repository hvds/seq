#include <stdlib.h>
#include <string.h>

#include "sym.h"
#include "group.h"

int bitcount[256];

/* For i in 0 .. 8, the set of distinct values in (0 .. 255) with i bits set.
 * These will be used as packed representations of the 8 squares neighbouring
 * a central square.
 */
pack_set_t pack_set[9];

/* What each packed representation maps to under each of the symmetries.
 */
int lookup[(MAXSYM + 1) * 256];

void init_sym(void) {
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < (1 << i); ++j)
            bitcount[j + (1 << i)] = bitcount[j] + 1;

    for (int i = 0; i < 256; ++i) {
        int k = bitcount[i];
        pack_set[k].set[ pack_set[k].count++ ] = i;
    }

    int transform[64] = {
        0, 1, 2, 3, 4, 5, 6, 7,     /* xy */
        2, 1, 0, 4, 3, 7, 6, 5,     /* xY */
        5, 6, 7, 3, 4, 0, 1, 2,     /* Xy */
        7, 6, 5, 4, 3, 2, 1, 0,     /* XY */
        0, 3, 5, 1, 6, 2, 4, 7,     /* yx */
        2, 4, 7, 1, 6, 0, 3, 5,     /* yX */
        5, 3, 0, 6, 1, 7, 4, 2,     /* Yx */
        7, 4, 2, 6, 1, 5, 3, 0      /* YX */
    };
    for (int i = 0; i < 256; ++i) {
        int t[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        for (int j = 0; j < 8; ++j)
            if (i & (1 << j))
                for (int k = 0; k < 8; ++k)
                    t[k] |= 1 << transform[ k * 8 + j ];
        for (int k = 0; k < 8; ++k)
            lookup[k * 256 + i] = t[k];
    }
}

void finish_sym(void) {
    return;
}

/*
    Return true if the dimensions of a grid are transposed by this symmetry.
*/
bool is_transpose(sym_t s) {
    return (s & 4) ? 1 : 0;
}

/*
    Return true if this symmetry is a reflection.
*/
bool is_reflect(sym_t s) {
    return (s == xY || s == Xy || s == yx || s == YX);
}

/*
    Return the symmetry out of the given (packed) symmetries for which
    the location lies on a line of reflection, or 0 if there is none.
*/
sym_t sym_reflect(int syms, int x, int y, loc_t l) {
    for (sym_t s = 1; s <= MAXSYM; ++s)
        if (is_reflect(s)
            && (syms & (1 << s))
            && sym_checkloc(s, x, y, l)
        )
            return s;
    return 0;
}
        
int *sym_lookup(sym_t s) {
    return &lookup[256 * s];
}

/*
    Return true if t is non-canonical, when copmosed with this symmetry.

    If s is a reflection, and s x t = t', this will return false for only
    one of t, t'; if not, it will always return false.
*/
bool sym_dup(sym_t s, sym_t t) {
    switch (s) {
        case xy:
        case XY:
        case yX:
        case Yx:
            return 0;
        case xY:
            return !(t == xy || t == Xy || t == yx || t == YX);
            /* opposing:  xY         XY         yX         Yx */
        case Xy:
            return !(t == xy || t == xY || t == yx || t == YX);
            /*            Xy         XY         Yx         yX */
        case yx:
            return !(t == xy || t == xY || t == Xy || t == XY);
            /*            yx         yX         Yx         YX */
        case YX:
            return !(t == xy || t == xY || t == Xy || t == yx);
            /*            YX         yX         Yx         XY */
    }
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
            return l.x + l.y == (x + y - 2) << 1;
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
