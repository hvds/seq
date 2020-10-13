#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "group.h"
#include "sym.h"

grouplist_t *cache_seed[9];
typedef struct seed_s {
    int count;
    int bits[13];  /* max needed */
} seed_t;
seed_t seed_base[9] = {
    { 0, {} },
    { 0, {} },
    { 6, { 0300, 0210, 0201, 0050, 0500, 0401 } },
    { 10, { 0700, 0610, 0601, 0602, 0604, 0640, 0504, 0502, 0412, 0250 } },
    { 13, { 0505, 0704, 0514, 0710, 0702, 0550, 0512, 0641, 0611, 0603,
            0642, 0342, 0252 } },
    { 10, { 0057, 0147, 0156, 0155, 0153, 0117, 0253, 0255, 0345, 0307 } },
    { 6, { 0457, 0547, 0556, 0707, 0257, 0356 } },
    { 2, { 0776, 0775 } },
    { 1, { 0777 } }
};

void init_group(void) {
    return;
}

group_t *new_group(int x, int y, int sym, int* vals) {
    group_t *g = malloc(sizeof(group_t));
    int maxsum = 0;

    g->x = x;
    g->y = y;
    g->sym = sym;
    g->vals = vals;
    /* caller will increment; freed on decrement to zero */
    g->refcount = 0;

    /* initialise to AVAIL = 0 */
    g->avail = calloc((x + 2) * (y + 2), sizeof(avail_t));
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            if (g->vals[i * y + j])
                g->avail[(i + 1) * (y + 2) + (j + 1)] = USED;

    g->sum_chains = calloc((x + 2) * (y + 2), sizeof(int));
    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j) {
            int v = g->vals[i * y + j];
            if (v == 0)
                continue;
            for (int di = 0; di < 3; ++di)
                for (int dj = 0; dj < 3; ++dj) {
                    int off = (i + di) * (y + 2) + (j + dj);
                    if (v > 1 && g->avail[off] == AVAIL)
                        g->avail[off] = RES;
                    if (g->avail[off] != USED) {
                        g->sum_chains[off] += v;
                        if (g->sum_chains[off] > maxsum)
                            maxsum = g->sum_chains[off];
                    }
                }
        }

    if (sym)
        for (int s = 1; s <= MAXSYM; ++s) {
            if (!(sym & (1 << s)))
                continue;
            for (int i = -1; i < x + 1; ++i)
                for (int j = -1; j < y + 1; ++j) {
                    loc_t l = sym_transloc(s, g, (loc_t){ i, j });
                    if (l.x * (y + 2) + l.y < i * (y + 2) + j)
                        g->sum_chains[(i + 1) * (y + 2) + (j + 1)] = 0;
                }
        }

    g->maxsum = maxsum;
    g->sum_heads = malloc((maxsum + 1) * sizeof(int));
    memset(g->sum_heads, 0xff, (maxsum + 1) * sizeof(int));

    /* Walk backwards turning the actual sums into linked lists headed by
     * sum_heads[sum]. We mask out USED and zero locations.
     * FIXME: should also mask out symmetries
     */
    for (int i = (x + 2) * (y + 2) - 1; i >= 0; --i) {
        int sum = g->sum_chains[i];
        if (sum == 0 || g->avail[i] == USED) {
            g->sum_chains[i] = -1;
        } else {
            g->sum_chains[i] = g->sum_heads[sum];
            g->sum_heads[sum] = i;
        }
    }

    return g;
}

void ref_group(group_t *g) {
    ++g->refcount;
}

void unref_group(group_t *g) {
    if (--g->refcount == 0) {
        free(g->vals);
        free(g->avail);
        free(g->sum_heads);
        free(g->sum_chains);
        free(g);
    }
}

grouplist_t *new_grouplist(int size) {
    grouplist_t *gl = malloc(sizeof(grouplist_t) + size * sizeof(group_t *));
    gl->count = size;
    return gl;
}

void free_grouplist(grouplist_t *gl) {
    for (int i = 0; i < gl->count; ++i) {
        unref_group(gl->g[i]);
    }
    free(gl);
}

void print_list(char *s, int x, int y, void *p, int size, bool nl) {
    printf("%s", s);
    for (int i = 0; i < x; ++i) {
        if (i)
            printf("; ");
        for (int j = 0; j < y; ++j) {
            if (j)
                printf(" ");
            printf("%d", (
                size == 4 ? ((int *)p)[i * y + j]
                : size == 2 ? ((short *)p)[i * y + j]
                : size == 1 ? ((char *)p)[i * y + j]
                : 0
            ));
        }
    }
    if (nl)
        printf("\n");
}

void print_group(group_t *g) {
    print_list("", g->x, g->y, g->vals, sizeof(int), 0);
}

void dprint_group(group_t *g) {
    printf("(%p) x=%d y=%d sym=%d maxsum=%d refcount=%d\n",
            g, g->x, g->y, g->sym, g->maxsum, g->refcount);
    print_list(" vals: ", g->x, g->y, g->vals, sizeof(int), 1);
    print_list(" avail: ", g->x + 2, g->y + 2, g->avail, sizeof(avail_t), 1);
    print_list(" heads: ", 1, g->maxsum + 1, g->sum_heads, sizeof(int), 1);
    print_list(" chain: ", g->x + 2, g->y + 2, g->sum_chains, sizeof(int), 1);
}

group_t *group_seedbits(int k, int bits) {
    int *vals = malloc(sizeof(int) * 9);
    int x = 3, y = 3, sym = 0, x0 = 0, y0 = 0;
    group_t *g;

    /* shrink to fit */
    if ((bits & 0b000000111) == 0)
        --x;
    if ((bits & 0b111000000) == 0)
        --x, ++x0;
    if ((bits & 0b001001001) == 0)
        --y;
    if ((bits & 0b100100100) == 0)
        --y, ++y0;

    for (int i = 0; i < x; ++i)
        for (int j = 0; j < y; ++j)
            vals[i * y + j] = (bits & (1 << (8 - ((i + x0) * 3 + j + y0))))
                    ? 1 : 0;
    vals[(1 - x0) * y + (1 - y0)] = k;

    for (int s = 1; s <= MAXSYM; ++s)
        if (sym_check(s, x, y, vals))
            sym |= (1 << s);

    g = new_group(x, y, sym, vals);
    ref_group(g);
    return g;
}

grouplist_t *group_seed(int k) {
    if (k < 2 || k > 8) {
        fprintf(stderr, "Error: group_seed(%d) called\n", k);
        exit(1);
    }
    if (!cache_seed[k]) {
        seed_t *seed = &seed_base[k];
        grouplist_t *gl = new_grouplist(seed->count);
        for (int i = 0; i < seed->count; ++i) {
            gl->g[i] = group_seedbits(k, seed->bits[i]);
            ref_group(gl->g[i]);
        }
        cache_seed[k] = gl;
    }
    return cache_seed[k];
}

sym_t next_sym(group_t *g, loc_t l) {
    int osym = g->sym, nsym = 0;

    if (osym)
        for (sym_t s = 1; s <= MAXSYM; ++s)
            if (osym & (1 << s)) {
                loc_t tl = sym_transloc(s, g, l);
                if (tl.x == l.x && tl.y == l.y)
                    nsym |= (1 << s);
            }
    return nsym;
}

group_t *group_place(group_t *g, loc_t loc, int k) {
    int x = g->x, y = g->y, x0 = 0, y0 = 0;
    int *vals, sym = next_sym(g, loc);

    if (loc.x < 0)
        ++x, ++x0;
    if (loc.x >= x)
        ++x;
    if (loc.y < 0)
        ++y, ++y0;
    if (loc.y >= y)
        ++y;

    vals = calloc(x * y, sizeof(int));
    for (int i = 0; i < g->x; ++i)
        for (int j = 0; j < g->y; ++j)
            vals[(i + x0) * y + (j + y0)] = g->vals[i * g->y + j];
    vals[(loc.x + x0) * y + (loc.y + y0)] = k;

    return new_group(x, y, sym, vals);
}

int _comb(int n, int d) {
    int x = 1;
    for (int i = 0; i < d; ++i)
        x = x * (n - i) / (i + 1);
    return x;
}

grouplist_t *group_place_with(group_t *g, loc_t loc, int k, int use) {
    loc_t avail[9];
    int availc = 0;
    int stack[9];
    int sp = 0;
    int count;
    grouplist_t *result;
    int ri = 0;
    int x, y, x0, y0, xmin, ymin, xmax, ymax;
    int *vals;

    /* put the centre location into the list to simplify bounding box checks */
    avail[availc++] = loc;

    /* find available places surrounding the location to put extra 1s */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == 1 && j == 1)
                continue;
            if (loc.x + i < 0 || loc.x +i >= g->x + 2
                || loc.y + j < 0 || loc.y + j >= g->y + 2
                || g->avail[(loc.x + i) * (g->y + 2) + (loc.y + j)] == AVAIL
            )
                avail[availc++] = (loc_t){ loc.x + i - 1, loc.y + j - 1 };
        }

    if (availc - 1 < use)
        count = 0;
    else
        count = _comb(availc - 1, use);
    result = new_grouplist(count);

    stack[sp++] = 0;    /* place central loc for bounding box check */
    stack[sp] = 0;      /* skip 0 otherwise */
    while (sp > 0) {    /* stop when we get back to central loc at stack[0] */
        ++stack[sp];
        if (stack[sp] >= availc) {
            --sp;
            continue;
        }
        if (sp < use) {
            stack[sp + 1] = stack[sp];
            ++sp;
            continue;
        }
        x = g->x;
        y = g->y;
        x0 = y0 = xmin = ymin = xmax = ymax = 0;
        for (int i = 0; i <= sp; ++i) {
            if (avail[stack[i]].x < xmin)
                xmin = avail[stack[i]].x;
            if (avail[stack[i]].x > xmax)
                xmax = avail[stack[i]].x;
            if (avail[stack[i]].y < ymin)
                ymin = avail[stack[i]].y;
            if (avail[stack[i]].y > ymax)
                ymax = avail[stack[i]].y;
        }
        if (xmin < 0)
            x -= xmin, x0 -= xmin;
        if (xmax >= g->x)
            x += xmax + 1 - g->x;
        if (ymin < 0)
            y -= ymin, y0 -= ymin;
        if (ymax >= g->y)
            y += ymax + 1 - g->y;
        vals = calloc(x * y, sizeof(int));
        for (int i = 0; i < g->x; ++i)
            for (int j = 0; j < g->y; ++j)
                vals[(i + x0) * y + (j + y0)] = g->vals[i * g->y + j];
        vals[(loc.x + x0) * y + (loc.y + y0)] = k;
        for (int i = 0; i < use; ++i) {
            loc_t loc1 = avail[stack[i + 1]];
            vals[(loc1.x + x0) * y + (loc1.y + y0)] = 1;
        }
        result->g[ri] = new_group(x, y, 0, vals);
        ref_group(result->g[ri++]);
    }
if (ri != count) {
    printf("found %d of %d (for comb(%d, %d))\n", ri, count, availc - 1, use);
    result->count = ri;
}
    return result;
}

grouplist_t *coalesce_group(
    group_t *ga, loc_t la, group_t *gb, loc_t lb, int k, int use
) {
    int ax = ga->x, ay = ga->y;
    int bx = gb->x, by = gb->y;
    int afree = 0, aneed = 0, bfree = 0, bneed = 0;
    /* not avail_t [], since we're passing to sym_transform() */
    int aavail[9], bavail[9];
    int maxavail, maxcombs;
    grouplist_t *result;
    int ri = 0;
    int axmin = -la.x, axmax = ax - la.x, aymin = -la.y, aymax = ay - la.y;

    /* first look at the 3x3 square around the common location, to see
     * if the two groups can in principle fit together around there
     * _and_ leave room for an additional 'use' 1s/
     */
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
            if (i == 1 && j == 1)
                continue;
            if (la.x + i < 0 || la.x + i >= ax + 2
                || la.y + j < 0 || la.y + j >= ay + 2
            ) {
                ++afree;
                aavail[i * 3 + j] = AVAIL;
            } else {
                avail_t a = ga->avail[(la.x + i) * (ay + 2) + (la.y + j)];
                if (a == AVAIL)
                    ++afree;
                else if (a == USED)
                    ++aneed;
                aavail[i * 3 + j] = a;
            }
            if (lb.x + i < 0 || lb.x + i >= bx + 2
                || lb.y + j < 0 || lb.y + j >= by + 2
            ) {
                ++bfree;
                bavail[i * 3 + j] = AVAIL;
            } else {
                avail_t a = gb->avail[(lb.x + i) * (by + 2) + (lb.y + j)];
                if (a == AVAIL)
                    ++bfree;
                else if (a == USED)
                    ++bneed;
                bavail[i * 3 + j] = a;
            }
        }

    if (afree < bneed + use || bfree < aneed + use)
        return new_grouplist(0);

    /* prepare a space big enough for the maximum possible number of results */
    maxavail = (afree - bneed < bfree - aneed)
        ? afree - bneed
        : bfree - aneed;
    maxcombs = _comb(maxavail, use);
    result = new_grouplist((MAXSYM + 1) * maxcombs);

    /* Now for each transform of gb, check first whether its 3x3 square can
     * be placed over ga's without conflict, and still leaving enough free
     * spots for us to add 'use' 1s. If it can, try the full monty.
     */
    for (int s = 0; s <= MAXSYM; ++s) {
        int *tavail = sym_transform(s, 3, 3, bavail);
        int free_c = 0;
        loc_t lfree[9];
        int ok = 1;

        lfree[free_c++] = la;

        for (int i = 0; ok && i <= 2; ++i)
            for (int j = 0; ok && j <= 2; ++j) {
                avail_t aa = aavail[i * 3 + j];
                avail_t ta = tavail[i * 3 + j];

                if (i == 1 && j == 1)
                    continue;

                if ((aa == USED && ta != AVAIL)
                    || (ta == USED && aa != AVAIL)
                ) {
                    ok = 0;
                    break;
                }
                if (aa == AVAIL && ta == AVAIL)
                    lfree[free_c++] = (loc_t){ la.x + i - 1, la.y + j - 1 };
            }
        if (ok && free_c - 1 < use)
            ok = 0;
        free(tavail);

        if (ok) {
            int *tvals = sym_transform(s, bx, by, gb->vals);
            loc_t lt = sym_transloc(s, gb, lb);
            bool trans = is_transpose(s);
            int tx = trans ? by : bx, ty = trans ? bx : by;
            int txmin = -lt.x, txmax = tx - lt.x;
            int tymin = -lt.y, tymax = ty - lt.y;
            int xmin = (axmin < txmin) ? axmin : txmin;
            int xmax = (axmax > txmax) ? axmax : txmax;
            int ymin = (aymin < tymin) ? aymin : tymin;
            int ymax = (aymax > tymax) ? aymax : tymax;
            int dax = axmin - xmin, day = aymin - ymin;
            int dtx = txmin - xmin, dty = tymin - ymin;
            int cx = xmax - xmin, cy = ymax - ymin;
            loc_t lc = (loc_t){ la.x + dax, la.y + day };
            int *cvals = calloc(cx * cy, sizeof(int));

            for (int i = 0; ok && i < cx; ++i)
                for (int j = 0; ok && j < cy; ++j) {
                    bool ba = (i >= dax && i <= dax + ax - 1
                        && j >= day && j <= day + ay - 1);
                    bool bt = (i >= dtx && i <= dtx + tx - 1
                        && j >= dty && j <= dty + ty - 1);
                    int va = ba ? ga->vals[(i - dax) * ay + (j - day)] : 0;
                    int vt = bt ? tvals[(i - dtx) * ty + (j - dty)] : 0;
                    if (va && vt) {
                        ok = 0;
                        break;
                    }
                    cvals[i * cy + j] = va ? va : vt;
                }

            int stack[9];
            int sp = 0;
            stack[sp++] = 0;    /* place central loc for bounding box check */
            stack[sp] = 0;      /* skip 0 otherwise */
            while (ok && sp > 0) { /* stop when we get back to stack[0] */
                ++stack[sp];
                if (stack[sp] >= free_c) {
                    --sp;
                    continue;
                }
                if (sp < use) {
                    stack[sp + 1] = stack[sp];
                    ++sp;
                    continue;
                }
                if (use == 0)
                    sp = 0;

                int fx = cx, fy = cy;
                int fx0 = 0, fy0 = 0;
                int fxmin = 0, fymin = 0, fxmax = 0, fymax = 0;
                for (int i = 0; i <= sp; ++i) {
                    if (lfree[stack[i]].x < fxmin)
                        fxmin = lfree[stack[i]].x;
                    if (lfree[stack[i]].x > fxmax)
                        fxmax = lfree[stack[i]].x;
                    if (lfree[stack[i]].y < fymin)
                        fymin = lfree[stack[i]].y;
                    if (lfree[stack[i]].y > fymax)
                        fymax = lfree[stack[i]].y;
                }
                if (fxmin < 0)
                    fx -= fxmin, fx0 -= fxmin;
                if (fxmax >= cx)
                    fx += fxmax + 1 - cx;
                if (fymin < 0)
                    fy -= fymin, fy0 -= fymin;
                if (fymax >= cy)
                    fy += fymax + 1 - cy;

                int *fvals = calloc(fx * fy, sizeof(int));
                for (int i = 0; i < cx; ++i)
                    for (int j = 0; j < cy; ++j)
                        fvals[(i + fx0) * fy + (j + fy0)] = cvals[i * cy + j];
                fvals[(la.x + fx0) * fy + (la.y + fy0)] = k;
                for (int i = 0; i < use; ++i) {
                    loc_t loc1 = lfree[stack[i + 1]];
                    fvals[(loc1.x + fx0) * fy + (loc1.y + fy0)] = 1;
                }
                result->g[ri] = new_group(fx, fy, 0, fvals);
                ref_group(result->g[ri++]);
            }
            free(cvals);
            free(tvals);
        }
    }
    result->count = ri;
    return result;
}
