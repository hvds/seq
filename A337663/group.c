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

/*
 * Packed representation for symmetry-reduced ways to arrange n 1s in the
 * 8 squares neighbouring a central point.
 * Bit numbers are [ 765 // 4*3 // 210 ] ie bit 7 represents { -1, -1 }
 * relative to the centre. Seed bits for k are the complement of those
 * for (8 - k).
 */
seed_t seed_base[9] = {
    { 0, {} },  /* actually { 1, { 0 } }, but we never need it */
    { 0, {} },  /* actually { 2, { 1, 2 } }, but we never need it */
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

void finish_group(void) {
    for (int i = 0; i <= 8; ++i)
        if (cache_seed[i])
            free_grouplist(cache_seed[i]);
}

/*
 * Allocate and initialize a new group struct with the supplied details.
 *
 * Initial refcount is zero; caller is expected to give it its initial
 * increment; on decrement to 0 both the group struct and 'vals' will
 * be freed.
 */
group_t *new_group(int x, int y, int sym, int* vals) {
    group_t *g = malloc(sizeof(group_t));
    int maxsum = 0;

    g->x = x;
    g->y = y;
    g->sym = sym;
    g->vals = vals;
    g->tvals = (int **)NULL;
    /* caller will increment; freed on decrement to zero */
    g->refcount = 0;

    /* initialise to AVAIL = 0 */
    g->avail = calloc((x + 2) * (y + 2), sizeof(avail_t));
    g->tavail = (avail_t**)NULL;
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
                    loc_t l = sym_transloc(s, x, y, (loc_t){ i, j });
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

/*
 * Increment the refcount of the group.
 */
void ref_group(group_t *g) {
    ++g->refcount;
}

/*
 * Decrement the refcount of the group; free it if it hits zero.
 */
void unref_group(group_t *g) {
    if (--g->refcount == 0) {
        free(g->vals);
        if (g->tvals) {
            for (int i = 1; i <= MAXSYM; ++i)
                free(g->tvals[i]);
            free(g->tvals);
        }
        free(g->avail);
        if (g->tavail) {
            for (int i = 1; i <= MAXSYM; ++i)
                free(g->tavail[i]);
            free(g->tavail);
        }
        free(g->sum_heads);
        free(g->sum_chains);
        free(g);
    }
}

/*
 * Allocate a grouplist of the specified size; set the size to 'size'.
 */
grouplist_t *new_grouplist(int size) {
    grouplist_t *gl = malloc(sizeof(grouplist_t) + size * sizeof(group_t *));
    gl->count = size;
    memset(&gl->g[0], 0, size * sizeof(group_t *));
    return gl;
}

/*
 * Free the grouplist, and dereference any groups within it.
 */
void free_grouplist(grouplist_t *gl) {
    for (int i = 0; i < gl->count; ++i)
        unref_group(gl->g[i]);
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

/*
 * Return an array of ints representing this group's vals transformed
 * under the specified symmetry 's'.
 *
 * The dimensions of the result will be swapped if is_transpose(s).
 * The array is cached in the group, and will be freed when the group is.
 */
int *trans_vals(group_t *g, sym_t s) {
    if (!g->tvals)
        g->tvals = calloc(MAXSYM + 1, sizeof(int *));
    if (!g->tvals[s])
        g->tvals[s] = sym_transform(s, g->x, g->y, g->vals);
    return g->tvals[s];
}

/*
 * Return an array of avail_t representing this group's avail transformed
 * under the specified symmetry 's'.
 *
 * The dimensions of the result will be swapped if is_transpose(s).
 * The array is cached in the group, and will be freed when the group is.
 */
avail_t *trans_avail(group_t *g, sym_t s) {
    if (!g->tavail)
        g->tavail = calloc(MAXSYM + 1, sizeof(avail_t *));
    if (!g->tavail[s])
        g->tavail[s] = (avail_t *)sym_transform(s, g->x, g->y, (int *)g->avail);
    return g->tavail[s];
}

/*
 * Construct and return a group consisting of a value k surrounded by
 * 0 or more 1s in the positions indicated by 'bits', using the
 * representation described for 'seed_base'.
 *
 * The group has an initial refcount of 1 for its copy in the cache,
 * which is freed via finish_group().
 */
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

/*
 * Return a grouplist with all groups (reduced by symmetry) consisting of
 * a central value k surrounded by k 1s in the neighbouring 8 squares.
 *
 * The grouplist is cached, and freed via finish_group().
 */
grouplist_t *group_seed(int k) {
    if (k < 2 || k > 8) {
        fprintf(stderr, "Error: group_seed(%d) called\n", k);
        exit(1);
    }
    if (!cache_seed[k]) {
        seed_t *seed = &seed_base[k];
        grouplist_t *gl = new_grouplist(seed->count);
        for (int i = 0; i < seed->count; ++i)
            gl->g[i] = group_seedbits(k, seed->bits[i]);
        cache_seed[k] = gl;
    }
    return cache_seed[k];
}

/*
 * Return a packed representation of the symmetries preserved by this
 * group when you add a new unique value at the specified location.
 */
sym_t next_sym(group_t *g, loc_t l) {
    int osym = g->sym, nsym = 0;

    if (osym)
        for (sym_t s = 1; s <= MAXSYM; ++s)
            if (osym & (1 << s))
                if (sym_checkloc(s, g->x, g->y, l))
                    nsym |= (1 << s);
    return nsym;
}

/*
 * Construct and return a new group formed by adding a new value k
 * at the specified location in this group.
 *
 * Equivalent to group_place_with(g, loc, k, 0).
 * The new group has an initial refcount of 0.
 */
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

/* nCd = n! / d! (n - d)! */
int _comb(int n, int d) {
    int x = 1;
    for (int i = 0; i < d; ++i)
        x = x * (n - i) / (i + 1);
    return x;
}

/*
 * Construct and return a grouplist of the new groups formed by adding
 * a new value k at the specified location in this group, along with
 * 'use' 1s in any of the 8 surrounding squares that are available.
 */
grouplist_t *group_place_with(group_t *g, loc_t loc, int k, int use) {
    int maxcount;
    grouplist_t *result;
    int ri = 0;
    int x, y, x0, y0, xmin, ymin, xmax, ymax;
    int packed = 0;
    int *vals, reflect, *refset;
    pack_set_t *maybes;

    /* create a packed representation of the available slots in which
     * to place 1s
     */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == 1 && j == 1)
                continue;
            packed <<= 1;
            if (loc.x + i < 0 || loc.x + i >= g->x + 2
                || loc.y + j < 0 || loc.y + j >= g->y + 2
                || g->avail[(loc.x + i) * (g->y + 2) + (loc.y + j)] == AVAIL
            )
                packed |= 1;
        }

    if (bitcount[packed] < use)
        return new_grouplist(0);

    maxcount = _comb(bitcount[packed], use);
    result = new_grouplist(maxcount);

    reflect = sym_reflect(g->sym, g->x, g->y, loc);
    refset = sym_lookup(reflect);
    maybes = &pack_set[use];

    for (int i = 0; i < maybes->count; ++i) {
        int maybe = maybes->set[i];
        /* skip if it requires spots we don't have */
        if (maybe & ~packed)
            continue;
        /* skip if it's a reflection of one we've already taken */
        if (refset[maybe] < maybe)
            continue;

        x = g->x;
        y = g->y;
        x0 = 0;
        y0 = 0;
        xmin = (maybe & 0b11100000) ? loc.x - 1 : loc.x;
        xmax = (maybe & 0b00000111) ? loc.x + 1 : loc.x;
        ymin = (maybe & 0b10010100) ? loc.y - 1 : loc.y;
        ymax = (maybe & 0b00101001) ? loc.y + 1 : loc.y;
        if (xmin < 0)
            x -= xmin, x0 -= xmin;
        if (xmax >= g->x)
            x += xmax + 1 - g->x;
        if (ymin < 0)
            y -= ymin, y0 -= ymin;
        if (ymax >= g->y)
            y += ymax + 1 - g->y;

        int *vals = calloc(x * y, sizeof(int));
        for (int i = 0; i < g->x; ++i)
            for (int j = 0; j < g->y; ++j)
                vals[(i + x0) * y + (j + y0)] = g->vals[i * g->y + j];
        vals[(loc.x + x0) * y + (loc.y + y0)] = k;
        if (maybe & 0b10000000)
            vals[(loc.x + x0 - 1) * y + (loc.y + y0 - 1)] = 1;
        if (maybe & 0b01000000)
            vals[(loc.x + x0 - 1) * y + (loc.y + y0    )] = 1;
        if (maybe & 0b00100000)
            vals[(loc.x + x0 - 1) * y + (loc.y + y0 + 1)] = 1;
        if (maybe & 0b00010000)
            vals[(loc.x + x0    ) * y + (loc.y + y0 - 1)] = 1;
        if (maybe & 0b00001000)
            vals[(loc.x + x0    ) * y + (loc.y + y0 + 1)] = 1;
        if (maybe & 0b00000100)
            vals[(loc.x + x0 + 1) * y + (loc.y + y0 - 1)] = 1;
        if (maybe & 0b00000010)
            vals[(loc.x + x0 + 1) * y + (loc.y + y0    )] = 1;
        if (maybe & 0b00000001)
            vals[(loc.x + x0 + 1) * y + (loc.y + y0 + 1)] = 1;

        result->g[ri] = new_group(x, y, 0, vals);
        ref_group(result->g[ri++]);
    }
    if (ri < maxcount)
        result->count = ri;
    return result;
}

/*
 * Construct and return a grouplist of the new groups formed by adding
 * a new value k at a common point connecting two existing groups, such
 * that the location of that point is as specified for each of the two
 * groups, along with 'use' 1s in any of the 8 surrounding squares that
 * are available.
 *
 * All relevant symmetries of the two groups are considered, and the
 * results include only those combinations that allowed the two groups
 * to be combined according to location availability checks.
 */
grouplist_t *coalesce_group(
    group_t *ga, loc_t la, group_t *gb, loc_t lb, int k, int use
) {
    int ax = ga->x, ay = ga->y;
    int bx = gb->x, by = gb->y;
    int afree = 0, aneed = 0, bfree = 0, bneed = 0;
    int maxavail, maxcombs;
    grouplist_t *result;
    int ri = 0;
    int axmin = -la.x, axmax = ax - la.x, aymin = -la.y, aymax = ay - la.y;
    sym_t reflect;

    /* first look at the 3x3 square around the common location, to see
     * if the two groups can in principle fit together around there
     * _and_ leave room for an additional 'use' 1s/
     */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == 1 && j == 1)
                continue;
            afree <<= 1;
            aneed <<= 1;
            bfree <<= 1;
            bneed <<= 1;
            if (la.x + i < 0 || la.x + i >= ax + 2
                || la.y + j < 0 || la.y + j >= ay + 2
            ) {
                afree |= 1;
            } else {
                avail_t a = ga->avail[(la.x + i) * (ay + 2) + (la.y + j)];
                if (a == AVAIL)
                    afree |= 1;
                else if (a == USED)
                    aneed |= 1;
            }
            if (lb.x + i < 0 || lb.x + i >= bx + 2
                || lb.y + j < 0 || lb.y + j >= by + 2
            ) {
                bfree |= 1;
            } else {
                avail_t a = gb->avail[(lb.x + i) * (by + 2) + (lb.y + j)];
                if (a == AVAIL)
                    bfree |= 1;
                else if (a == USED)
                    bneed |= 1;
            }
        }

    if (bitcount[afree] < bitcount[bneed] + use
        || bitcount[bfree] < bitcount[aneed] + use
    )
        return new_grouplist(0);

    /* prepare a space big enough for the maximum possible number of results */
    maxavail = (
        bitcount[afree] - bitcount[bneed] < bitcount[bfree] - bitcount[aneed]
    )
        ? bitcount[afree] - bitcount[bneed]
        : bitcount[bfree] - bitcount[aneed];
    maxcombs = _comb(maxavail, use);
    result = new_grouplist((MAXSYM + 1) * maxcombs);

    /* check if the first group's location is on a line of reflection */
    reflect = sym_reflect(ga->sym, ax, ay, la);

    /* Now for each transform of gb, check first whether its 3x3 square can
     * be placed over ga's without conflict, and still leaving enough free
     * spots for us to add 'use' 1s. If it can, try the full monty.
     */
    for (int s = 0; s <= MAXSYM; ++s) {
        /* if this transform is a symmetry of the group _and_ the location
         * is invariant under the symmetry, we've already checked it.
         */
        if ((gb->sym & (1 << s)) && sym_checkloc(s, bx, by, lb))
            continue;

        /* if the shared location is on a line of symmetry of the other group,
         * check this symmetry on this group only if it is canonical under
         * that symmetry on the other group.
         */
        if (reflect && sym_dup(reflect, s))
            continue;

        /* if the shared location is on the same line of reflection of
         * both groups, could deduplicate asymmetric placement of 1s using
         * the same 'if (refset[maybe] < maybe)' logic as group_place_with.
         * TODO - need something like:
         *    breflect = sym_reflect(gb->sym, bx, by, lb)
         *    treflect = sym_combine(breflect, s)
         *    refset = sym_lookup(areflect == treflect ? areflect : 0)
         */

        int *slookup = sym_lookup(s);
        int btfree = slookup[bfree];
        int btneed = slookup[bneed];

        /* must be able to superimpose without conflict */
        if ((aneed & ~btfree) || (btneed & ~afree))
            continue;

        int tfree = afree & btfree & ~aneed & ~btneed;
        /* superimposed, must have 'use' spare spots available */
        if (bitcount[tfree] < use)
            continue;

        int *tvals = trans_vals(gb, s);
        loc_t lt = sym_transloc(s, bx, by, lb);
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
        int ok = 1;

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

        if (ok) {
            pack_set_t *maybes = &pack_set[use];
            for (int i = 0; i < maybes->count; ++i) {
                int maybe = maybes->set[i];
                /* skip if it requires spots we don't have */
                if (maybe & ~tfree)
                    continue;

                int fx = cx, fy = cy;
                int fx0 = 0, fy0 = 0;
                int fxmin = (maybe & 0b11100000) ? lc.x - 1 : lc.x;
                int fxmax = (maybe & 0b00000111) ? lc.x + 1 : lc.x;
                int fymin = (maybe & 0b10010100) ? lc.y - 1 : lc.y;
                int fymax = (maybe & 0b00101001) ? lc.y + 1 : lc.y;
                if (fxmin < 0)
                    fx -= fxmin, fx0 -= fxmin;
                if (fxmax >= cx)
                    fx += fxmax + 1 - cx;
                if (fymin < 0)
                    fy -= fymin, fy0 -= fymin;
                if (fymax >= cy)
                    fy += fymax + 1 - cy;
                loc_t lf = (loc_t){ lc.x + fx0, lc.y + fy0 };

                int *fvals = calloc(fx * fy, sizeof(int));
                for (int i = 0; i < cx; ++i)
                    for (int j = 0; j < cy; ++j)
                        fvals[(i + fx0) * fy + (j + fy0)] = cvals[i * cy + j];
                fvals[lf.x * fy + lf.y] = k;
                if (maybe & 0b10000000)
                    fvals[(lf.x - 1) * fy + (lf.y - 1)] = 1;
                if (maybe & 0b01000000)
                    fvals[(lf.x - 1) * fy + (lf.y    )] = 1;
                if (maybe & 0b00100000)
                    fvals[(lf.x - 1) * fy + (lf.y + 1)] = 1;
                if (maybe & 0b00010000)
                    fvals[(lf.x    ) * fy + (lf.y - 1)] = 1;
                if (maybe & 0b00001000)
                    fvals[(lf.x    ) * fy + (lf.y + 1)] = 1;
                if (maybe & 0b00000100)
                    fvals[(lf.x + 1) * fy + (lf.y - 1)] = 1;
                if (maybe & 0b00000010)
                    fvals[(lf.x + 1) * fy + (lf.y    )] = 1;
                if (maybe & 0b00000001)
                    fvals[(lf.x + 1) * fy + (lf.y + 1)] = 1;

                result->g[ri] = new_group(fx, fy, 0, fvals);
                ref_group(result->g[ri++]);
            }
        }
        free(cvals);
    }
    result->count = ri;
    return result;
}
