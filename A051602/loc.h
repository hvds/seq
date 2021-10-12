#ifndef LOC_H
#define LOC_H

#ifdef DEBUG
#   include <assert.h>
#else
#   define assert(x) /* noop */
#endif

typedef struct {
    int x;
    int y;
} loc_t;

typedef struct {
    loc_t p;
    int power;
} loc2p_t;

typedef struct {
    loc_t min;
    loc_t max;
} span_t;

typedef struct {
    loc_t p[2];
} pair_t;

/* We use two different loclists: the global point[], which is initialized
 * to be big enough for all points, and the per-context seen[], which
 * grows as needed.
 * For point[], knowledge of the recursion level is used to determine
 * which entries are relevant; this is mostly captured in the variable
 * 'points', so the *_lim() functions below are used to interrogate it.
 * For seen[], new points are appended, and seen->used tracks how many
 * entries it currently has. The bare (not *_lim()) functions are used
 * to interrogate them.
 */
typedef struct {
    int size;
    int used;
    loc_t *list;
} loclist_t;

typedef struct {
    int size;
    int used;
    loc2p_t *list;
} loc2plist_t;

typedef struct {
    int size;
    int used;
    pair_t *list;
} pairlist_t;

loclist_t *new_loclist(int size);
void free_loclist(loclist_t *ll);
void resize_loclist(loclist_t *ll, int size);
loclist_t *dup_loclist(loclist_t *ll);

loc2plist_t *new_loc2plist(int size);
void free_loc2plist(loc2plist_t *l2l);

pairlist_t *new_pairlist(int size);
void free_pairlist(pairlist_t *pl);
void resize_pairlist(pairlist_t *pl, int size);
pairlist_t *dup_pairlist(pairlist_t *pl);

static loc_t loc_diff(loc_t s1, loc_t s2) {
    loc_t d;
    d.x = s2.x - s1.x;
    d.y = s2.y - s1.y;
    return d;
}

static loc_t loc_sum(loc_t s1, loc_t s2) {
    loc_t d;
    d.x = s1.x + s2.x;
    d.y = s1.y + s2.y;
    return d;
}

static loc_t loc_rot90(loc_t src, loc_t diff) {
    loc_t d;
    d.x = src.x - diff.y;
    d.y = src.y + diff.x;
    return d;
}

static loc_t loc_rot270(loc_t src, loc_t diff) {
    loc_t d;
    d.x = src.x + diff.y;
    d.y = src.y - diff.x;
    return d;
}

static int loc_eq(loc_t s1, loc_t s2) {
    return s1.x == s2.x && s1.y == s2.y;
}

/* Primarily used for choosing one of a symmetric pair as canonical */
static int loc_lt(loc_t s1, loc_t s2) {
    return s1.x < s2.x || (s1.x == s2.x && s1.y < s2.y);
}

static void list_set(loclist_t *ll, int i, loc_t val) {
    if (i >= ll->size)
        resize_loclist(ll, i * 3 / 2);
    ll->list[i] = val;
    if (ll->used <= i)
        ll->used = i + 1;
}

static void list_append(loclist_t *ll, loc_t val) {
    list_set(ll, ll->used, val);
}

static loc_t list_get(loclist_t *ll, int i) {
    assert(i < ll->used);
    return ll->list[i];
}

/* Return the index of this point in the loclist, or -1 */
static int list_find_lim(loclist_t *ll, int lim, loc_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc_eq(ll->list[i], val))
            return i;
    return -1;
}

static int list_find(loclist_t *ll, loc_t val) {
    return list_find_lim(ll, ll->used, val);
}

/* Return TRUE if this point is in the loclist, else FALSE */
static int list_exists_lim(loclist_t *ll, int lim, loc_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc_eq(ll->list[i], val))
            return 1;
    return 0;
}

static int list_exists(loclist_t *ll, loc_t val) {
    return list_exists_lim(ll, ll->used, val);
}

static void list_remove(loclist_t *ll, loc_t val) {
    int i = list_find(ll, val);
    if (i >= 0)
        ll->list[i] = ll->list[--ll->used];
}

static int loc2p_eq(loc2p_t s1, loc2p_t s2) {
    int adjust = s2.power - s1.power;
    loc_t p1 = s1.p, p2 = s2.p;
    if (adjust > 0) {
        return (p1.x << adjust) == p2.x && (p1.y << adjust) == p2.y;
    } else {
        return p1.x == (p2.x << -adjust) && p1.y == (p2.y << -adjust);
    }
}

static loc_t raise2p(loc_t p) {
    return (loc_t){ p.x << 1, p.y << 1 };
}

static loc_t convert2p(loc2p_t p2, int power) {
    int adjust = power - p2.power;

    if (adjust >= 0) {
        return (loc_t){ p2.p.x << adjust, p2.p.y << adjust };
    } else {
        loc_t p = p2.p;
        while (adjust < 0) {
            assert((p.x & 1) == 0 && (p.y & 1) == 0);
            p.x >>= 1;
            p.y >>= 1;
        }
        return p;
    }
}

static loc2p_t lift2p(loc2p_t p, int power) {
    return (loc2p_t){ convert2p(p, power), power };
}

static int is_loc2p_odd(loc2p_t p, int power) {
    assert(p.power <= power);
    if (p.power < power)
        return 0;
    return (p.p.x & 1) || (p.p.y & 1);
}

static void list2p_set(loc2plist_t *l2l, int i, loc2p_t val) {
    assert(i < l2l->size);
    l2l->list[i] = val;
    if (l2l->used <= i)
        l2l->used = i + 1;
}

static loc2p_t list2p_get2p(loc2plist_t *l2l, int i) {
    assert(i < l2l->used);
    return l2l->list[i];
}

static loc_t list2p_get(loc2plist_t *l2l, int i, int power) {
    loc2p_t p2 = list2p_get2p(l2l, i);
    return convert2p(p2, power);
}

/* Return the index of this point in the loc2plist, or -1 */
static int list2p_find_lim(loc2plist_t *l2l, int lim, loc2p_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc2p_eq(l2l->list[i], val))
            return i;
    return -1;
}

static int list2p_exists_lim(loc2plist_t *l2l, int lim, loc2p_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc2p_eq(l2l->list[i], val))
            return 1;
    return 0;
}

static int pair_eq(pair_t pair1, pair_t pair2) {
    return (loc_eq(pair1.p[0], pair2.p[0]) && loc_eq(pair1.p[1], pair2.p[1]))
        || (loc_eq(pair1.p[0], pair2.p[1]) && loc_eq(pair1.p[1], pair2.p[0]));
}

static void pair_set(pairlist_t *pl, int i, pair_t pair) {
    if (i >= pl->size)
        resize_pairlist(pl, i * 3 / 2);
    pl->list[i] = pair;
    if (pl->used <= i)
        pl->used = i + 1;
}

static void pair_append(pairlist_t *pl, pair_t pair) {
    pair_set(pl, pl->used, pair);
}

static pair_t pair_get(pairlist_t *pl, int i) {
    assert(i < pl->used);
    return pl->list[i];
}

static int pair_exists(pairlist_t *pl, pair_t pair) {
    for (int i = 0; i < pl->used; ++i)
        if (pair_eq(pl->list[i], pair))
            return 1;
    return 0;
}

#endif
