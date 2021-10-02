#ifndef LOC_H
#define LOC_H

#include <assert.h>

typedef struct {
    int x;
    int y;
} loc_t;

typedef struct {
    int size;
    int used;
    loc_t *list;
} loclist_t;

loclist_t *new_loclist(int size);
void free_loclist(loclist_t *ll);
void resize_loclist(loclist_t *ll, int size);

static loc_t loc_diff(loc_t s1, loc_t s2) {
    loc_t d;
    d.x = s2.x - s1.x;
    d.y = s2.y - s1.y;
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

static int list_find_lim(loclist_t *ll, int lim, loc_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc_eq(ll->list[i], val))
            return i;
    return -1;
}

static int list_find(loclist_t *ll, loc_t val) {
    return list_find_lim(ll, ll->used, val);
}

static int list_exists_lim(loclist_t *ll, int lim, loc_t val) {
    for (int i = 0; i < lim; ++i)
        if (loc_eq(ll->list[i], val))
            return 1;
    return 0;
}

static int list_exists(loclist_t *ll, loc_t val) {
    return list_exists_lim(ll, ll->used, val);
}

#endif
