#include <stdlib.h>
#include <string.h>
#include "loc.h"

loclist_t *new_loclist(int size) {
    loclist_t *ll = (loclist_t *)malloc(sizeof(loclist_t));
    ll->size = size;
    ll->used = 0;
    ll->list = (loc_t *)malloc(sizeof(loc_t) * size);
    return ll;
}

void free_loclist(loclist_t *ll) {
    free(ll->list);
    free(ll);
}

void resize_loclist(loclist_t *ll, int size) {
    ll->size = size;
    ll->list = realloc(ll->list, sizeof(loc_t) * size);
}

loclist_t *dup_loclist(loclist_t *ll) {
    loclist_t *ll2 = new_loclist(ll->used + 10);
    ll2->used = ll->used;
    memcpy(ll2->list, ll->list, sizeof(loc_t) * ll->used);
    return ll2;
}

pairlist_t *new_pairlist(int size) {
    pairlist_t *pl = (pairlist_t *)malloc(sizeof(pairlist_t));
    pl->size = size;
    pl->used = 0;
    pl->list = (pair_t *)malloc(sizeof(pair_t) * size);
    return pl;
}

void free_pairlist(pairlist_t *pl) {
    free(pl->list);
    free(pl);
}

void resize_pairlist(pairlist_t *pl, int size) {
    pl->size = size;
    pl->list = realloc(pl->list, sizeof(pair_t) * size);
}

pairlist_t *dup_pairlist(pairlist_t *pl) {
    pairlist_t *pl2 = new_pairlist(pl->used + 10);
    pl2->used = pl->used;
    memcpy(pl2->list, pl->list, sizeof(pair_t) * pl->used);
    return pl2;
}
