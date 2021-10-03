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
