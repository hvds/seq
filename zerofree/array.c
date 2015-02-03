#include <stdlib.h>
#include "array.h"

/*
 * Simple support for resizable arrays of a given element size.
 */

void init_array(array_t *a, uint objsize, uint size) {
    a->objsize = objsize;
    a->size = size;
    a->count = 0;
    a->array = malloc(objsize * size);
}

void free_array(array_t *a) {
    free(a->array);
}

inline void resize_array(array_t *a, uint size) {
    if (size > a->size) {
        uint newsize = (a->size * 3) / 2;
        ulong actual = (ulong)newsize * (ulong)a->objsize;
        a->array = (void*)realloc(a->array, actual);
        a->size = newsize;
    }
}

inline void *array_element(array_t *a, uint index) {
    return a->array + index * a->objsize;
}

