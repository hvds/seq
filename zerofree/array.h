#ifndef ARRAY_H

#include <stdio.h>

typedef unsigned int uint;

/* Simple resizable array support */
typedef struct {
    void* array;
    uint objsize;
    uint size;
    uint count;
} array_t;

extern void init_array(array_t *a, uint objsize, uint size);
extern void free_array(array_t *a);

extern inline void resize_array(array_t *a, uint size) {
    if (size > a->size) {
        uint newsize = (a->size * 3) / 2;
        ulong actual = (ulong)newsize * (ulong)a->objsize;
        a->array = (void*)realloc(a->array, actual);
        a->size = newsize;
    }
}

extern inline void *array_element(array_t *a, uint index) {
    return a->array + index * a->objsize;
}

#endif
