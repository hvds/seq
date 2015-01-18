#ifndef VECTOR_H
#define VECTOR_H 1

#include <stdlib.h>

typedef struct vector_s {
	unsigned char *v;
	size_t size;
} vector_t;

extern vector_t *new_vector(size_t size);
extern void fill_vector(vector_t *v, int value);
extern void set_vector(vector_t *v, size_t offset, int value);
extern int test_vector(vector_t *v, size_t offset);
extern void free_vector(vector_t *v);
extern vector_t *merge_vector(vector_t *vin, ulong fac);
extern vector_t *copy_vector(vector_t *vin);

#endif
