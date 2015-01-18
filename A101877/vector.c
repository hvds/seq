#include "vector.h"
#include <stdio.h>
#include <string.h>

#define MSIZE(z) ((z + 8) >> 3)
#define MBYTE(z) (z >> 3)
#define MOFF(z) (z & 7)

typedef unsigned char uchar;

vector_t *new_vector(size_t size) {
	vector_t *v = (vector_t *)malloc(sizeof(vector_t));
	v->v = (uchar *)malloc(MSIZE(size));
	v->size = size;
	return v;
}

void fill_vector(vector_t *v, int value) {
	size_t i;
	value = (value ? 0xff : 0);
	for (i = 0; i < MSIZE(v->size); ++i) v->v[i] = value;
}

void set_vector(vector_t *v, size_t offset, int value) {
	if (offset > v->size) {
		fprintf(stderr, "offset %zu out of range for vector size %zu\n", offset, v->size);
		exit(-1);
	}
	if (value) {
		v->v[MBYTE(offset)] |= (1 << MOFF(offset));
	} else {
		v->v[MBYTE(offset)] &= ~(1 << MOFF(offset));
	}
}

int test_vector(vector_t *v, size_t offset) {
	if (offset > v->size) {
		fprintf(stderr, "offset %zu out of range for vector size %zu\n", offset, v->size);
		exit(-1);
	}
	return (v->v[MBYTE(offset)] & (1 << MOFF(offset))) ? 1 : 0;
}

void free_vector(vector_t *v) {
	free(v->v);
	free(v);
}

/*
	given input vector vin, construct and return a new vector vout
	that is (vin | T(vin)) where T(vin) is vin rotated by mod bits.
*/
vector_t *merge_vector(vector_t *vin, ulong mod) {
	size_t p = vin->size;
	ulong size = MSIZE(p);
	ulong off, value, byte, inbyte;
	ulong tail = vin->v[0] + ((size > 1) ? (vin->v[1] << 8) : 0);
	vector_t *vout = new_vector(p);

	for (byte = 0; byte < size; ++byte) {
		off = (byte * 8 + p - mod) % p;
		inbyte = off / 8;
		value = vin->v[inbyte] + (
			(inbyte + 1 == size) ? 0 : (vin->v[inbyte + 1] << 8) 
		) + (
			(inbyte + 2 >= size) ? (
				tail << ((p & 7) + ((inbyte + 1 == size) ? 0 : 8))
			) : 0
		);
		vout->v[byte] = vin->v[byte] | (uchar)((value >> (off % 8)) & 255);
	}
	vout->v[size - 1] &= (1 << (p & 7)) - 1;
	return vout;
}

vector_t *copy_vector(vector_t *vin) {
	vector_t *vout = new_vector(vin->size);
	memcpy((void *)vout->v, (void *)vin->v, MSIZE(vin->size));
	return vout;
}
