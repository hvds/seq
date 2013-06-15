#ifndef _SETH_H
#define _SETH_H 1

#include <stdio.h>
#include <assert.h>
#include "part.h"
#include "vec.h"

typedef struct set_s {
	unsigned char p[NODES];
} set_t;

#ifndef IS_SET_C
#define SET_INLINE extern inline
#else
#define SET_INLINE inline
#endif

SET_INLINE void set_zero(set_t* set) {
	memset(set, 0, sizeof(set_t));
}

SET_INLINE void set_copy(set_t* source, set_t* dest) {
	memcpy(dest, source, sizeof(set_t));
}

SET_INLINE void set_merge(set_t* source, set_t* dest, uint offset) {
	uint i;
	for (i = 0; i < NODES; ++i)
		if (source->p[i])
			dest->p[i] = source->p[i] + offset;
}

SET_INLINE void set_init(set_t* set, vec_t* v) {
	uint i;
	for (i = 0; i < NODES; ++i)
		set->p[i] = vec_testbit(v, i);
}

SET_INLINE void set_append(set_t* dest, set_t* source, vec_t* v, uint size) {
	uint i;
	set_copy(source, dest);
	for (i = 0; i < NODES; ++i)
		if (vec_testbit(v, i)) {
			assert(dest->p[i] == 0);
			dest->p[i] = size;
		}
}

void fprint_set(FILE* stream, set_t* set);

#endif /* seth.h */
