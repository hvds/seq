#ifndef MATCH_H
#define MATCH_H 1

#include "vector.h"
#include <gmp.h>

typedef struct match_s {
	mpq_t sum;
	vector_t* v;
} match_t;

extern match_t *new_match(void);
extern void free_match(match_t *m);
extern match_t *copy_match(match_t *m);

#endif
