#define IS_SET_C
#include "set.h"

typedef unsigned int uint;

/* assume 32 is the max we need for now */
static char* output_map = "*abcdefghijklmnopqrstuvwxyzABCDEF";

void fprint_set(FILE* stream, set_t* set) {
	uint i, j;
	vec_t* v;

	for (i = 0; i < NODES; ++i)
		fprintf(stream, "%c", output_map[set->p[i]]);
}

