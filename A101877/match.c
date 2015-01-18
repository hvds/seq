#include "match.h"
#ifdef DEBUG
#  include <stdio.h>
#endif

match_t *new_match(void) {
	match_t *m = (match_t *)malloc(sizeof(match_t));
	mpq_init(m->sum);
	m->v = (vector_t *)NULL;
	return m;
}

void free_match(match_t *m) {
	mpq_clear(m->sum);
	if (m->v) free_vector(m->v);
	free(m);
}

match_t *copy_match(match_t *m) {
	match_t *m2 = new_match();
	mpq_set(m2->sum, m->sum);
	if (m->v) m2->v = copy_vector(m->v);
	return m2;
}
