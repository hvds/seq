#ifndef TOTAL_H
#define TOTAL_H 1

#include <gmp.h>

typedef struct total_s {
	mpq_t *q;
	int rows;
	int size;
} total_t;

extern total_t *new_total(void);
extern void free_total(total_t *tp);
extern mpq_t *get_total(total_t *tp, int row, int col);
extern void expand_total(total_t *tp, mpq_t base);

#endif
