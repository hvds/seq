#include "total.h"
#include <stdlib.h>
#include <stdio.h>

#define ROW(r) ({ int _r = r; (_r * (_r + 1)) >> 1; })

total_t *new_total(void) {
	total_t *tp = (total_t *)malloc(sizeof(total_t));
	tp->q = (mpq_t *)NULL;
	tp->rows = 0;
	tp->size = 0;
	return tp;
}

void free_total(total_t *tp) {
	int i;
	for (i = 0; i < ROW(tp->rows); ++i) {
		mpq_clear(tp->q[i]);
	}
	free(tp->q);
	free(tp);
}

mpq_t *get_total(total_t *tp, int row, int col) {
	if (row < 0 || row >= tp->rows || col < 0 || col > row) {
		fprintf(stderr, "Illegal offset (%d:%d) in totals[%d]\n", row, col, tp->rows);
		exit(-1);
	}
	return &(tp->q[ROW(row) + col]);
}

void expand_total(total_t *tp, mpq_t base) {
	int row = tp->rows;
	int size = ROW(++tp->rows);
	int i;
	if (size > tp->size) {
		if (size < 10) size = 10;
		if (size < tp->size * 2) size = tp->size * 2;
		tp->q = (mpq_t *)realloc(tp->q, size * sizeof(mpq_t));
		tp->size = size;
	}
	for (i = 0; i <= row; ++i) {
		mpq_init(tp->q[ROW(row) + i]);
		if (i > 0)
			mpq_add(tp->q[ROW(row) + i], base, tp->q[ROW(row - 1) + i - 1]);
		else
			mpq_set(tp->q[ROW(row) + i], base);
	}
}
