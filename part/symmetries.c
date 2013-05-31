#define IS_SYMMETRIES_C
#include "symmetries.h"

typedef unsigned int uint;
uint sym_count;
sym_t* symmetries;

int map_cmp(sym_t* ma, sym_t* mb) {
	uint i;
	int c;
	for (i = 0; i < NODES; ++i) {
		c = ma->map[i] - mb->map[i];
		if (c > 0)
			return 1;
		if (c < 0)
			return -1;
	}
	return 0;
}

void next_perm(uint* p, uint size) {
	uint last = size - 1;
	uint j, k;
	uint temp;
	int i = (int)last - 1;

	while (i >= 0 && p[i] > p[i+1])
		--i;
	if (i < 0)
		return;
	for (j = i + 1, k = last; j < k; ++j, --k) {
		temp = p[j];
		p[j] = p[k];
		p[k] = temp;
	}
	for (j = i + 1; p[j] < p[i]; ++j)
		;
	temp = p[i];
	p[i] = p[j];
	p[j] = temp;
}

void setup_symmetries(void) {
	uint fac = 1;
	uint i, j, k, value;
	sym_t *m, *m2;
	uint perm[NBASE];

	for (i = 2; i <= NBASE; ++i)
		fac *= i;
	sym_count = NODES * fac;
	symmetries = (sym_t*)malloc(sym_count * sizeof(sym_t));

	for (i = 0; i < NBASE; ++i) {
		perm[i] = i;
	}
	for (i = 0; i < fac; ++i, next_perm(perm, NBASE)) {
		m = sym_map(i);
		for (j = 0; j < NODES; ++j) {
			value = 0;
			for (k = 0; k < NBASE; ++k) {
				if (j & (1 << perm[k]))
					value |= 1 << k;
			}
			m->map[j] = value;
		}
	}

	m = sym_map(0);
	for (i = 1; i < NODES; ++i) {
		m2 = sym_map(i * fac);
		for (j = 0; j < NODES * fac; ++j)
			m2->map[j] = m->map[j] ^ i;
	}

	qsort(symmetries, sym_count, sizeof(sym_t), (__compar_fn_t)map_cmp);
}

void teardown_symmetries(void) {
	free(symmetries);
}
