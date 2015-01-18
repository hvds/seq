#include <stdio.h>
#include <stdlib.h>

#include "inverse.h"

uint g_fail = 0;
uint g_test = 0;

#define MARKER 0xffffffffu

void test_euclid(uint p) {
	uint i, j, k;
	for (i = 1; i < p; ++i) {
		j = inveuclid(i, p);
		k = (j * i) % p;
		if (k != 1) {
			printf("Error: inveuclid(%u, %u) = %u\n", i, p, j);
			++g_fail;
		}
		++g_test;
	}
}

void test_table(uint p) {
	uint i, j, k;
	uint* up;

	up = (uint*)malloc(sizeof(uint) * (p + 2));
	up[0] = up[p + 1] = MARKER;
	invtable(p, &(up[1]));
	if (up[0] != MARKER || up[p + 1] != MARKER) {
		printf("Error: array bounds overwrite for invtable(%u)\n", p);
		++g_fail;
	}
	++g_test;
	for (i = 1; i < p; ++i) {
		j = up[i + 1];
		k = (j * i) % p;
		if (k != 1) {
			printf("Error: invtable[%u][%u] = %u\n", p, i, j);
			++g_fail;
		}
		++g_test;
	}
	free(up);
}

void test_cache(uint p) {
	uint i, j, k;

	inverse_table(p);
	for (i = 1; i < p; ++i) {
		j = invfast(i, p);
		k = (j * i) % p;
		if (k != 1) {
			printf("Error: invfast(%u, %u) = %u\n", i, p, j);
			++g_fail;
		}
		++g_test;
	}
}

void test_all(uint p) {
	test_euclid(p);
	test_table(p);
	test_cache(p);
}

int main(int argc, char** argv) {
	setup_inverse();
	test_all(2);
	test_all(3);
	test_all(5);
	test_all(7);
	test_all(11);
	test_all(17);
	test_all(257);
	clear_inverse();
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}

