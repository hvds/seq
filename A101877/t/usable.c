#include "usable.h"
#include <stdio.h>

int g_fail = 0;
int g_test = 0;

void test_usable(int p, int pp, int k, int expect) {
	int result = is_usable(p, pp, k);
	if (result != expect) {
		printf("Error: for is_usable(%d, %d, %d) expected %d, got %d\n",
				p, pp, k, expect, result);
		++g_fail;
	}
	++g_test;
}

void test_cusp(int p, int pp, int m) {
	test_usable(p, pp, pp * m, 1);
	test_usable(p, pp, pp * m - 1, 0);
}

int main(int argc, char** argv) {
	int i;
	init_usable(1250);
	for (i = 1; i <= 9; ++i) {
		test_usable(2, 1 << i, 1250, 1);
	}
	test_usable(2, 1 << 10, 1250, 0);
	test_cusp(137, 137, 5);
	test_cusp(139, 139, 6);
	free_usable();
	init_usable(4096);	/* try to force realloc */
	test_cusp(2, 2048, 2);
	free_usable();
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
