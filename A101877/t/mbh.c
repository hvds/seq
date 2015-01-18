#include <stdio.h>
#include <stdlib.h>

#include "mbh.h"

uint g_fail = 0;
uint g_test = 0;

int mbh_compare(bhp h, bhsize_t left, bhsize_t right) {
	int il = P2I(BHP(h)->heap[left]);
	int ir = P2I(BHP(h)->heap[right]);
	return (il < ir) ? -1 : (il == ir) ? 0 : +1;
}

void test_cycle(bhp h) {
	int i, j;
	void* v;

	i = mbh_size(h);
	if (i != 0) {
		++g_fail;
		printf("Error: expected empty heap to have size 0, got %d\n", i);
	}
	++g_test;
	for (i = 0; i < 100; ++i) {
		j = (i * 3) % 100;
		mbh_insert(h, I2P(j));
	}
	i = mbh_size(h);
	if (i != 100) {
		++g_fail;
		printf("Error: expected filled heap to have size 100, got %d\n", i);
	}
	++g_test;
	for (i = 0; i < 100; ++i) {
		v = mbh_shift(h);
		if (P2I(v) != i) {
			++g_fail;
			printf("Error: from mod 3 heap expected %d got %d\n", i, P2I(v));
		}
		++g_test;
	}
	i = mbh_size(h);
	if (i != 0) {
		++g_fail;
		printf("Error: expected emptied heap to have size 0, got %d\n", i);
	}
	++g_test;
}

void test_b(void) {
	bhp b;
	b = mbh_new((void*)NULL);
	test_cycle(b);
}

void test_a(void) {
	bhp a;
	int i;
	void* v;
	a = mbh_new((void*)NULL);
	test_cycle(a);
	for (i = 0; i < 3; ++i) {
		mbh_insert(a, I2P(i));
	}
	test_b();
	i = mbh_size(a);
	if (i != 3) {
		++g_fail;
		printf("Error: expected suspended heap to have size 3, got %d\n", i);
	}
	for (i = 0; i < 3; ++i) {
		v = mbh_shift(a);
		if (P2I(v) != i) {
			++g_fail;
			printf("Error: heap a corrupted by heap b: expected %d, got %d\n",
					i, P2I(v));
		}
		++g_test;
	}
}

int main(int argc, char** argv) {
	bhp a, b;
	int i, j;
	setup_mbh();
	test_a();
	teardown_mbh();

	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}

