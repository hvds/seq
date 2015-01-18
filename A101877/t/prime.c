#include "prime.h"
#include <stdio.h>

int g_fail = 0;
int g_test = 0;

void testgcd(int a, int b, int expect) {
	int result = gcd(a, b);
	if (result != expect) {
		printf("Error: for gcd(%d, %d) expected %d, got %d\n",
				a, b, expect, result);
		++g_fail;
	}
	++g_test;
}

void test_gpp(int n, int expectpower, int expectprime) {
	int result = greatest_prime_power(n, (int*)NULL);
	int prime = -1;
	if (result != expectpower) {
		printf("Error: for greatest_prime_power(%d) expected %d, got %d\n",
				n, expectpower, result);
		++g_fail;
	}
	++g_test;
	result = greatest_prime_power(n, &prime);
	if (result != expectpower || prime != expectprime) {
		printf("Error: for greatest_prime_power(%d) expected (%d, %d), got (%d, %d)\n",
				n, expectpower, expectprime, result, prime);
		++g_fail;
	}
	++g_test;
}

void test_gp(int n, int expect) {
	int result = greatest_prime(n);
	if (result != expect) {
		printf("Error: for greatest_prime(%d) expected %d, got %d\n",
				n, expect, result);
		++g_fail;
	}
	++g_test;
}

void test_gppair(int n, int expectgp, int expectgpp, int expectgpprime) {
	test_gp(n, expectgp);
	test_gpp(n, expectgpp, expectgpprime);
}

int main(int argc, char** argv) {
	int i, j;
	int fib0 = 1, fib1 = 2, fib2;
	for (i = 1; i < 10; ++i) {
		testgcd(1, i, 1);
		testgcd(i, 1, 1);
		testgcd(i, 2 * i + 1, 1);
		testgcd(2, 2 * i, 2);
		testgcd(i, 2 * i, i);
		fib2 = fib0 + fib1;
		fib0 = fib1;
		fib1 = fib2;
		testgcd(fib0, fib1, 1);
	}

	for (j = 2; j <= 3; ++j) {
		int pow = 1;
		for (i = 1; i < 10; ++i) {
			pow *= j;
			test_gppair(pow, j, pow, j);
		}
	}
	test_gppair(121, 11, 121, 11);
	test_gppair(1237, 1237, 1237, 1237);
	test_gppair(78, 13, 13, 13);
	test_gppair(234, 13, 13, 13);
	test_gppair(702, 13, 27, 3);
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
