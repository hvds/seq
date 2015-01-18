#include <stdio.h>
#include "prime.h"

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

void test_gppz(mpz_t z, int expectpower, int expectprime) {
	int result = z_greatest_prime_power(z, (int*)NULL);
	int prime = -1;
	if (result != expectpower) {
		printf("Error: for greatest_prime_power(");
		mpz_out_str(stdout, 10, z);
		printf(") expected %d, got %d\n", expectpower, result);
		++g_fail;
	}
	++g_test;
	result = z_greatest_prime_power(z, &prime);
	if (result != expectpower || prime != expectprime) {
		printf("Error: for greatest_prime_power(");
		mpz_out_str(stdout, 10, z);
		printf(") expected (%d, %d), got (%d, %d)\n",
				expectpower, expectprime, result, prime);
		++g_fail;
	}
	++g_test;
}

int main(int argc, char** argv) {
	int i, j;
	int fib0 = 1, fib1 = 2, fib2;
	mpz_t z;
	ZINIT(&z, "test temp");
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
			test_gpp(pow, pow, j);
		}
	}
	test_gpp(121, 121, 11);
	test_gpp(1237, 1237, 1237);
	test_gpp(78, 13, 13);
	test_gpp(234, 13, 13);
	test_gpp(702, 27, 3);
	mpz_set_ui(z, 2 * 3 * 5 * 7 * 11 * 13 * 17);
	mpz_mul_ui(z, z, 19 * 23 * 29 * 31);
	test_gppz(z, 31, 31);
	mpz_mul_ui(z, z, 19 * 19 * 31);
	test_gppz(z, 19 * 19 * 19, 19);

	ZCLEAR(&z, "test temp");
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
