#include "prime.h"

int gcd(int a, int b) {
	int temp;
	if (a > b) {
		temp = a;
		a = b;
		b = temp;
	}
	while (a > 0) {
		temp = b % a;
		b = a;
		a = temp;
	}
	return b;
}

/* assuming this gets called once for each n: 1 <= n <= 1250, a simplistic
 * approach is amply fast enough */
int greatest_prime_power(int n, int* prime) {
	int d, q, best = 0, bestp = 0, power;
	for (d = 2; d * d <= n; ++d) {
		q = n / d;
		if (d * q != n)
			continue;
		power = d;
		n = q;
		q = n / d;
		while (d * q == n) {
			power *= d;
			n = q;
			q = n / d;
		}
		if (best < power) {
			best = power;
			bestp = d;
		}
	}
	/* what's left is prime */
	if (best < n)
		best = bestp = n;
	if (prime)
		*prime = bestp;
	return best;
}

/* assuming this gets called once for each n: 1 <= n <= 1250, a simplistic
 * approach is amply fast enough */
int greatest_prime(int n) {
	int d, q, best = 0;
	for (d = 2; d * d <= n; ++d) {
		q = n / d;
		while (d * q == n) {
			best = d;
			n = q;
			q = n / d;
		}
	}
	/* what's left is prime */
	return (n > best) ? n : best;
}
