#include "prime.h"
#include "mygmp.h"

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

/*
 * Find the greatest prime power p^k that divides n.
 * Input:
 *   int n: the number to test
 *   int* prime: optionally pointer to an int to store p, else NULL
 * Returns:
 *   int p^k
 * Notes:
 *   If int* prime is not NULL, the base prime p will be stored there.
 *   The approach is very simplistic, this routine should not be used in
 * any tight loop.
 */
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

/*
 * Find the greatest prime power p^k that divides z.
 * Input:
 *   mpz_t z: the number to test.
 *   int* prime: optional pointer to an int, else NULL.
 * Returns:
 *   int p^k, the greatest prime power that divides z.
 * Notes:
 *   It is required that p^k fit in an int.
 *   If int* prime is not NULL, the base prime p will be stored there.
 *   The approach is very simplistic, this routine should not be used in
 * any tight loop.
 */
int z_greatest_prime_power(mpz_t z, int* prime) {
	mpz_t d, q, r, n;
	int best = 0, bestp = 0, power;

	ZINIT(&n, "z_greatest_prime_power n");
	ZINIT(&d, "z_greatest_prime_power d");
	ZINIT(&q, "z_greatest_prime_power q");
	ZINIT(&r, "z_greatest_prime_power r");
	mpz_set(n, z);
	mpz_set_ui(d, 2);
	while (1) {
		mpz_fdiv_qr(q, r, n, d);
		if (mpz_get_ui(r) == 0) {
			power = mpz_get_ui(d);
			mpz_set(n, q);
			mpz_fdiv_qr(q, r, n, d);
			while (mpz_get_ui(r) == 0) {
				power *= mpz_get_ui(d);
				mpz_set(n, q);
				mpz_fdiv_qr(q, r, n, d);
			}
			if (best < power) {
				best = power;
				bestp = mpz_get_ui(d);
			}
		}
		mpz_nextprime(d, d);
		mpz_mul(q, d, d);
		if (mpz_cmp(q, n) > 0)
			break;
	}
	/* what's left is prime */
	if (best < mpz_get_ui(n))
		best = bestp = mpz_get_ui(n);
	if (prime)
		*prime = bestp;
	ZCLEAR(&r, "z_greatest_prime_power r");
	ZCLEAR(&q, "z_greatest_prime_power q");
	ZCLEAR(&d, "z_greatest_prime_power d");
	ZCLEAR(&n, "z_greatest_prime_power n");
	return best;
}
