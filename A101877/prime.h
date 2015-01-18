#ifndef PRIME_H
#define PRIME_H

#include "mygmp.h"

extern int gcd(int a, int b);
extern int greatest_prime_power(int n, int* prime);
extern int z_greatest_prime_power(mpz_t z, int* prime);

#endif /* PRIME_H */
