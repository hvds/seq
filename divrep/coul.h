#ifndef COUL_H
#define COUL_H

/* 'divisors[i].div' is a list of divisors of 'i' in descending order of
 * highest prime factor, then ascending. 'high' is the highest prime
 * factor of 'i'; 'alldiv' is the number of factors; 'highdiv' is the
 * number of factors that are a multiple of 'high'; 'sumpm' is sum{p_j - 1}
 * of the primes dividing 'i', with multiplicity; 'gcddm' is gcd{d_j - 1}
 * of the divisors of i.
 * Eg divisors[18] = {
 *   alldiv = 6, highdiv = 4, high = 3, sumpm = 5 = (3-1)+(3-1)+(2-1),
 *   gcddm = 1, div = = [3, 6, 9, 18, 2, 1]
 * }
 * whereas divisors[65].gcddm = gcd(1-1, 5-1, 13-1, 65-1) = 4.
 */
typedef struct s_divisors {
    uint alldiv;    /* number of divisors of i */
    uint highdiv;   /* number of divisors that are multiples of 'high' */
    uint high;      /* highest prime dividing i */
    uint sumpm;     /* sum{p_j - 1} of primes dividing i /*/
    uint gcddm;     /* gcd{d_j - 1} of divisors of i */
    uint *div;      /* array of divisors of i */
} t_divisors;
extern t_divisors *divisors;

#endif
