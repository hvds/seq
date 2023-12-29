#include <stdlib.h>
#include <string.h>
#include "sieve.h"

uint *prime = NULL;
uint nprime = 0;
uint primesize = 0;
uint maxsieve = 0;

char *sieve = NULL;
uint sievesize = 0;

static inline void SETBIT(uint off) {
    sieve[off >> 3] |= 1 << (off & 7);
}

static inline bool HASBIT(uint off) {
    return (sieve[off >> 3] & (1 << (off & 7))) ? 1 : 0;
}

void store_prime(uint p) {
    if (nprime >= primesize) {
        primesize += 4096;
        prime = realloc(prime, sizeof(uint) * primesize);
    }
    prime[nprime++] = p;
}

void extend_sieve(uint limit) {
    if (limit <= maxsieve)
        return;

    uint sievez = (limit - maxsieve + 7) / 8;
    if (sievez > sievesize) {
        sieve = realloc(sieve, sievez);
        sievesize = sievez;
    }
    memset(sieve, 0, sievez);

    bool large = 0;
    for (uint pi = 0; pi < nprime; ++pi) {
        uint p = prime[pi];
        uint t = p * p;
        if (t < maxsieve)
            t = maxsieve - (maxsieve % p) + p;
        if (t > limit) {
            large = 1;
            break;
        }
        while (t <= limit) {
            SETBIT(t - maxsieve);
            t += p;
        }
    }

    for (uint try = (maxsieve == 0) ? 2 : maxsieve + 1; try <= limit; ++try) {
        if (HASBIT(try - maxsieve))
            continue;
        store_prime(try);
        if (large)
            continue;
        /* sieve the newly found prime */
        uint t = try * try;
        if (t > limit) {
            large = 1;
            continue;
        }
        while (t <= limit) {
            SETBIT(t - maxsieve);
            t += try;
        }
    }

    maxsieve = limit;
}
