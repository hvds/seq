#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "clock.h"

/* we rely on mask_t to be an integer type
 * with at least 2 ** MAX_DIMENSION bits
 */
#define MAX_DIMENSION (6)
typedef unsigned long long mask_t;

mask_t fullmask[MAX_DIMENSION + 1] = {
    0x0000000000000001,
    0x0000000000000003,
    0x000000000000000f,
    0x00000000000000ff,
    0x000000000000ffff,
    0x00000000ffffffff,
    0xffffffffffffffff
};
unsigned int step = 0;
clock_t t0;

mpz_t sum[MAX_DIMENSION + 1];
mpz_t scratch[MAX_DIMENSION + 1];

void init_gmp(void) {
    int i;
    for (i = 0; i <= MAX_DIMENSION; ++i) {
        mpz_init(sum[i]);
        mpz_init(scratch[i]);
    }
}

void clear_gmp(void) {
    int i;
    for (i = 0; i <= MAX_DIMENSION; ++i) {
        mpz_clear(sum[i]);
        mpz_clear(scratch[i]);
    }
}

/*
 * count_part(d, mask): set sum[d] to the number of ways the specified
 * shape can be partitioned into polyominoes of size at most 2.
 * The shape is some or all of a I<d>-dimensional hypercube of side 2;
 * I<mask> is a vector of C<2 ** d> bits specifying the unit I<d>-cubes
 * that are included in the shape.
 *
 * Overwrites scratch[0..d] and sum[0..d] during calculation.
 */
void count_part(int d, mask_t mask) {
    mask_t left, right, shared, u;
    if (d == 0) {
        mpz_set_ui(sum[d], 1);
        return;
    }
    left = mask & fullmask[d - 1];
    right = mask >> (1 << (d - 1));
    shared = left & right;

    u = 0;
    mpz_set_ui(sum[d], 0);
    while (1) {
        count_part(d - 1, left ^ u);
        mpz_set(scratch[d], sum[d - 1]);
        count_part(d - 1, right ^ u);
        mpz_mul(scratch[d], scratch[d], sum[d - 1]);
        mpz_add(sum[d], sum[d], scratch[d]);
        u = (u - shared) & shared;
        if (u == 0)
            break;
    }
}

/*
 * count_full(d): same as count_part(d, fullmask[d])
 */
void count_full(int d) {
    mask_t lim, u;
    if (d == 0) {
        mpz_set_ui(sum[d], 1);
        return;
    }
    lim = fullmask[d - 1];
    mpz_set_ui(sum[d], 0);
    for (u = 0; u <= lim; ++u) {
        count_part(d - 1, u);
        mpz_mul(scratch[d], sum[d - 1], sum[d - 1]);
        mpz_add(sum[d], sum[d], scratch[d]);
        if ((++step & 0xfffff) == 0) {
            gmp_printf("300 at (%d, %llx) sum is %Zu (%.2f)\n",
                    d, u, sum[d], difftime(t0, curtime()));
        }
    }
}

#define LIM4 (1 << (1 << 4))
unsigned int cache4[LIM4];  /* the sum of these is 41025, int is ample room */
void init_cache4(void) {
    int u;
    for (u = 0; u < LIM4; ++u) {
        count_part(4, (mask_t)u);
        cache4[u] = mpz_get_ui(sum[4]);
    }
}

/*
 * count_5_part(mask): same as count_part(5, mask)
 *
 * Uses a prepared cache of results for d=4, to minimize recursion.
 */
void count_5_part(mask_t mask) {
    mask_t left, right, shared, u;
    left = mask & fullmask[4];
    right = mask >> (1 << 4);
    shared = left & right;
    u = 0;
    mpz_set_ui(sum[5], 0);
    while (1) {
        mpz_add_ui(sum[5], sum[5], cache4[left ^ u] * cache4[right ^ u]);
        u = (u - shared) & shared;
        if (u == 0)
            break;
    }
}

/*
 * count_6_full(): same as count_full(6)
 *
 * Uses a prepared cache of results for d=4 (via count_5_part()) to reduce
 * recursion at a reasonable memory cost.
 */
void count_6_full(void) {
    mask_t lim, u;
    init_cache4();
    lim = fullmask[5];
    mpz_set_ui(sum[6], 0);
    for (u = 0; u <= lim; ++u) {
        count_5_part(u);
        mpz_mul(scratch[6], sum[5], sum[5]);
        mpz_add(sum[6], sum[6], scratch[6]);
        if ((++step & 0x3ffffff) == 0) {
            gmp_printf("300 at (%d, %llx) sum is %Zu (%.2f)\n",
                    6, u, sum[6], difftime(t0, curtime()));
        }
    }
}

/*
 * Calculates the number of ways a I<d>-dimensional hypercube of side 2
 * can be partitioned into polyominoes of size at most 2.
 */
int main(int argc, char** argv) {
    int d;

    if (argc != 2) {
        fprintf(stderr, "usage: %s <dimension>\n", argv[0]);
        return 1;
    }
    d = atoi(argv[1]);
    if (d < 0 || d > MAX_DIMENSION) {
        fprintf(stderr, "dimension must be in the range 0 to %u\n",
                MAX_DIMENSION);
        return 1;
    }

    setup_clock();
    t0 = curtime();
    init_gmp();
    (d == 6) ? count_6_full() : count_full(d);
    gmp_printf("200 a(2, %d) = %Zu (%.2f)\n",
            d, sum[d], difftime(t0, curtime()));
    clear_gmp();
    return 0;
}
