#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pell.h"
#include "rootmod.h"
#include "coultau.h"
#include "gmp_main.h"

mpz_t zA;
mpz_t zD;
int in;
mpz_t zlimit;
mpz_t zx;
mpz_t zy;

/* needed for coultau.c; but not used, since we init_tau(0) */
#include "coul.h"
t_divisors *divisors = NULL;

void ston(mpz_t targ, char *s) {
    char *t = strchr(s, 'e');
    if (t) {
        *t = 0;
        mpz_set_str(targ, s, 10);
        ulong exp = strtoul(&t[1], NULL, 10);
        mpz_t temp;
        mpz_init(temp);
        mpz_ui_pow_ui(temp, 10, exp);
        mpz_mul(targ, targ, temp);
        mpz_clear(temp);
        *t = 'e';
    } else {
        mpz_set_str(targ, s, 10);
    }
}

/* Run this with 4 arguments (A, D, N, lim) to test the Pell solver.
 * We will search for solutions to Ax^2 - Dy^2 = N with 0 <= x <= lim.
 * A, D, lim are each positive integers, of arbitrary size;
 * N is a signed int, so must typically fit in 32 bits.
 * lim can be given in exponential format (without decimal), eg "1e20".
 */
int main(int argc, char **argv) {
    if (argc < 5)
        exit(1);
    mpz_init_set_str(zA, argv[1], 10);
    mpz_init_set_str(zD, argv[2], 10);
    in = strtol(argv[3], NULL, 10);
    mpz_init(zlimit);
    ston(zlimit, argv[4]);
    gmp_printf("try new_pell(%Zu, %Zu, %d, %Zu)\n", zA, zD, in, zlimit);
    mpz_init(zx);
    mpz_init(zy);
    _GMP_init();
    init_tau(0);
    init_rootmod(1);
    init_pell();
    new_pell(zA, zD, in, zlimit);
    while (next_pell(zx, zy)) {
        gmp_printf("solution (%Zd, %Zd)\n", zx, zy);
    }
    done_pell();
    done_rootmod();
    done_tau();
}
