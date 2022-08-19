#include <stdlib.h>
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

int main(int argc, char **argv) {
    if (argc < 5)
        exit(1);
    mpz_init_set_str(zA, argv[1], 10);
    mpz_init_set_str(zD, argv[2], 10);
    in = strtol(argv[3], NULL, 10);
    mpz_init_set_str(zlimit, argv[4], 10);
    gmp_printf("try new_pell(%Zu, %Zu, %d, %Zu)\n", zA, zD, in, zlimit);
    mpz_init(zx);
    mpz_init(zy);
    _GMP_init();
    init_tau();
    init_rootmod();
    init_pell();
    new_pell(zA, zD, in, zlimit);
    while (next_pell(zx, zy)) {
        gmp_printf("solution (%Zd, %Zd)\n", zx, zy);
    }
    done_pell();
    done_rootmod();
    done_tau();
}
