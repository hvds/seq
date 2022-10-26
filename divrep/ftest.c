#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "coul.h"
#include "coultau.h"

#include "factor.h"
#include "gmp_main.h"
#include "utility.h"
#include "primality.h"

extern bool tau_single_try(uint i);
t_divisors *divisors = NULL;
double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

double seconds(double t1) {
    return t1 - t0;
}

int main(int argc, char **argv) {
    _GMP_init();
    init_tau(0);
    alloc_taum(1);
    t0 = utime();

    t_tm *tm = &taum[0];
    /* some tests assume this already set */
    /* TODO: calculate based on n */
    tm->B1 = 160000;

    if (argc > 1)
        mpz_set_str(tm->n, argv[1], 10);
    else {
        gmp_fprintf(stderr, "Usage: %s <n> <test> <seed> [count]\n");
        return 1;
    }

    uint test = 31;
    if (argc > 2)
        test = strtoul(argv[2], NULL, 0);

    ulong randseed = 1;
    if (argc > 3)
        randseed = strtoul(argv[3], NULL, 10);
    uint count = 1;
    if (argc > 4)
        count = strtoul(argv[4], NULL, 10);

    tm->t = 4;
    tau_multi_prep(0);
    for (uint i = 1; i <= count; ++i) {
        gmp_printf("Test %u on %Zu seed %lu (%.2f)\n", test, tm->n, randseed,
                seconds(utime()));
        clear_randstate();
        init_randstate(randseed++);
        tm->bits = 1UL << test;
        if (!tau_single_try(0))
            continue;
        gmp_printf("Found factor %Zu (%.2f)\n", tm->n, seconds(utime()));
        return 0;
    }
    gmp_printf("Nothing found (%.2f)\n", seconds(utime()));
    return 1;
}
