#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <sys/times.h>

#include "coul.h"
#include "coultau.h"

#include "factor.h"
#include "gmp_main.h"
#include "utility.h"
#include "primality.h"

extern bool tau_single_try(uint i);
t_divisors *divisors = NULL;
long ticks_per_second;
clock_t ticks = 0;
struct tms time_buf;
static inline clock_t utime(void) {
    times(&time_buf);
    return time_buf.tms_utime;
}

double seconds(clock_t t1) {
    return (double)(t1 - ticks) / ticks_per_second;
}

int main(int argc, char **argv) {
    _GMP_init();
    init_tau(0);
    alloc_taum(1);
    ticks_per_second = sysconf(_SC_CLK_TCK);
    ticks = utime();

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
