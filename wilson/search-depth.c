#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include "mygmp.h"
#include "depth.h"

long clock_tick;

mpq_t r;        /* the rational we're testing */
mpq_t rone;     /* handy 1 */
mpq_t limit;    /* r + 1/r */

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

void init(int p, int q) {
    QINIT(&r, "r");
    mpq_set_ui(r, (ulong)p, (ulong)q);
    QINIT(&rone, "rone");
    mpq_set_ui(rone, (ulong)1, (ulong)1);
    QINIT(&limit, "limit");
    mpq_inv(limit, r);
    mpq_add(limit, limit, r);
}

void finish(void) {
    QCLEAR(&limit, "limit");
    QCLEAR(&rone, "rone");
    QCLEAR(&r, "r");
}

int main(int argc, char** argv) {
    int p, q, d;

    if (argc != 4) {
        fprintf(stderr, "Usage: search-depth <p> <q> <depth>\n");
        return 1;
    }
    p = atoi(argv[1]);
    q = atoi(argv[2]);
    d = atoi(argv[3]);
    if (!(p > 0 && q > p)) {
        fprintf(stderr, "Value error: need 0 < p/q < 1\n");
        return 1;
    }
    if (!(d > 1)) {
        fprintf(stderr, "Value error: need depth > 1\n");
        return 1;
    }

    clock_tick = sysconf(_SC_CLK_TCK);
    setlinebuf(stdout);
    init(p, q);
    init_depth((ulong)d);
    if (!search_depth())
        gmp_printf("%Qd - g > %d (%.2fs)\n", r, d, timing());
    finish_depth();
    finish();
    return 0;
}
