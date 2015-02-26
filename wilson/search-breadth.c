#include <unistd.h>
#include <sys/times.h>
#include "mygmp.h"
#include "breadth.h"

long clock_tick;

mpq_t r;        /* the rational we're testing */
mpq_t rone;     /* handy 1 */
mpq_t report_calc;  /* scratch space for reporting solution */

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
    QINIT(&report_calc, "report_calc");
}

void finish(void) {
    QCLEAR(&report_calc, "report_calc");
    QCLEAR(&rone, "rone");
    QCLEAR(&r, "r");
}

void report_breadth(ulong gen, ulong count) {
    mpq_set(report_calc, r);
    mpz_mul_ui(mpq_numref(report_calc), mpq_numref(report_calc), count);
    mpq_canonicalize(report_calc);
    mpq_add(report_calc, report_calc, rone);
    gmp_printf("%Qd is solved in %lu steps starting with %lu -> %Qd\n",
            r, gen, count, report_calc);
}

extern ulong best_bits_num, best_bits_den; /* simplest in a generation */
void search_breadth(void) {
    ulong gen = 0;
    init_breadth(r);
    /* Each value processed will queue at least one new value,
     * so we don't need to check for exhaustion.
     */
    while (1) {
        rat_array_t *cur;
        ++gen;
        cur = choose_cur(gen);
        gmp_printf("%Qd - g > %lu: %lu (%.2fs) size %lu, best_bits %lu/%lu\n",
            r, gen, (ulong)cur->actual, timing(),
            (ulong)cur->count, best_bits_num, best_bits_den
        );
        if (breadth_one(gen))
            break;
    }
    finish_breadth();
}

int main(int argc, char** argv) {
    int p, q;

    if (argc != 3) {
        fprintf(stderr, "Usage: search-breadth <p> <q>\n");
        return 1;
    }
    p = atoi(argv[1]);
    q = atoi(argv[2]);
    if (!(p > 0 && q > p)) {
        fprintf(stderr, "Value error: need 0 < p/q < 1\n");
        return 1;
    }

    clock_tick = sysconf(_SC_CLK_TCK);
    init(p, q);
    search_breadth();
    finish();
    return 0;
}
