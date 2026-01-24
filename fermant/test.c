/* needed on FreeBSD for getline() to be exported from stdio.h */
#define _WITH_GETLINE

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <signal.h>
#include <time.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "int.h"
#include "path.h"
#include "frag.h"
#include "calc.h"
#include "source.h"
#include "../inc/diag.h"

/* primary parameters - we are searching for g(na, nb), the expected cost
 * of traversing an (a x b) grid where the cost of each segment varies
 * uniformly between 0 and 1.
 */
uint na, nb;
/* there are nv = na(nb-1) + nb(na-1) variables to integrate over */
uint nv;

int debug_split = 0;
int debug_integrate = 1;
int debug_suppress_write = 1;
uint sync_count = 0;

double t0 = 0;
volatile bool need_diag = 0, need_log = 0;

void report(char *format, ...) {
    keep_diag();
    va_list ap;
    va_start(ap, format);
    gmp_vfprintf(stdout, format, ap);
    va_end(ap);
}

double seconds(double t1) {
    return (t1 - t0);
}

double elapsed(void) {
    return seconds(utime());
}

void done(void) {
    done_calc();
    done_frags();
    done_paths();
}

void fail_silent(void) {
    exit(0);
}
void fail(char *format, ...) {
    va_list ap;
    va_start(ap, format);
    gmp_vfprintf(stderr, format, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

void init_pre(void) {
    t0 = utime();
}

void init_post(void) {
    nv = 2 * na * nb - na - nb;

    /* generate list of paths and resolutions
     * TODO: add strategies
     */
    init_paths(0);
    init_frags();
    init_calc(npaths);

    init_diag();    /* ignore result: worst case we lose ^Z handling */
}

void parse_limit(fid_t f, uint vi, bool high, char *s, limitz_t num, limitz_t den) {
    limit_t *lp = (high)
        ? range_high(frag_range(f, vi))
        : range_low(frag_range(f, vi));
    limitp_num_set(lp, num);
    limitp_den_set(lp, den);
    lincom_t *lc = limitp_lc(lp);
    memset(lc, 0, lc_size());
    while (*s) {
        int sign = 1;
        int coeff = 1;
        int var = 0;
        while (*s == '+' || *s == '-') {
            if (*s == '-')
                sign = -sign;
            ++s;
        }
        if (isdigit(*s))
            coeff = strtol(s, &s, 10);
        if (isalpha(*s)) {
            var = *s++ - ('a' - 1);
            if (var < 1 || var > nv)
                fail("parse fail: c = '%c'\n", var + ('a' - 1));
        }
        lc_set(lc, var, (lincomz_t)sign * coeff);
        if (*s && (*s != '+' && *s != '-'))
            fail("parse fail at '%s'\n", s);
    }
}

void init_test(fid_t f) {
    frag_ps_set(f, (pathset_t)1);
    parse_limit(f, 1, 0, "1", 0, 1);
    parse_limit(f, 1, 1, "1", 1, 2);
    parse_limit(f, 2, 0, "1", 0, 1);
    parse_limit(f, 2, 1, "1-2a", 1, 2);
    parse_limit(f, 3, 0, "1-2a-2b", 1, 2);
    parse_limit(f, 3, 1, "1-a-b", 1, 1);
    parse_limit(f, 4, 0, "b+c", 1, 1);
    parse_limit(f, 4, 1, "1-a", 1, 1);
    parse_limit(f, 5, 0, "1-2a-b-c-d", -1, 2);
    parse_limit(f, 5, 1, "1-3a-2b-2c-d", -1, 3);
    parse_limit(f, 6, 0, "1-8194a-2b-2c+2e", 1, 1);
    parse_limit(f, 6, 1, "1-a-2b-2c+d+e", 1, 2);
    parse_limit(f, 7, 0, "1-f", 1, 1);
    parse_limit(f, 7, 1, "1+a+2b+2c-d-e-f", 1, 2);
    parse_limit(f, 8, 0, "1-f-g", 1, 1);
    parse_limit(f, 8, 1, "1-g", 1, 1);
    parse_limit(f, 9, 0, "b+c-d+h", 1, 1);
    parse_limit(f, 9, 1, "a+b+c-e-f+h", 1, 1);
    parse_limit(f, 10, 0, "1", 0, 1);
    parse_limit(f, 10, 1, "a+b+c-e-g+i", 1, 1);
    parse_limit(f, 11, 0, "a+b+c-e-g+i-j", 1, 1);
    parse_limit(f, 11, 1, "f+i-j", 1, 1);
    parse_limit(f, 12, 0, "a+b+c-e+h-j-k", 1, 1);
    parse_limit(f, 12, 1, "g+h-i", 1, 1);
    return;
}
    
int main(int argc, char **argv, char **envp) {
    int i = 1;
    na = 3;
    nb = 3;
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        else
            fail("unknown option '%s'", arg);
    }
    if (i == argc) {
        ;
    } else
        fail("wrong number of arguments");

    init_post();
    fid_t f = new_frag();
    init_test(f);
    frag_dump(f);
extern void integrate(fid_t fi);
    integrate(f);

    done();
    return 0;
}
