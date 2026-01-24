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
int debug_integrate = 0;
int debug_suppress_write = 0;
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

int main(int argc, char **argv, char **envp) {
    int i = 1;
    uint ri, pi, rec;
    bool resolved = 1;
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        else if (strcmp(arg, "-r") == 0)
            resolved = 1 - resolved;
        else
            fail("unknown option '%s'", arg);
    }
    if (i + 5 == argc) {
        na = strtoul(argv[i++], NULL, 10);
        if (na < 1)
            fail("require a >= 1, not %lu", na);
        nb = strtoul(argv[i++], NULL, 10);
        if (nb < 1)
            fail("require b >= 1, not %lu", nb);
        ri = strtoul(argv[i++], NULL, 10);
        pi = strtoul(argv[i++], NULL, 10);
        rec = strtoul(argv[i++], NULL, 10);
    } else
        fail("wrong number of arguments");

    init_post();
    mmfrag_init();
    uint count = mmfrag_open(resolved, ri, pi);
    if (rec >= count)
        fail("cannot fetch record %u of %u\n", rec, count);
    frag_t *fp = mmfrag_get(rec);
    if (fp == NULL)
        fail("could not read frag\n");
    fid_t f = new_frag();
    memcpy(frag_p(f), fp, frag_size());
    frag_dump(f);
    mmfrag_done();

    done();
    return 0;
}
