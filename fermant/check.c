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
#include "diag.h"

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

/* set to utime at start of run, minus last timestamp of recovery file */
double t0 = 0;
struct rusage rusage_buf;
extern inline double utime(void);

timer_t diag_timerid;
volatile bool need_diag = 0, need_log = 0;
bool clock_is_realtime = 0;
#define DIAG 1
double diag_delay = DIAG, diagt;

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

void handle_sig(int sig) {
    if (sig == SIGUSR1)
        need_diag = 1;
    else
        need_log = 1;
}

void init_time(void) {
    struct sigaction sa;
    struct sigevent sev;
    struct itimerspec diag_timer, log_timer;

    sa.sa_handler = &handle_sig;
    sa.sa_flags = SA_RESTART;
    sigemptyset(&sa.sa_mask);

    if (diag_delay) {
        if (sigaction(SIGUSR1, &sa, NULL))
            fail("Could not set USR1 handler: %s\n", strerror(errno));
        sev.sigev_notify = SIGEV_SIGNAL;
        sev.sigev_signo = SIGUSR1;
        sev.sigev_value.sival_ptr = &diag_timerid;
        if (timer_create(CLOCK_PROCESS_CPUTIME_ID, &sev, &diag_timerid)) {
            /* guess that the CPUTIME clock is not supported */
            if (timer_create(CLOCK_REALTIME, &sev, &diag_timerid))
                fail("Could not create diag timer: %s\n", strerror(errno));
            clock_is_realtime = 1;
        }
        diag_timer.it_value.tv_sec = diag_delay;
        diag_timer.it_value.tv_nsec = 0;
        diag_timer.it_interval.tv_sec = diag_delay;
        diag_timer.it_interval.tv_nsec = 0;
        if (timer_settime(diag_timerid, 0, &diag_timer, NULL))
            fail("Could not set diag timer: %s\n", strerror(errno));
    }
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

    diagt = diag_delay;
    init_time();
    init_diag();    /* ignore result: worst case we lose ^Z handling */
}

typedef struct {
    limitz_t val;
    uint ri;
    uint pi;
    uint rec;
} vcheck_t;
typedef struct {
    vcheck_t num;
    vcheck_t den;
} check_t;
check_t *checked = NULL;
uint cri, cpi;

void do_check(uint count) {
    for (uint fi = 0; fi < count; ++fi) {
        if (need_diag) {
            need_diag = 0;
            diag("%u %u %u/%u", cri, cpi, fi, count);
        }
        frag_t *fp = mmfrag_get(fi);
        for (uint vi = 0; vi < nv; ++vi) {
            check_t *cp = &checked[vi];
            range_t *rp = (range_t *)add_p(fp, sizeof(frag_t)
                    + vi * range_size());
            limit_t *lp = range_low(rp);
            if (abs(limitp_num(lp)) > cp->num.val) {
                cp->num.val = abs(limitp_num(lp));
                cp->num.ri = cri;
                cp->num.pi = cpi;
                cp->num.rec = fi;
            }
            if (abs(limitp_den(lp)) > cp->den.val) {
                cp->den.val = abs(limitp_den(lp));
                cp->den.ri = cri;
                cp->den.pi = cpi;
                cp->den.rec = fi;
            }
            lp = range_high(rp);
            if (abs(limitp_num(lp)) > cp->num.val) {
                cp->num.val = abs(limitp_num(lp));
                cp->num.ri = cri;
                cp->num.pi = cpi;
                cp->num.rec = fi;
            }
            if (abs(limitp_den(lp)) > cp->den.val) {
                cp->den.val = abs(limitp_den(lp));
                cp->den.ri = cri;
                cp->den.pi = cpi;
                cp->den.rec = fi;
            }
        }
    }
}

void check_all(void) {
    checked = calloc(nv, sizeof(check_t));
    mmfrag_init();
    for (cri = 0; cri < nresolve; ++cri) {
        cpi = resolves[cri].pi;
        uint count = mmfrag_open(0, cri, cpi);
        do_check(count);
        cpi = resolves[cri].pj;
        count = mmfrag_open(cri, cpi);
        do_check(count);
    }
    diag("");
    for (uint vi = 0; vi < nv; ++vi) {
        check_t *cp = &checked[vi];
        vcheck_t *vcp = &cp->num;
        printf("var %c: num %d at %u.%u.%u\n", 'a' + vi,
                vcp->val, vcp->ri, vcp->pi, vcp->rec);
        vcp = &cp->den;
        printf("var %c: den %d at %u.%u.%u\n", 'a' + vi,
                vcp->val, vcp->ri, vcp->pi, vcp->rec);
    }
    mmfrag_done();
    free(checked);
}

int main(int argc, char **argv, char **envp) {
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        if (strncmp("-Ls", arg, 3) == 0)
            diag_delay = strtoul(&arg[3], NULL, 10);
        else
            fail("unknown option '%s'", arg);
    }
    if (i + 2 == argc) {
        na = strtoul(argv[i++], NULL, 10);
        if (na < 1)
            fail("require a >= 1, not %lu", na);
        nb = strtoul(argv[i++], NULL, 10);
        if (nb < 1)
            fail("require b >= 1, not %lu", nb);
    } else
        fail("wrong number of arguments");

    init_post();

    check_all();

    done();
    return 0;
}
