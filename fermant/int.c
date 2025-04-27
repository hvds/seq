/* needed on FreeBSD for getline() to be exported from stdio.h */
#define _WITH_GETLINE

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
/* used to set process title */
#include <sys/prctl.h>
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

/* set to utime at start of run, minus last timestamp of recovery file */
double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}
timer_t diag_timerid;
volatile bool need_diag;
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
    need_diag = 1;
}

void init_time(void) {
    struct sigaction sa;
    struct sigevent sev;
    struct itimerspec diag_timer;

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

    diagt = diag_delay;
    init_time();
    init_diag();    /* ignore result: worst case we lose ^Z handling */

    /* generate list of paths and resolutions
     * TODO: add strategies
     */
    init_paths(0);
    init_frags();
    init_calc(npaths);

}

void run(uint recover) {
    uint num = split_all(recover);
    int fdi = resolve_reader(nresolve);
    uint count = 0;
    while (read_frag(fdi)) {
        if (need_diag) {
            diag("int %u/%u", count, num);
            need_diag = 0;
        }
        ++count;
        integrate(nfrags - 1);
        reset_frags();
    }
    if (fdi)
        close(fdi);
    diag("");
    report_total();
}

typedef enum {
    IS_DEEPER = 0,
    IS_NEXTX,
    IS_NEXT,
    IS_MIDP
} e_is;
/* On recovery, set up the recursion stack to the point we had reached.
 * Returns IS_DEEPER if we should continue by recursing deeper from this
 * point; returns IS_NEXTX if we should continue by advancing the power
 * applied at the current position; returns IS_NEXT if we should continue
 * by advancing the current level; and returns IS_MIDP if we should continue
 * via walk_midp().
 */
e_is insert_stack(void) {
    e_is jump = IS_DEEPER;
    /* ... */
    return jump;
}

/* we emulate recursive calls via the levels[] array */
void recurse(e_is jump_continue) {
    /* ... */
}
int main(int argc, char **argv, char **envp) {
    int i = 1;
    uint run_recover = 0;
    prctl(PR_SET_NAME, argv[0]);
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        if (strncmp("-r", arg, 2) == 0)
            run_recover = strtoul(&arg[2], NULL, 10);
        else if (strncmp("-Ls", arg, 3) == 0)
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

    run(run_recover);
    keep_diag();

    double tz = utime();
    report("int(%u, %u) = %Qd: (%.2fs)\n", na, nb, totalq, seconds(tz));
    done();
    return 0;
}
