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
int debug_suppress_write = 0;
uint sync_count;    /* force sync every n frag_write() calls */
bool sync_stderr = 0;

/* set to utime at start of run, minus last timestamp of recovery file */
double t0 = 0;
struct rusage rusage_buf;
extern inline double utime(void);

timer_t diag_timerid, log_timerid;
volatile bool need_diag, need_log;
bool clock_is_realtime = 0;
#define DIAG 1
#define LOG 600
double diag_delay = DIAG, log_delay = LOG, diagt, logt;
char *rpath = NULL;     /* path to log file */
FILE *rfp = NULL;       /* file handle to log file */
bool skip_recover = 0;  /* true if we should not attempt recovery */

void report(char *format, ...) {
    keep_diag();
    va_list ap;
    va_start(ap, format);
    gmp_vfprintf(stdout, format, ap);
    va_end(ap);
    if (rfp) {
        va_start(ap, format);
        gmp_vfprintf(rfp, format, ap);
        va_end(ap);
        fflush(rfp);
    }
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
    if (rfp)
        fclose(rfp);
    free(rpath);
}

void fail_silent(void) {
    /* we accept leaks on fatal error, but should close the log file */
    if (rfp)
        fclose(rfp);
    exit(0);
}
void fail(char *format, ...) {
    va_list ap;
    va_start(ap, format);
    gmp_vfprintf(stderr, format, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    /* we accept leaks on fatal error, but should close the log file */
    if (rfp)
        fclose(rfp);
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
    if (log_delay) {
        if (sigaction(SIGUSR2, &sa, NULL))
            fail("Could not set USR2 handler: %s\n", strerror(errno));
        sev.sigev_notify = SIGEV_SIGNAL;
        sev.sigev_signo = SIGUSR2;
        sev.sigev_value.sival_ptr = &log_timerid;
        if (timer_create(CLOCK_PROCESS_CPUTIME_ID, &sev, &log_timerid)) {
            /* guess that the CPUTIME clock is not supported */
            if (timer_create(CLOCK_REALTIME, &sev, &log_timerid))
                fail("Could not create log timer: %s\n", strerror(errno));
            clock_is_realtime = 1;
        }
        log_timer.it_value.tv_sec = log_delay;
        log_timer.it_value.tv_nsec = 0;
        log_timer.it_interval.tv_sec = log_delay;
        log_timer.it_interval.tv_nsec = 0;
        if (timer_settime(log_timerid, 0, &log_timer, NULL))
            fail("Could not set log timer: %s\n", strerror(errno));
    }
}

void init_pre(void) {
    t0 = utime();
}

void parse_checkpoint(char *s) {
    double dtime;
    int off = 0, extra;
    uint rid, pid;
    off_t fdi_off, fdor_off, fdop1_off, fdop2_off;

    if (s[0] == 'r') {
        if (EOF == sscanf(s, "resolve %u %lu %lu %lu %lu (%lfs)%n",
            &rid, &fdi_off, &fdor_off, &fdop1_off, &fdop2_off, &dtime, &extra)
        )
            fail("error parsing 305 line '%s'", s);
        set_resolve_cp(rid, fdi_off, fdor_off, fdop1_off, fdop2_off);
    } else {
        if (EOF == sscanf(s, "integrate %u %u %lu%n",
            &rid, &pid, &fdi_off, &extra)
        )
            fail("error parsing 305 line '%s'", s); 
        off += extra;
        for (uint pi = 0; pi < npaths; ++pi) {
            if (EOF == gmp_sscanf(&s[off], " %Qd%n", path_total[pi], &extra))
                fail("error parsing 305 line '%s'", s);
            off += extra;
        }
        if (EOF == sscanf(&s[off], " (%lfs)%n", &dtime, &extra))
            fail("error parsing 305 line '%s'", s);
        set_integrate_cp(rid, pid, fdi_off);
    }
    t0 -= dtime;
}

void recover(FILE *fp) {
    char *last_cp = NULL;
    char *curbuf = NULL;
    size_t len = 120, len_cp = 0;

    while (1) {
        ssize_t nread = getline(&curbuf, &len, fp);
        if (nread <= 0) {
            if (errno == 0)
                break;
            fail("error reading %s: %s", rpath, strerror(errno));
        }
        if (curbuf[nread - 1] != '\n'
                || memchr(curbuf, 0, nread) != NULL) {
            /* corrupt line, file should be truncated */
            off_t offset = ftello(fp);
            if (offset == -1)
                fail("could not ask offset: %s", strerror(errno));
            /* not ftruncate(), we are open only for reading */
            if (truncate(rpath, offset - nread) != 0)
                fail("could not truncate %s to %lu: %s", rpath, offset - nread,
                        strerror(errno));
            break;
        }
        if (strncmp("305 ", curbuf, 4) == 0) {
            char *t = last_cp;
            last_cp = curbuf;
            curbuf = t;
            size_t lt = len_cp;
            len_cp = len;
            len = lt;
        } else if (strncmp("001 ", curbuf, 4) == 0)
            ;   /* init */
        else if (strncmp("000 ", curbuf, 4) == 0)
            ;   /* comment */
        else
            fail("unexpected log line %.3s in %s", curbuf, rpath);
    }
    if (last_cp)
        parse_checkpoint(last_cp + 4);
    free(curbuf);
    free(last_cp);
}

void init_post(void) {
    nv = 2 * na * nb - na - nb;

    /* generate list of paths and resolutions
     * TODO: add strategies
     */
    init_paths(0);
    init_frags();
    init_calc(npaths);

    if (rpath) {
        printf("path %s\n", rpath);
        if (!skip_recover) {
            FILE *fp = fopen(rpath, "r");
            if (fp) {
                recover(fp);
                fclose(fp);
            }
        }

        rfp = fopen(rpath, "a");
        if (rfp == NULL)
            fail("%s: %s", rpath, strerror(errno));
        setlinebuf(rfp);
    }
    sync_stderr = isatty(STDERR_FILENO) ? 0 : 1;

    diagt = diag_delay;
    if (rfp)
        logt = log_delay;
    else
        logt = log_delay = 0;
    init_time();
    init_diag();    /* ignore result: worst case we lose ^Z handling */
}

int main(int argc, char **argv, char **envp) {
    int i = 1;
    prctl(PR_SET_NAME, argv[0]);
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        if (arg[1] == 'r') {
            rpath = (char *)malloc(strlen(&arg[2]) + 1);
            strcpy(rpath, &arg[2]);
        } else if (arg[1] == 'R')
            skip_recover = 1;
        else if (strncmp("-Ls", arg, 3) == 0)
            diag_delay = strtoul(&arg[3], NULL, 10);
        else if (strncmp("-Lf", arg, 3) == 0)
            log_delay = strtoul(&arg[3], NULL, 10);
        else if (strncmp("-ds", arg, 3) == 0)
            debug_split = 1;
        else if (strncmp("-di", arg, 3) == 0)
            debug_integrate = 1;
        else if (strncmp("-dw", arg, 3) == 0)
            debug_suppress_write = 1;
        else if (strncmp("-w", arg, 2) == 0)
            sync_count = strtoul(&arg[2], NULL, 10);
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

    split_all();
    for (uint pi = 0; pi < npaths; ++pi)
        integrate_path(pi);
    diag("");
    report_total();

    double tz = utime();
    report("int(%u, %u) = %Qd: (%.2fs)\n", na, nb, totalq, seconds(tz));
    done();
    return 0;
}
