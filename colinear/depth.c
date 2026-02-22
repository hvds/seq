#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

#include "poly.h"
#include "../inc/diag.h"

/* arbitrary limit; we expect to need less than 120 for n=5 */
#define MAXDEPTH 150

uint opt_n;
typedef struct s_level {
    t_shape *s;
    t_iter *i;
    uint cur;
    uint max;
} t_level;
t_level *levels = NULL;
uint maxdepth = 0;
uint level;

double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

double diag_delay = 1, log_delay = 0, diagt, logt;
timer_t diag_timerid, log_timerid;
volatile bool need_work, need_diag, need_log;

double seconds(double t1) {
    return (t1 - t0);
}

double elapsed(void) {
    return seconds(utime());
}

void fail(char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vprintf(format, ap);
    printf("\n");
    va_end(ap);
    exit(1);
}

char buf[4096];
void diag_prog(uint level) {
    size_t off = 0;
    for (uint i = 0; i < level; ++i) {
        if (i)
            off += snprintf(&buf[off], sizeof(buf) - off, " ");
        off += snprintf(&buf[off], sizeof(buf) - off, "%u",
                levels[i].max - levels[i].cur);
    }
    diag(buf);
    need_diag = 0;
}

void handle_sig(int sig) {
    need_work = 1;
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
        }
        log_timer.it_value.tv_sec = log_delay;
        log_timer.it_value.tv_nsec = 0;
        log_timer.it_interval.tv_sec = log_delay;
        log_timer.it_interval.tv_nsec = 0;
        if (timer_settime(log_timerid, 0, &log_timer, NULL))
            fail("Could not set log timer: %s\n", strerror(errno));
    }
}

void resize_levels(uint depth) {
    levels = realloc(levels, depth * sizeof(t_level));
    for (uint i = maxdepth; i < depth; ++i) {
        levels[i].s = shape_new((t_dim){ MAXDIM, MAXDIM });
        levels[i].i = malloc(sizeof(t_iter) + max_poly_size());
    }
    maxdepth = depth;
}

void run(void) {
    shape_empty(levels[0].s);
    level = 0;
    uint max_level = 0;
    t_level *cur, *next;
    goto get_iter;

    while (1) {
      next_neighbour:
        if (need_diag)
            diag_prog(level);
        cur = &levels[level];
        t_point p = shape_iter(cur->i);
        if (p.x == 0 && p.y == 0) {
            /* derecurse */
            if (level == 0)
                break;
            --level;
            continue;
        }
        ++cur->cur;
        shape_mark_disallowed(cur->s, p);
        if (level + 1 >= maxdepth) {
            resize_levels(maxdepth * 3 / 2);
            cur = &levels[level];
        }
        ++level;
        next = &levels[level];
        shape_append(next->s, cur->s, p);
        if (level > max_level) {
            keep_diag();
            fprintf(stderr, "size %u", level);
            diag_shape(next->s, ":");
            max_level = level;
        }

      get_iter:
        cur = &levels[level];
        shape_neighbours(cur->i, cur->s);
        cur->max = 0;
        while (1) {
            t_point p = shape_iter(cur->i);
            if (p.x == 0 && p.y == 0)
                break;
            if (test_colinear(cur->s, p, opt_n)) {
                /* ok */
                ++cur->max;
            } else {
                shape_mark_disallowed(cur->s, p);
                iter_remove(cur->i, p);
            }
        }
        iter_reset(cur->i);
        cur->cur = 0;
    }
}

void init(void) {
    init_diag();
    init_time();
    init_poly();
    resize_levels(100);
}

void done(void) {
    for (uint i = 0; i < maxdepth; ++i) {
        shape_free(levels[i].s);
        free(levels[i].i);
    }
    done_poly();
}

int main(int argc, char **argv) {
    uint argi = 1;
    while (argi < argc && *argv[argi] == '-') {
        char *arg = argv[argi++];
        if (strcmp(arg, "--") == 0)
            break;
        fprintf(stderr, "Unknown option '%s'\n", arg);
        exit(1);
    }
    if (argi + 1 != argc) {
        fprintf(stderr, "Usage: %s <options> <n>\n", argv[0]);
        exit(1);
    }
    opt_n = strtoul(argv[argi], 0, 10);
    if (opt_n < 2 || opt_n > 8) {
        fprintf(stderr, "Expected 2 <= n <= 8, not %u\n", opt_n);
        exit(1);
    }
    init();
    run();
    keep_diag();
    done();
    return 0;
}
