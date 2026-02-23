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
uint maxsize = 0;
uint *rstack = NULL;
uint rstacksize = 0;

double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

char *rpath = NULL; /* path to log file */
FILE *rfp = NULL;   /* file handle to log file */
double diag_delay = 1, log_delay = 600;
timer_t diag_timerid, log_timerid;
volatile bool need_work, need_diag, need_log;

double seconds(double t1) {
    return (t1 - t0);
}

double elapsed(void) {
    return seconds(utime());
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
    vprintf(format, ap);
    printf("\n");
    va_end(ap);
    /* we accept leaks on fatal error, but should close the log file */
    if (rfp)
        fclose(rfp);
    exit(1);
}

void report(char *format, ...) {
    keep_diag();
    va_list ap;
    va_start(ap, format);
    vfprintf(stdout, format, ap);
    va_end(ap);

    if (rfp) {
        va_start(ap, format);
        vfprintf(rfp, format, ap);
        va_end(ap);
        fflush(rfp);
    }
}

char buf[4096];
void diag_prog(uint level) {
    double t1 = utime();
    size_t off = 0;
    for (uint i = 0; i < level; ++i) {
        if (i)
            off += snprintf(&buf[off], sizeof(buf) - off, " ");
        off += snprintf(&buf[off], sizeof(buf) - off, "%u",
                levels[i].max - levels[i].cur);
    }
    if (need_diag) {
        diag("%s", buf);
        need_diag = 0;
    }
    if (rfp && need_log) {
        fprintf(rfp, "305 %s (%.2fs)\n", buf, seconds(t1));
        need_log = 0;
    }
    need_work = 0;
}

void parse_305(char *s) {
    uint i = 0;
    uint stacksize = 0;
    while (1) {
        if (*s == '(')
            break;
        if (i >= stacksize) {
            stacksize += 100;
            rstack = realloc(rstack, stacksize * sizeof(uint));
        }
        int off = 0;
        sscanf(s, "%u %n", &rstack[i++], &off);
        if (off == 0)
            fail("501 error parsing 305 line '%s'", s);
        s += off;
    }
    rstacksize = i;
    int off = 0;
    double dtime;
    sscanf(s, "(%lfs)%n", &dtime, &off);
    if (off == 0)
        fail("501 could not parse 305 line '%s'", s);
    t0 -= dtime;
}

void report_best(uint size, t_shape *shape, bool running) {
    double t1 = utime();
    if (running) {
        need_diag = 1;
        need_log = 1;
        diag_prog(size);
    }
    keep_diag();
    printf("size %u (%.2f)", size, seconds(t1));
    diag_shape(shape, ":");
    maxsize = size;

    if (rfp) {
        size_t rz = report_size(shape->d);
        char buf[rz];
        report_shape(buf, rz, shape);
        fprintf(rfp, "202 %u:%s (%.2fs)\n", size, buf, seconds(t1));
    }
}

void parse_202(char *s, uint *best, t_shape **shape) {
    uint size;
    int off = 0;
    sscanf(s, "%u:%n", &size, &off);
    if (off == 0)
        fail("501 error parsing 202 line '%s'", s);
    if (size <= *best)
        return;
    *best = size;
    if (*shape == NULL)
        *shape = shape_new((t_dim){ MAXDIM, MAXDIM });
    t_shape *sh = *shape;
    s += off;
    off = 0;
    sscanf(s, "%u:%u:%n", &sh->d.x, &sh->d.y, &off);
    if (off == 0)
        fail("501 error parsing 202 line '%s'", s);
    s += off - 1;
    shape_reset(sh);
    t_point p;
    for (p.x = 0; p.x <= sh->d.x + 1; ++p.x) {
        ++s;
        for (p.y = 0; p.y <= sh->d.y + 1; ++p.y) {
            char c = *s++;
            switch (c) {
              case '*':
                shape_mark_point(sh, p);
                /* fall through */
              case 'o':
                shape_mark_disallowed(sh, p);
                /* fall through */
              case '.':
                break;
              default:
                fail("501 error parsing 202 line '%s'", s);
            }
        }
    }
    /* ignore tail " (%.2fs)\n" */
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

void recover(FILE *fp) {
    char *curbuf = NULL, *last305 = NULL;
    size_t len = 120, len305 = 0;
    uint best = 0;
    t_shape *shape = NULL;

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
            char *t = last305;
            last305 = curbuf;
            curbuf = t;
            size_t lt = len305;
            len305 = len;
            len = lt;
        } else if (strncmp("202 ", curbuf, 4) == 0)
            parse_202(curbuf + 4, &best, &shape);
        else if (strncmp("000 ", curbuf, 4) == 0)
            ;   /* comment */
        else
            fail("unexpected log line %.3s in %s", curbuf, rpath);
    }
    /* parse this first, to get t0 */
    if (last305)
        parse_305(last305 + 4);
    if (best) {
        report_best(best, shape, 0);
        shape_free(shape);
    }
    free(curbuf);
    free(last305);
}

void resize_levels(uint depth) {
    levels = realloc(levels, depth * sizeof(t_level));
    for (uint i = maxdepth; i < depth; ++i) {
        levels[i].s = shape_new((t_dim){ MAXDIM, MAXDIM });
        levels[i].i = malloc(sizeof(t_iter) + max_poly_size());
    }
    maxdepth = depth;
}

void init_neighbours(uint l) {
    t_level *cur = &levels[l];
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

void run(void) {
    shape_empty(levels[0].s);
    level = 0;
    if (rstack) {
        if (rstacksize + 1 >= maxdepth)
            resize_levels(rstacksize * 3 / 2);
        for (uint i = 0; i < rstacksize; ++i) {
            init_neighbours(i);
            t_level *cur = &levels[i];
            while (cur->cur + rstack[i] + 1 < cur->max) {
                t_point p = shape_iter(cur->i);
                shape_mark_disallowed(cur->s, p);
                ++cur->cur;
            }
            if (i + 1 < rstacksize) {
                t_point p = shape_iter(cur->i);
                shape_mark_disallowed(cur->s, p);
                shape_append(levels[i + 1].s, cur->s, p);
                ++cur->cur;
            }
        }
        level = rstacksize;
    }

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
        if (level > maxsize)
            report_best(level, next->s, 1);

      get_iter:
        init_neighbours(level);
    }
}

void init(void) {
    t0 = utime();
    init_diag();
    init_time();
    init_poly();
    resize_levels(100);
    if (rpath) {
        printf("path %s\n", rpath);
        FILE *fp = fopen(rpath, "r");
        if (fp) {
            recover(fp);
            fclose(fp);
        }
        rfp = fopen(rpath, "a");
        if (rfp == NULL)
            fail("%s: %s", rpath, strerror(errno));
        setlinebuf(rfp);
    }
}

void done(void) {
    for (uint i = 0; i < maxdepth; ++i) {
        shape_free(levels[i].s);
        free(levels[i].i);
    }
    done_poly();
    if (rfp)
        fclose(rfp);
    free(rpath);
}

int main(int argc, char **argv) {
    uint argi = 1;
    while (argi < argc && *argv[argi] == '-') {
        char *arg = argv[argi++];
        if (strcmp(arg, "--") == 0)
            break;
        if (arg[1] == 'r') {
            rpath = malloc(strlen(&arg[2]) + 1);
            strcpy(rpath, &arg[2]);
        } else {
            fprintf(stderr, "Unknown option '%s'\n", arg);
            exit(1);
        }
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
