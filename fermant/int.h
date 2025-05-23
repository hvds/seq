#ifndef INT_H
#define INT_H

#include <gmp.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "types.h"

extern uint na, nb; /* we are calculating g(na, nb) */
extern uint nv;     /* ... which is an integration over nv variables */
extern int debug_split;
extern int debug_integrate;
extern int debug_suppress_write;
extern uint sync_count;
extern bool sync_stderr;
extern volatile bool need_diag, need_log;
extern FILE *rfp;
extern struct rusage rusage_buf;

double elapsed(void);
double seconds(double t1);
void fail(char *format, ...);
void fail_silent(void);
void report(char *format, ...);

__inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

#endif
