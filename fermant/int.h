#ifndef INT_H
#define INT_H

#include <gmp.h>
#include "types.h"

extern uint na, nb; /* we are calculating g(na, nb) */
extern uint nv;     /* ... which is an integration over nv variables */
extern int debug_split;
extern int debug_integrate;

double elapsed(void);
double seconds(double t1);
void fail(char *format, ...);
void fail_silent(void);
void report(char *format, ...);

#endif
