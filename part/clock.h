#ifndef CLOCK_H
#define CLOCK_H

#include <sys/times.h>

void setup_clock(void);
void reset_clock(void);
void teardown_clock(void);

extern int clock_tick;
extern int gtime;
int curtime(void);
double difftime(clock_t t0, clock_t t1);

#define TIMETHIS(code) ({ \
	clock_t t0, t1; \
	t0 = curtime(); \
	({ code }); \
	t1 = curtime(); \
	difftime(t0, t1); \
})

#define GTIME difftime(gtime, curtime())

#endif
