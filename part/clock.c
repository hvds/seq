#include "clock.h"
#include <unistd.h>

/* clock functions */
int clk_tck;
int gtime;

void setup_clock(void) {
	clk_tck = sysconf(_SC_CLK_TCK);
	gtime = curtime();
}

void reset_clock(void) {
	gtime = curtime();
}

void teardown_clock(void) {
	/* nothing to do */
}

double difftime(clock_t t0, clock_t t1) {
	return ((double)t1 - t0) / clk_tck;
}

int curtime(void) {
	struct tms t;
	times(&t);
	return (int) t.tms_utime;
}
