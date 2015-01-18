#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <errno.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>
#include "pp.h"

long clock_tick;

char* use_str =
	"Usage: %s n [kstart [kend]]\n"
	"Search for a(n) from k = (kstart or 1) to k = (kend or kstart or \\inf)\n";

void usage(char* prog) {
	fprintf(stderr, use_str, prog);
	exit(0);
}

void diag_signal(int signo) {
	diag_signal_seen = 1;
}

void setup_signals(void) {
	struct sigaction sa;

	diag_signal_seen = 0;
	sa.sa_handler = &diag_signal;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART;
	if (sigaction(SIGUSR1, &sa, (struct sigaction*)NULL) != 0) {
		fprintf(stderr, "Attempt to set signal handler for SIGUSR1: %s\n", strerror(errno));
		exit(1);
	}
	printf("Diagnostics available with 'kill -USR1 %d'\n", (int)getpid());
}

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

int main(int argc, char** argv) {
	int n = 0, kstart = 0, kend;
	int i, j, success;

	if (argc > 1)
		n = atoi(argv[1]);
	else
		usage(argv[0]);
	if (argc > 2)
		kstart = atoi(argv[2]);
	if (argc > 3)
		kend = atoi(argv[3]);
	else
		kend = kstart;
	setup_signals();
    clock_tick = sysconf(_SC_CLK_TCK);
	for (i = n ? n : 1; n ? (i <= n) : 1; ++i) {
		for (j = kstart ? kstart : 1; kend ? (j <= kend) : 1; ++j) {
			/* printf("Try n=%d, k=%d\n", i, j); */
			setup_pp(j);
			success = pp_find(i);
			teardown_pp();
			if (success)
				break;
		}
	}
	return 0;
}
