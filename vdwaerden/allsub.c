#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "cf.h"

/* Could extend to 66 if we store actual subsets */
#define MAX_ALLSUB 64

typedef struct allsub_s {
    /* just store counts for now */
    uint count;
} allsub_t;

bool inited = false;
allsub_t allsub[MAX_ALLSUB + 1];
uint pat_count[MAX_ALLSUB + 1];
ullong pat[(MAX_ALLSUB >> 1) + 1];

double t0 = 0;      /* base for timing calculations */
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

void init_allsub(void) {
    memset(allsub, 0, sizeof(allsub));
    allsub[0].count = 1;
    for (uint i = 0; i <= MAX_ALLSUB >> 1; ++i)
        pat[i] = (1ULL << i) | (1ULL << ((i << 1) + 1));
    for (uint i = 0; i <= MAX_ALLSUB; ++i)
        pat_count[i] = (i + 1) >> 1;
    inited = true;
}

static inline void store(uint n, ullong val) {
    ++allsub[n].count;
}

/* Find all non-balanced subsets of {1 .. n} that include both 1 and n.
 * Supports 1 <= n <= 64, if there is enough memory; could in principle
 * be extended to n <= 66.
 */
void find_allsub(uint n) {
    if (!inited)
        init_allsub();
    ullong max = (n < 2) ? 0 : ((1ULL << (n - 2)) - 1);
    ullong mask = 1ULL | (1ULL << (n - 1));
    for (ullong i = 0; i <= max; ++i) {
        ullong val = (i << 1) | mask;
        while (val) {
            uint bi = lsb64(val);
            val >>= bi + 1;
            if (val == 0)
                break;
            uint bj = msb64(val);
            uint count = pat_count[bj];
            for (uint ji = 0; ji < count; ++ji)
                if ((val & pat[ji]) == pat[ji])
                    goto fail;
        }
        store(n, i);
      fail:
        ;
    }
}

int main(int argc, char **argv) {
    for (uint n = 1; n <= MAX_ALLSUB; ++n) {
        find_allsub(n);
        uint f = 1;
        for (uint i = 1; i <= n; ++i)
            f += allsub[i].count * (n + 1 - i);
        printf("%u %u %u (%.2fs)\n", n, f, allsub[n].count, utime());
    }
    return 0;
}
