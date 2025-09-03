#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "cf.h"

/* Could extend to 66 if we store actual subsets */
#define MAX_ALLSUB 64

typedef struct allsub_s {
    ullong *sets;
    size_t sets_size;
    size_t sets_count;
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

void sets_realloc(allsub_t *asp) {
    size_t new_size = asp->sets_size ? (asp->sets_size * 3 / 2) : 256;
    asp->sets = realloc(asp->sets, new_size * sizeof(ullong));
    asp->sets_size = new_size;
}

static inline void store(uint n, ullong val) {
    allsub_t *asp = &allsub[n];
    if (asp->sets_count >= asp->sets_size)
        sets_realloc(asp);
    asp->sets[asp->sets_count++] = val;
}

void init_allsub(void) {
    memset(allsub, 0, sizeof(allsub));
    store(0, 0);
    for (uint i = 0; i <= MAX_ALLSUB >> 1; ++i)
        pat[i] = (1ULL << i) | (1ULL << ((i << 1) + 1));
    for (uint i = 0; i <= MAX_ALLSUB; ++i)
        pat_count[i] = (i + 1) >> 1;
    inited = true;
}

/* Find all non-balanced subsets of {1 .. n} that include both 1 and n.
 * Supports 1 <= n <= 64, if there is enough memory; could in principle
 * be extended to n <= 66.
 */
void find_allsub(uint n) {
    store(n, 0);
    if (n <= 2)
        return;
    ullong midmask = (n & 1) ? (1ULL << ((n - 3) >> 1)) : 0;
    for (uint n2 = 1; n2 <= n - 2; ++n2) {
        uint rot = n - 1 - n2;
        ullong incmask = 1ULL | (1ULL << (n2 - 1));
        allsub_t *asp = &allsub[n2];
        for (size_t i = 0; i < asp->sets_count; ++i) {
            ullong v0 = (asp->sets[i] << 1) | incmask;
            for (uint j = 0; j < rot; ++j) {
                ullong v = v0 << j;
                if (v & midmask)
                    goto fail_v;
                ullong w = v;
                while (w) {
                    uint bi = lsb64(w);
                    uint bj = (bi << 1) + 1;
                    if (bj <= n - 2 && (w & (1ULL << bj)))
                        goto fail_v;
                    if (((bi ^ n) & 1) == 0) {
                        bj = bi + ((n - 2 - bi) >> 1);
                        if (w & (1ULL << bj))
                            goto fail_v;
                    }
                    w ^= (1ULL << bi);
                }
                store(n, v);
              fail_v:
                ;
            }
        }
    }
}

int main(int argc, char **argv) {
    init_allsub();
    for (uint n = 1; n <= MAX_ALLSUB; ++n) {
        find_allsub(n);
        ullong f = 1;
        for (uint i = 1; i <= n; ++i)
            f += (ullong)allsub[i].sets_count * (n + 1 - i);
        printf("%u %llu %llu (%.2fs)\n",
                n, f, (ullong)allsub[n].sets_count, utime());
    }
    return 0;
}
