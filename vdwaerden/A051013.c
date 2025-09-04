#include "cf.h"
#include "allsub.h"

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

/* A051013 [https://oeis.org/A051013] is the "number of nonaveraging subsets
 * on {1,2,...,n}", ie the number of subsets in which there are no three
 * elements in arithmetic progression.
 *
 * We calculate this using a helper function g(n), the number of such subsets
 * that include both 1 and n. Then A051013(n) = 1 + \sum_1^n{ (n+1-i)g(n) }.
 *
 * We calculate g(n) itself by storing all the subsets found, and then
 * iterating over the possible ways a previously found subset can be
 * included with an outer (1, n) pair.
 *
 * Note that this takes a lot of memory: I was able to calculate up to n=62
 * in 1080s using 50GB. Each subsequent step would take approximately
 * 1.35^(n-62) * 14GB, with a peak of 1.5x that requirement.
 *
 * In any case, as written this code can only calculate up to n=64. It would
 * be possible to extend to n=66 by avoiding storage for the last two values
 * of n, but beyond that we would need to go to 128-bit integers (which would
 * also double memory use).
 */

int main(int argc, char **argv) {
    init_allsub();
    for (uint n = 1; n <= MAX_ALLSUB; ++n) {
        find_allsub(n);
        sets_t f = 1;
        for (uint i = 1; i <= n; ++i)
            f += (sets_t)allsub[i].sets_count * (n + 1 - i);
        printf("%u %llu %llu (%.2fs)\n",
                n, (ullong)f, (ullong)allsub[n].sets_count, utime());
    }
    return 0;
}
