/* See http://www.pixelbeat.org/programming/gcc/static_assert.html if this
   doesn't give you static_assert(). */
#define __USE_ISOC11
#include <assert.h>

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;
typedef unsigned char uchar;

/* we may search for f(n): n <= MAXN */
#define MAXN 1022
typedef ushort n_t;
static_assert(MAXN < (1 << (sizeof(n_t) * 8)), "n_t too small for MAXN");

typedef uint scratch_t;
static_assert(MAXN <= (UINT_MAX >> 1), "scratch_t too small for 2 * MAXN");

/* we may search for f(n): f(n) <= MAXF */
#define MAXF 255
typedef uchar f_t;
static_assert(MAXF < (1 << (sizeof(f_t) * 8)), "f_t too small for MAXF");

n_t n;              /* current target n */
f_t max;            /* cardinality of largest subset found so far */
                    /* NB f(n) must be either max or max+1 */
f_t f3[MAXN + 1];   /* known values of f() */
n_t s[MAXF + 1];    /* the integers in the current subset */
f_t s_i;            /* offset in s[] being worked on */
ulong iter;         /* number of iterations, for debugging/stats only */

double t0 = 0;      /* base for timing calculations */
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

/* block_t is a bit vector representing n_t blocked by current subset; we
 * need to allow room up to MAXN + 1. */
/* TODO: there may be value in making the size adaptive - it can be smaller
 * for lower n, so we can move fewer bytes around.
 */
#define SIZEB (((MAXN + 1) + 7) >> 3)
typedef struct block_s {
    uchar b[SIZEB];
} block_t;
block_t blocks[MAXF + 1];

static inline void copy_block(block_t *dest, block_t *src) {
    memcpy(dest, src, sizeof(block_t));
}

static inline bool is_blocked(block_t *bp, n_t val) {
    n_t off = val >> 3;
    n_t bit = 1 << (val & 7);
    return (bp->b[off] & bit) ? true : false;
}

static inline void set_blocked(block_t *bp, n_t val) {
    n_t off = val >> 3;
    n_t bit = 1 << (val & 7);
    bp->b[off] |= bit;
}

void reset_data(void) {
    s[0] = 1;   /* always fixed */
    s[1] = 1;   /* ready to increment to 2 */
    bzero(&blocks[0], sizeof(block_t) * 2);
    s_i = 1;
}

void disp_result(bool ok) {
    printf("%u (%lu %.2fs): ", n, iter, utime() - t0);
    if (ok) {
        printf("[");
        for (f_t s_j = 0; s_j < s_i; ++s_j) {
            if (s_j)
                printf(" ");
            printf("%u", s[s_j]);
        }
        printf("]");
    } else {
        printf("no improvement");
    }
    printf("\n");
}

#define ULBITS (8 * sizeof(unsigned long))

static inline n_t bits_avail(block_t *bp, n_t after) {
    if (after >= n)
        return 0;
    n_t sum = 1;    /* n is always allowed */
    for (n_t i = (after + 1) & ~(ULBITS - 1); i < n; i += ULBITS) {
        ulong word = *(ulong *)(&bp->b[i >> 3]);
        if (i <= after)
            word |= ((1UL << (after + 1 - i)) - 1);
        if (i + ULBITS > n)
            word |= ~((1UL << (n - i)) - 1);
        sum += ULBITS - __builtin_popcountl(word);
    }
    return sum;
}

void findmax(void) {
    n = 2;
    f3[0] = 0;
    f3[1] = 1;
    max = 1;
    iter = 0;
    reset_data();
  TRY_NEXT:
    while (s_i > 0) {
        ++iter;
        n_t cur = ++s[s_i];
        if (cur > n || s_i + f3[n - cur + 1] <= max)
            goto derecurse;

        block_t *bprev = &blocks[s_i - 1];
        if (is_blocked(bprev, cur))
            continue;

        /* ok to continue with this subset */

        /* block the third element of any new 2-element APs */
        block_t *bcur = &blocks[s_i];
        copy_block(bcur, bprev);
        scratch_t cur2 = (scratch_t)cur + cur;
        bool fail = false;
        for (f_t s_j = s_i; s_j > 0;) {
            --s_j;
            scratch_t targ = cur2 - s[s_j];
            if (targ >= (scratch_t)n) {
                if (targ == (scratch_t)n) {
                    /* if subset cannot contain n, it cannot be a new max */
                    fail = true;
                } else if (targ == (scratch_t)n + 1) {
                    /* mark this (in case of continuation) */
                    set_blocked(bcur, (n_t)targ);
                }
                /* any further values are beyond bounds */
                break;
            }
            set_blocked(bcur, (n_t)targ);
        }
        if (fail)
            continue;
        if (s_i + 1 + bits_avail(bcur, cur) <= max)
            continue;

        /* advance for next element */
        ++s_i;
        s[s_i] = cur;   /* ready to increment */

        if (s_i > max) {
            /* we have a solution: f(n) = f(n-1) + 1 */
            disp_result(1);
            ++max;
            f3[n] = max;
            ++n;
            iter = 0;
            /* preserve invariant "n is never blocked" */
            while (s_i > 0 && is_blocked(&blocks[s_i - 1], n))
                --s_i;
        }
        continue;
      derecurse:
        --s_i;
    }
    /* no solution: f(n) = f(n-1) */
    disp_result(0);
    reset_data();
    f3[n] = max;
    ++n;
    iter = 0;
    goto TRY_NEXT;
}

int main(int argc, char **argv) {
    findmax();
    return 0;
}
