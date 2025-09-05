/* See http://www.pixelbeat.org/programming/gcc/static_assert.html if this
   doesn't give you static_assert(). */
#define __USE_ISOC11
#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "cf.h"

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
n_t debug_n = 0;    /* show debug info when n == debug_n */
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

/* Given a block_t telling us which values are available from cur+1 to n,
 * we want to know the cardinality of the largest 3AP-free subset of those
 * values that includes n. We will make a 16-bit lookup table of two values:
 * one that requires the inclusion of the (virtual) msb, and one that gives
 * a free choice. We can then stitch together a total using the latter for
 * the last 17 bits, and the former for any remaining segments.
 */
#define LOOKUP_BYTES 3
#define LOOKUP_BITS (LOOKUP_BYTES << 3)
#define LOOKUP_SIZE (1 << LOOKUP_BITS)
#define LOOKUP_MASK ((1 << LOOKUP_BITS) - 1)
f_t lookup_sub[LOOKUP_SIZE];        /* free choice */
f_t lookup_sub_force[LOOKUP_SIZE];  /* (LOOKUP_BITS+1)th bit forced */

bool better_sub(uint avail, uint have, f_t need) {
    while (1) {
        if (lookup_sub[avail] < need)
            return false;
        uint bi = msb32(avail);
        uint b = 1 << bi;
        uint have2 = have | b;
        uint avail2 = avail ^ b;
        f_t need2 = need - 1;
        for (uint j = 1; bi + j <= UINT32_HIGHBIT && bi >= j; ++j)
            if (have2 & (1 << (bi + j)))
                avail2 &= ~(1 << (bi - j));
        if (need2 <= 1) {
            if (avail2 >= need2)
                return true;
            goto skip2;
        }
        if (avail2 == 0)
            goto skip2;
        if (better_sub(avail2, have2, need2))
            return true;
      skip2:
        avail ^= b;
    }
}
    
void init_lookup_sub(void) {
    lookup_sub[0] = 0;
    lookup_sub_force[0] = 1;
    for (uint i = 1; i < LOOKUP_SIZE; ++i) {
        /* free choice */
        int bits = __builtin_popcount(i);
        if (bits <= 2) {
            lookup_sub[i] = bits;
        } else if ((i & 1) == 0) {
            lookup_sub[i] = lookup_sub[i >> 1];
        } else {
            uint b = 1 << msb32(i);
            uint j = i ^ b;
            f_t prev = lookup_sub[j];
            lookup_sub[i] = better_sub(j, b, prev)
                    ? prev + 1 : prev;
        }

        /* force 1 << LOOKUP_BITS */
        if (bits <= 1) {
            lookup_sub_force[i] = bits + 1;
        } else {
            uint b = 1 << msb32(i);
            uint j = i ^ b;
            f_t prev = lookup_sub_force[j];
            lookup_sub_force[i] = better_sub(i, 1 << LOOKUP_BITS, prev)
                    ? prev + 1 : prev;
        }
    }

    /* our actual data has bits set for _disallowed_ values, so invert */
    for (uint i = 0; i < LOOKUP_SIZE; ++i) {
        uint j = (~i) & LOOKUP_MASK;
        if (i < j) {
            f_t temp = lookup_sub[i];
            lookup_sub[i] = lookup_sub[j];
            lookup_sub[j] = temp;
            temp = lookup_sub_force[i];
            lookup_sub_force[i] = lookup_sub_force[j];
            lookup_sub_force[j] = temp;
        }
    }
}

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

/* Extract the specified range of LOOKUP_BITS bits from the given pointer.
 * The result is shifted such that the end bit is bit (LOOKUP_BITS-1) of
 * the result. Bits with negative indices are permitted, and set to 1.
 */
static inline uint extract_LUB(uchar *cp, signed int start_bit) {
    uint word = 0;
    for (uint i = 0; i <= LOOKUP_BYTES; ++i) {
        signed int target = start_bit + (i << 3);
        if (target >= 0)
            word |= cp[target >> 3] << (i << 3);
    }
    word >>= ((uint)start_bit) & 7;
    if (start_bit < 0)
        word |= (1 << -start_bit) - 1;
    return word & LOOKUP_MASK;
}

/* eg after = 12, n = 57 we want:
 *   w1 = (41..47 >> 25) | (48..56 << 7)  (force)
 *   w2 = (25..31 >> 25) | (32..40 << 7)
 *   w3 =                  (13..24 << 7)
 * eg after = 12, n = 26 we want:
 *   w1 =                  (13..25 << 6)  (force)
 */
static inline n_t bits_avail(block_t *bp, n_t after) {
    signed int from = (signed int)after + 1;
    if (from >= (signed int)n)
        return (from == (signed int)n) ? 1 : 0;
    n_t sum = 0;
    bool force = true;
    for (signed int i = (signed int)n; i > from; ) {
        i -= LOOKUP_BITS;
        uint word = extract_LUB(&bp->b[0], i);
        if (i < from)
            word |= (1 << (from - i)) - 1;
        sum += (force) ? lookup_sub_force[word] : lookup_sub[word];
        force = false;
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
        if (s_i + 1 + bits_avail(bcur, cur) <= max) {
#ifdef DEBUG
            if (n == debug_n) {
                printf("[");
                for (uint dsi = 0; dsi <= s_i; ++dsi) {
                    if (dsi)
                        printf(" ");
                    printf("%d", s[dsi]);
                }
                printf("] avail = %d\n", bits_avail(bcur, cur));
            }
#endif
            continue;
        }

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
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (arg[1] == '-')
            break;
        else if (arg[1] == 'd') {
            debug_n = (n_t)strtoul(&arg[2], NULL, 10);
        } else {
            printf("unknown option '%s'", arg);
            exit(1);
        }
    }
    setlinebuf(stdout);
    init_lookup_sub();
    t0 = utime();
    findmax();
    return 0;
}
