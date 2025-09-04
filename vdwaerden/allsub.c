#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cf.h"
#include "allsub.h"

#if MAX_ALLSUB > 32
#   define ONE 1ULL
#   define set_lsb lsb64
#else
#   define ONE 1U
#   define set_lsb lsb32
#endif

bool inited = false;
/* allsub[i] represents a resizable array of sets_t, representing the
 * subsets S_i of {2 .. n-1} for which S_i U { 1, n } is nonaveraging.
 * Each subset is stored as a bitvector in which bit 0 represents 2.
 */
allsub_t allsub[MAX_ALLSUB + 1];

/* pat[] lists values with two bits set; the first pat_count[i] of them
 * cover all such values with no bit set above 2^i.
 */
uint pat_count[MAX_ALLSUB + 1];
sets_t pat[(MAX_ALLSUB >> 1) + 1];

static inline void sets_realloc(allsub_t *asp, size_t new_size) {
    asp->sets = realloc(asp->sets, new_size * sizeof(sets_t));
    if (!asp->sets) {
        fprintf(stderr, "realloc failed resizing from %lu to %lu\n",
                asp->sets_size * sizeof(sets_t), new_size * sizeof(sets_t));
        exit(1);
    }
    asp->sets_size = new_size;
}

/* store a new subset, resizing the array as needed */
static inline void store(uint n, sets_t val) {
    allsub_t *asp = &allsub[n];
    if (asp->sets_count >= asp->sets_size)
        sets_realloc(asp, asp->sets_size ? (asp->sets_size * 3 / 2) : 256);
    asp->sets[asp->sets_count++] = val;
}

/* mandatory initialization */
void init_allsub(void) {
    if (!inited) {
        memset(allsub, 0, sizeof(allsub));
        store(0, 0);
        for (uint i = 0; i <= MAX_ALLSUB >> 1; ++i)
            pat[i] = (ONE << i) | (ONE << ((i << 1) + 1));
        for (uint i = 0; i <= MAX_ALLSUB; ++i)
            pat_count[i] = (i + 1) >> 1;
        inited = true;
    }
}

/* Find all non-balanced subsets of {1 .. n} that include both 1 and n,
 * storing the results in allsub[n] as bitvectors of the elements
 * excluding 1 and n.
 */
void find_allsub(uint n) {
    store(n, 0);
    if (n <= 2)
        return;

    /* When n is odd, there is an element halfway between 1 and n which
     * must always be disallowed.
     */
    sets_t midmask = (n & 1) ? (ONE << ((n - 3) >> 1)) : 0;

    /* Try building the interior of the subset using the known valid
     * subsets: we then need to exclude only those that interact with
     * with the additional 1 and n.
     */
    for (uint n2 = 1; n2 <= n - 2; ++n2) {
        /* this combines with the stored value to create the effective subset */
        sets_t incmask = ONE | (ONE << (n2 - 1));
        /* each value can be shifted into this many positions */
        uint shifts = n - 1 - n2;
        allsub_t *asp = &allsub[n2];
        for (size_t i = 0; i < asp->sets_count; ++i) {
            sets_t v0 = (asp->sets[i] << 1) | incmask;
            for (uint j = 0; j < shifts; ++j) {
                sets_t v = v0 << j;
                if (v & midmask)
                    goto exclude_v;
                /* For each bit set in v, we are excluded if it makes an AP
                 * with 1 and a later bit, or if it makes an AP with a later
                 * bit and n.
                 */
                sets_t w = v;
                while (w) {
                    uint bi = set_lsb(w);
                    uint bj = (bi << 1) + 1;
                    /* (1, bi, bj) would form an AP */
                    if (bj <= n - 2 && (w & (ONE << bj)))
                        goto exclude_v;
                    if ((bi & 1) == (n & 1)) {
                        bj = bi + ((n - 2 - bi) >> 1);
                        /* (bi, bj, n) would form an AP */
                        if (w & (ONE << bj))
                            goto exclude_v;
                    }
                    w ^= (ONE << bi);
                }
                store(n, v);
              exclude_v:
                ;
            }
        }
    }
    /* give back any spare memory */
    sets_realloc(&allsub[n], allsub[n].sets_count);
}
