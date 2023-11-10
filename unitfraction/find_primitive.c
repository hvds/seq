#include "unit.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/times.h>

uint limit;
uint count_prim = 0, count_qual = 0;
uint *cur_set;

long clock_tick;    /* ticks per second */

mpq_t *minq, *primq, best;

void init(void) {
    init_unit(limit + 1);
    minq = malloc((limit + 1) * sizeof(mpq_t));
    primq = malloc((limit + 1) * sizeof(mpq_t));
    for (uint i = 0; i <= limit; ++i) {
        QINIT(minq[i]);
        QINIT(primq[i]);
    }

    QINIT(best);
    cur_set = malloc(limit * sizeof(uint));

    mpq_t running_min, recip;
    QINIT(running_min);
    QINIT(recip);
    mpq_set_ui(running_min, 2, 1);
    for (uint i = limit; i > 0; --i) {
        mpq_set_ui(recip, 1, i);
        mpq_sub(running_min, running_min, recip);
        if (mpq_sgn(running_min) <= 0)
            break;
        mpq_set(minq[i], running_min);
    }
    QCLEAR(recip);
    QCLEAR(running_min);
}

void done(void) {
    for (uint i = 0; i <= limit; ++i) {
        QCLEAR(minq[i]);
        QCLEAR(primq[i]);
    }
    free(minq);
    free(primq);
    QCLEAR(best);
    free(cur_set);
    done_unit();
}

/* time used so far (CPU seconds) */
double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

#define diag(s) for (uint j = 0; j < c2; ++j) { if (j) printf(" "); printf("%u", cur_set[j]); } printf(": %s\n", s);

/* Recursively search for primitive S, q that may qualify. Called initially
 * with primq[0] = 0, min = 0, c = 0, we call recursively with primq[c] the
 * sum of reciprocals in the set so far, 'min' the largest element used so far,
 * 'c' the number of elements used so far.
 */
void try_prim(uint min, uint c) {
    const uint c2 = c + 1;
    for (uint n = min + 1; n <= limit; ++n) {
        /* set primq[c2] to be the new sum */
        mpq_set_ui(primq[c2], 1, n);
        mpq_add(primq[c2], primq[c2], primq[c]);
        if (mpq_cmp(primq[c2], minq[n]) < 0)
            break;

        /* record the current set, for diagnostics */
        cur_set[c] = n;

        /* if primq[c2] > 2, this is a primitive set: check if it satisfies
         * our constraints */
        int cmp = mpq_cmp_ui(primq[c2], 2, 1);
        if (cmp > 0) {
            /* if we add a "seen" check on primq[c2], must check in the order:
             * - continue if better_set
             * - continue if n < lim
             * - continue if seen
             * - mark as seen
             */
            /* constraint (a): max(S) = limit */
            if (n < limit)
                continue;
            /* constraint (c): S is optimal */
            if (better_set(primq[c2], c2))
                continue;

            /* show some progress */
            ++count_prim;
            if ((count_prim % 100) == 0) {
                gmp_printf("%lu %lu %Qu [", count_prim, count_qual, primq[c2]);
                for (uint i = 0; i <= c; ++i) {
                    if (i) printf(" ");
                    printf("%u", cur_set[i]);
                }
                printf("] (%.2fs)\n", timing());
            }

            /* constraint (d): S is flat */
            if (better_multi(primq[c2], c2))
                continue;

            /* solution found: record and report it */
            ++count_qual;
            gmp_printf("** found %Qu [", primq[c2]);
            for (uint i = 0; i <= c; ++i) {
                if (i) printf(" ");
                printf("%u", cur_set[i]);
            }
            printf("] (%.2fs)\n", timing());
            if (mpq_cmp(best, primq[c2]) < 0)
                mpq_set(best, primq[c2]);
        } else if (cmp < 0) {
            /* primq[c2] < 2, so check optimality, and recurse */
            if (better_set(primq[c2], c2))
                continue;
            try_prim(n, c2);
        } else {
            /* primq[c2] = 0, nothing to do */
            continue;
        }
    }
}

/* Given sets S and multisets M of positive integers,
 * find S and q = sum{1/s_i} such that:
 * a) max(S) = limit (specified)
 * b) S is primitive, i.e. 2 < q < 2 + 1/max(S)
 * c) S is optimal, i.e. no S' exists with |S'| < |S|, sum{1/s'_i} = q
 * d) S is flat, i.e. no M exists with |M| < |S|, sum{q/m_i} = q
 *
 * The intent is, among other things, to find max(q) with the above
 * constraints. We know S = {1, 2, 3, 7, 43, 47} qualifies with q > 2 + 1/49,
 * so for that purpose we need check only up to limit=48.
 *
 * Constraint (d) implies (c), but if S is optimal every subset of S must
 * be optimal; so we take advantage of that to prune the search space.
 *
 */
int main(int argc, char** argv) {
    uint arg = 1;

    while (arg < argc && argv[arg][0] == '-') {
        char* s = argv[arg++];
        if (strcmp(s, "--") == 0)
            break;
        fprintf(stderr, "Unknown option '%s'\n", s);
        argc = -1;  /* force usage message */
        break;
    }

    if (argc - arg != 1) {
        fprintf(stderr, "Usage: %s <limit>\n", argv[0]);
        return(-1);
    }
    limit = atoi(argv[arg++]);

    init();
    clock_tick = sysconf(_SC_CLK_TCK);
    setvbuf(stdout, NULL, _IOLBF, 0);

    mpq_set_ui(best, 2, 1);
    try_prim(0, 0);
    gmp_printf(
        "{\n  best => %Qu,\n  prim => %u,\n  qual => %u,\n  time => %.2fs\n}\n",
        best, count_prim, count_qual, timing()
    );
    done();
    return 0;
}
