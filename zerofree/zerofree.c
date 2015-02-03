#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>
#include <gmp.h>
#include "array.h"
#include "seen.h"

#define DEBUG 0

long clock_tick;

uint base = 10;     /* expand expressions in this base */
ulong total = 0;    /* number of distinct values seen */
array_t curpend;    /* pending array for the current generation */
array_t newpend;    /* pending array for the next generation */
uint pendoff;       /* offset into current generation */
mpz_t expanded, expand_new; /* big integers used in expand() */

/* default calculation: n -> n^2 */
#define CALCULATE() mpz_mul(expanded, n, n)

#if 0
/* alternative calculation: n -> 2n */
#   define CALCULATE() mpz_mul_ui(expanded, n, 2)
#endif

double timing(void) {
    struct tms ttd;
    times(&ttd);
    return ((double)ttd.tms_utime) / clock_tick;
}

void dumpstats(void) {
    size_t i;
    printf("pend %u; ", newpend.count);
    dump_seen();
    printf(", total: %lu [%.2fs]\n", total, timing());
}

void final(void) {
    mpz_t z;

    printf("****** finished ******\n");
    dumpstats();
    mpz_init(z);
    max_seen(z);
    gmp_printf("last value: %Zd\n", z);
    mpz_clear(z);
    free_array(&newpend);
    free_array(&curpend);
    free_seen();
    mpz_clear(expand_new);
    mpz_clear(expanded);
}

void init(void) {
    clock_tick = sysconf(_SC_CLK_TCK);
    init_seen();
    init_array(&curpend, sizeof(mp_limb_t), 1 << 20);
    init_array(&newpend, sizeof(mp_limb_t), 1 << 20);
    pendoff = 0;
    mpz_init(expanded);
    mpz_init(expand_new);
}

/*
 * The specified integer hasn't been seen before, so put it on the next
 * generation's pending list
 */
void pend(mpz_t n) {
    ulong limbs = mpz_size(n);
    mp_limb_t *np;
    resize_array(&newpend, newpend.count + limbs + 1);
    np = (mp_limb_t *)array_element(&newpend, newpend.count);
    *np++ = (mp_limb_t)limbs;
    mpz_export((void*)np, NULL, 1, sizeof(mp_limb_t), 0, 0, n);
    newpend.count += limbs + 1;
}

/*
 * Pick the next integer from this generation's pending list and store it
 * into the supplied integer.
 *
 * If the list is empty, advance the generation (showing stats), and store
 * the first integer of the new generation's pending list instead.
 *
 * Returns 1 if it managed to find an integer to store, 0 if there's nothing
 * left to do.
 */
int nextpend(mpz_t n) {
    if (pendoff < curpend.count) {
        mp_limb_t *np = (mp_limb_t *)array_element(&curpend, pendoff);
        ulong size = (ulong)(*np++);
        mpz_import(n, (size_t)size, 1, sizeof(mp_limb_t), 0, 0, (void *)np);
        pendoff += size + 1;
        return 1;
    }
    dumpstats();
    {
        array_t tmp = curpend;
        curpend = newpend;
        newpend = tmp;
        newpend.count = 0;
        pendoff = 0;
    }
    if (pendoff < curpend.count) {
        mp_limb_t *np = (mp_limb_t *)array_element(&curpend, pendoff);
        ulong size = (ulong)(*np++);
        mpz_import(n, (size_t)size, 1, sizeof(mp_limb_t), 0, 0, (void *)np);
        pendoff += size + 1;
        return 1;
    }
    return 0;
}

/*
 * Process a new integer: apply the calculation, write it in the requested
 * base, split it on zero digits, and for each of the resulting substrings
 * check if we've seen it before and if not, mark it as seen and store it
 * on the pending list for the next generation.
 *
 * TODO: do the base expansion manually, to avoid GMP's limit of base 62.
 */
void expand(mpz_t n) {
    char* bn2;
    char* bnc;
    char* bnt;
    int done = 0;

    /* apply the calculation on <n>, leaving the result in <expanded> */
    CALCULATE();

    /* this malloc()s a string we'll need to free */
    bnc = bn2 = mpz_get_str((char*)NULL, (int)base, expanded);
    if (DEBUG) {
        char* s = mpz_get_str((char*)NULL, 10, n);
        char* t = mpz_get_str((char*)NULL, 10, expanded);
        printf("expand %s to %s_%u (%s_10)\n", s, bn2, base, t);
        free(t);
        free(s);
    }
    while (!done) {
        bnt = bnc;
        while (*bnt != '0' && *bnt) ++bnt;
        if (*bnt == '0')
            *bnt = '\0';
        else
            done = 1;
        if (*bnc) {
            /*
             * The meat: we've found a substring, so evaluate it, check if
             * we've seen it before, and account for it if we haven't.
             */
            mpz_set_str(expand_new, bnc, (int)base);
            if (!seen(expand_new)) {
                ++total;
                pend(expand_new);
            }
        }
        bnc = ++bnt;
    }
    free(bn2);
}

int main(int argc, char** argv) {
    mpz_t start;
    mpz_t next;
    if (argc > 1) {
        base = (uint)atoi(argv[1]);
    }

    init();
    setvbuf(stdout, (char*)NULL, _IONBF, 0);
    mpz_init_set_ui(start, 2ul);
    seen(start);
    ++total;
    pend(start);
    mpz_init(next);
    while (nextpend(next)) {
        expand(next);
    }
    final();
    mpz_clear(next);
    mpz_clear(start);
    return 0;
}

/*
Results for calculation s->s^2:
n=2: count 2, max 2
n=3: count 18, max 1849
n=4: count 2, max 2
n=5: count 3050, max 266423227914725931
n=6: count 34762, max 3100840870711697060720215047
n=7: count 3087549, max 845486430620513036335402848567278325780455810752216401
n=8: count 2, max 4
n=9: after 74 generations:
  - 434,155,947 distinct values seen
  - 37,990,284 new values on the next generation's pending list
  - pending list growing by 6.1% per generation
  - growth rate reducing by about 0.4 percentage points per generation
  - second derivative is also slowly reducing
  - number of values in each hash, by number of limbs:
      1:262,230,247 2:53,981,017 3:801,061 4:9,405 5:87
  - the rest (117,134,130) are in the bit vector for small ints
*/
