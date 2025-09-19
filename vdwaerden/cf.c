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
#include "allsub.h"

/* we may search for f(n): n <= MAXN */
#define MAXN 254
typedef uint16_t n_t;
static_assert(MAXN < (1 << (sizeof(n_t) * 8)), "n_t too small for MAXN");

typedef uint32_t scratch_t;
static_assert(MAXN <= (UINT32_MAX >> 1), "scratch_t too small for 2 * MAXN");

/* we may search for f(n): f(n) <= MAXF */
#define MAXF 255
typedef uint8_t f_t;
static_assert(MAXF < (1 << (sizeof(f_t) * 8)), "f_t too small for MAXF");

n_t n;              /* current target n */
n_t debug_n = 0;    /* show debug info when n == debug_n */
f_t max;            /* cardinality of largest subset found so far */
                    /* NB f(n) must be either max or max+1 */
f_t f3[MAXN + 1];   /* known values of f() */
ulong iter;         /* number of iterations, for debugging/stats only */

double t0 = 0;      /* base for timing calculations */
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

uint8_t drev8[1 << 8];

#define CHUNK_BYTES 1
#define CHUNK_BITS (CHUNK_BYTES << 3)
#define CHUNK_MASK ((1 << CHUNK_BITS) - 1)
#if CHUNK_BYTES == 1
    typedef uint8_t chunk_t;
#   define revchunk(x) rev8(x)
#elif CHUNK_BYTES == 2
    typedef uint16_t chunk_t;
#   define revchunk(x) rev16(x)
#elif (CHUNK_BYTES == 3) | (CHUNK_BYTES == 4)
    typedef uint32_t chunk_t;
#   define revchunk(x) rev32(x)
#else
#   error "invalid CHUNK_BYTES"
#endif

typedef struct chunkset_s {
    size_t count;
    chunk_t *cp;
} chunkset_t;
chunkset_t chunks[CHUNK_BITS + 1];
chunkset_t chunks_bottom[CHUNK_BITS + 1];
chunkset_t chunks_top[CHUNK_BITS + 1];
chunk_t self_prod[1 << CHUNK_BITS];
chunk_t self_prod_reverse[1 << CHUNK_BITS];
chunk_t cross_prod[3 * (1 << (CHUNK_BITS * 2))];

typedef struct stack_s {
    chunkset_t *csp;
    size_t chunk_off;
    uint bitcount;
    chunk_t exclude[MAXN / (sizeof(chunk_t) << 3) + 1];
} stack_t;
stack_t stack[MAXN / CHUNK_BITS + 1];
signed int spi;

#define LOOKUP_BYTES CHUNK_BYTES
#define LOOKUP_BITS (LOOKUP_BYTES << 3)
#define LOOKUP_SIZE (1 << LOOKUP_BITS)
#define LOOKUP_MASK ((1 << LOOKUP_BITS) - 1)
f_t lookup_sub[LOOKUP_SIZE];        /* free choice */

void init_rev(void) {
    memset(drev8, 0, sizeof(drev8));
    for (uint bi = 0; bi < 8; ++bi) {
        uint b = 1 << bi;
        uint br = 1 << (7 - bi);
        for (uint j = 0; j < (1 << 8); ++j)
            if (j & b)
                drev8[j] |= br;
    }
}

static inline uint8_t rev8(uint8_t x) {
    return drev8[x];
}

static inline uint16_t rev16(uint16_t x) {
    return (rev8(x & 0xff) << 8) | rev8(x >> 8);
}

static inline uint32_t rev24(uint32_t x) {
    return (rev8(x & 0xff) << 16) | (rev8((x >> 8) & 0xff) << 8)
            | rev8(x >> 16);
}

static inline uint32_t rev32(uint32_t x) {
    return (rev16(x & 0xffff) << 16) | rev16(x >> 16);
}

int cmp_lexical(const void *va, const void *vb) {
    chunk_t ra = revchunk(*(chunk_t *)va);
    chunk_t rb = revchunk(*(chunk_t *)vb);
    return (ra < rb) ? 1 : -1;
}

static inline void cross_product(
    chunk_t *result, chunk_t left, chunk_t right
) {
    ullong rleft = (ullong)revchunk(left);
    while (right) {
        uint bi = lsb32((uint32_t)right);
        right ^= 1 << bi;
        uint shift = (bi << 1) + 1;
        result[0] |= (chunk_t)(
            (rleft << shift) & CHUNK_MASK
        );
        if (shift > CHUNK_BITS)
            result[1] |= (chunk_t)(
                (rleft << (shift - CHUNK_BITS)) & CHUNK_MASK
            );
        else
            result[1] |= (chunk_t)(
                (rleft >> (CHUNK_BITS - shift)) & CHUNK_MASK
            );
        result[2] |= (chunk_t)(
            (rleft >> ((CHUNK_BITS << 1) - shift)) & CHUNK_MASK
        );
    }
    return;
}

static inline chunk_t self_product(chunk_t right) {
    ullong rleft = (ullong)revchunk(right);
    chunk_t result = 0;
    while (right) {
        uint bi = lsb32((uint32_t)right);
        right ^= 1 << bi;
        uint shift = (bi << 1) + 1;
        uint revi = CHUNK_BITS - bi;
        ullong masked = rleft & ~((1ULL << revi) - 1);
        result |= (chunk_t)(
            (masked >> ((CHUNK_BITS << 1) - shift)) & CHUNK_MASK
        );
    }
    return result;
}

void init_chunks(void) {
    for (uint i = 0; i <= CHUNK_BITS; ++i) {
        size_t count = (size_t)count_allsub(i);
        chunks[i].count = count;
        chunks[i].cp = malloc(count * sizeof(chunk_t));
        /* chunks with bottom bit set all come lexically before those
         * without, so we share the array just with a lower count.
         */
        chunks_bottom[i].cp = chunks[i].cp;
        /* By symmetry, there are as many with bottom bit as with top bit
         * set, and by definition count(n with top bit set) is
         * count(n) - count(n - 1).
         */
        if (i) {
            size_t new = count - chunks[i - 1].count;
            chunks_bottom[i].count = new;
            chunks_top[i].count = new;
            chunks_top[i].cp = malloc(new * sizeof(chunk_t));
        }
    }
    size_t ci[CHUNK_BITS + 1], cti[CHUNK_BITS + 1];
    memset(ci, 0, sizeof(ci));
    memset(cti, 0, sizeof(cti));

    allsub_iter_t *aip = init_iter(CHUNK_BITS);
    chunkset_t *csp = &chunks[CHUNK_BITS];
    while (1) {
        chunk_t val = (chunk_t)do_iter(aip);
        csp->cp[ ci[CHUNK_BITS]++ ] = val;
        if (val == 0)
            break;
    }
    done_iter(aip);

    /* sort once, then fill all the derived sets in sorted order */
    qsort(csp->cp, csp->count, sizeof(chunk_t), cmp_lexical);

    for (size_t off = 0; off < csp->count; ++off) {
        chunk_t val = csp->cp[off];
        uint bi = val ? msb32((uint)val) + 1 : 0;
        for (uint i = bi; i < CHUNK_BITS; ++i)
            chunks[i].cp[ ci[i]++ ] = val;
        if (bi)
            chunks_top[bi].cp[ cti[bi]++ ] = val;
    }

    /* sanity check */
    for (uint i = 0; i <= CHUNK_BITS; ++i) {
        if (ci[i] != chunks[i].count) {
            fprintf(stderr, "for size %u found %llu expecting %llu\n",
                    i, (ullong)ci[i], (ullong)chunks[i].count);
            exit(1);
        }
        if (cti[i] != chunks_top[i].count) {
            fprintf(stderr, "for top size %u found %llu expecting %llu\n",
                    i, (ullong)cti[i], (ullong)chunks_top[i].count);
            exit(1);
        }
    }

    static_assert(sizeof(uint) > sizeof(chunk_t), "our loop must end");
    for (uint i = 0; i <= (1 << CHUNK_BITS); ++i) {
        self_prod[i] = self_product((chunk_t)i);
        for (uint j = 0; j <= (1 << CHUNK_BITS); ++j) {
            chunk_t *cp = &cross_prod[3 * ((i << CHUNK_BITS) + j)];
            cross_product(cp, (chunk_t)i, (chunk_t)j);
        }
    }
    for (uint i = 0; i <= (1 << CHUNK_BITS); ++i)
        self_prod_reverse[i] = revchunk(self_prod[revchunk((chunk_t)i)]);
}

bool better_sub(uint32_t avail, uint32_t have, f_t need) {
    while (1) {
        if (lookup_sub[avail] < need)
            return false;
        uint bi = msb32(avail);
        uint32_t b = 1 << bi;
        uint32_t have2 = have | b;
        uint32_t avail2 = avail ^ b;
        f_t need2 = need - 1;
        uint32_t exclude = rev32(have >> (bi + 1))
                >> (UINT32_HIGHBIT - (bi - 1));
        avail2 &= ~exclude;
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
    for (uint32_t i = 1; i < LOOKUP_SIZE; ++i) {
        /* free choice */
        int bits = __builtin_popcount(i);
        if (bits <= 2) {
            lookup_sub[i] = bits;
        } else if ((i & 1) == 0) {
            lookup_sub[i] = lookup_sub[i >> 1];
        } else {
            uint32_t b = 1 << msb32(i);
            uint32_t j = i ^ b;
            f_t prev = lookup_sub[j];
            lookup_sub[i] = better_sub(j, b, prev)
                    ? prev + 1 : prev;
        }
    }

    /* our actual data has bits set for _disallowed_ values, so invert */
    for (uint32_t i = 0; i < LOOKUP_SIZE; ++i) {
        uint32_t j = (~i) & LOOKUP_MASK;
        if (i < j) {
            f_t temp = lookup_sub[i];
            lookup_sub[i] = lookup_sub[j];
            lookup_sub[j] = temp;
        }
    }
}

void reset_data(uint n) {
    uint final = (n - 1) / CHUNK_BITS;
    /* will overwrite below if this chunk is also final */
    stack[0].csp = &chunks_bottom[CHUNK_BITS];
    for (uint i = 1; i < final; ++i)
        stack[i].csp = &chunks[CHUNK_BITS];

    /* For n <= chunk bits, we set stack[0] to fix bottom bit rather than
     * top bit to save a test in the main loop. It does not seem worth an
     * extra chunks_bottom_top array just to save a few iterations for
     * small n.
     */
    uint fragchunk = ((n - 1) % CHUNK_BITS) + 1;
    stack[final].csp = (final == 0)
            ? &chunks_bottom[fragchunk] : &chunks_top[fragchunk];

    stack[0].chunk_off = 0;
    stack[0].bitcount = 0;
    bzero(&stack[0].exclude, sizeof(stack[0].exclude));
    spi = 0;
}

char bufstack[4096];
char *show_stack(uint sp) {
    uint off = 0;
    for (uint si = 0; si <= sp; ++si) {
        stack_t *sp = &stack[si];
        chunkset_t *csp = sp->csp;
        chunk_t val = csp->cp[ sp->chunk_off ];
        while (val) {
            uint bi = lsb32((uint32_t)val);
            val ^= (1 << bi);
            if (off)
                off += snprintf(&bufstack[off], sizeof(bufstack) - off, " ");
            off += snprintf(&bufstack[off], sizeof(bufstack) - off,
                    "%u", si * CHUNK_BITS + bi + 1);
        }
    }
    return &bufstack[0];
}

char bufshort[4096];
char *show_short(uint sp) {
    uint off = 0;
    for (uint si = 0; si <= sp; ++si) {
        stack_t *sp = &stack[si];
        chunkset_t *csp = sp->csp;
        chunk_t val = csp->cp[ sp->chunk_off ];
        if (off)
            off += snprintf(&bufshort[off], sizeof(bufshort) - off, " ");
        off += snprintf(&bufshort[off], sizeof(bufshort) - off,
                "%0*x", CHUNK_BITS / 4, val);
    }
    return &bufshort[0];
}

char bufexcl[4096];
char *show_excl(uint sp) {
    uint off = 0;
    chunk_t *excl = &(stack[sp].exclude[0]);
    uint lim = (n - 1) / CHUNK_BITS;
    for (uint ci = 0; ci <= lim; ++ci) {
        chunk_t c = excl[ci];
        if (off)
            off += snprintf(&bufexcl[off], sizeof(bufexcl) - off, " ");
        off += snprintf(&bufexcl[off], sizeof(bufexcl) - off, 
                "%0*x", CHUNK_BITS / 4, c);
    }
    return &bufexcl[0];
}

void disp_result(bool ok) {
    printf("%u (%lu %.2fs): ", n, iter, utime() - t0);
    if (ok)
        printf("[%s]\n", show_stack((n - 1) / CHUNK_BITS));
    else
        printf("no improvement\n");
}

static inline n_t bits_avail(chunk_t v) {
    return lookup_sub[v];
}

void findmax(void) {
    n = 2;
    f3[0] = 0;
    f3[1] = 1;
    max = 1;
  TRY_NEXT:
#ifdef DEBUG
    if (debug_n && n > debug_n)
        return;
#endif
    iter = 0;
    reset_data(n);
    uint final = (n - 1) / CHUNK_BITS;
    uint finalbit = (n - 1) % CHUNK_BITS;
    uint finalmask = (uint)((1ULL << (finalbit + 1)) - 1);

    while (spi >= 0) {
        ++iter;
        stack_t *sp = &stack[spi];
        chunkset_t *csp = sp->csp;
        size_t chunk_off = sp->chunk_off;

        if (chunk_off >= csp->count)
            goto derecurse;
        chunk_t val = csp->cp[chunk_off];
#ifdef DEBUG
        if (debug_n && n == debug_n)
            fprintf(stderr, "> <%s> [%s] {%s}\n",
                    show_short(spi), show_stack(spi), show_excl(spi));
#endif

        if (val & sp->exclude[spi])
            goto next_iter;
        /* In the first extent we want to force 1 in the set; lexical
         * ordering ensures those are grouped at the start.
         */
        if (spi == 0 && !(val & 1))
            goto derecurse;

        uint bitcount = sp->bitcount + __builtin_popcount((uint32_t)val);
        if (spi == final) {
            if (bitcount <= max)
                goto next_iter;
            if (spi > 0) {
                /* Our exclude mask told us if any preceding pair of bits
                 * disallowed a bit in val, now check if a pair of bits in
                 * val conflicts with a single bit preceding.
                 */
                stack_t *prev = &stack[spi - 1];
                chunk_t mask = self_prod_reverse[val];
                if (mask & prev->csp->cp[prev->chunk_off])
                    goto next_iter;
            }
            /* we have a solution */
            disp_result(1);
            ++max;
            f3[n] = max;
            ++n;
            /* FIXME: don't try to continue from where we are, for now */
            goto TRY_NEXT;
        }

        if (bitcount + f3[n - (spi + 1) * CHUNK_BITS] <= max)
            goto next_iter;

        /* ok to continue with this subset */

        /* block the third element of any new 2-element APs */
        stack_t *nextsp = &stack[spi + 1];
        chunkset_t *nextcsp = nextsp->csp;
        chunk_t *exclude = nextsp->exclude;
/* FIXME:
 * copy only up to n;
 * stop loop when low bit would exceed n;
 * ensure we have room for 3 chunks
 */
        memcpy(exclude, &sp->exclude, sizeof(sp->exclude));
        exclude[spi + 1] |= self_prod[val];
        for (uint rel = 1; rel <= spi; ++rel) {
            stack_t *relsp = &stack[spi - rel];
            chunkset_t *relcsp = relsp->csp;
            chunk_t *mask = &cross_prod[
                3 * ((relcsp->cp[relsp->chunk_off] << CHUNK_BITS) + val)
            ];
            for (uint off = 0; off < 3; ++off)
                exclude[spi + rel - 1 + off] |= mask[off];
        }
        if (val & exclude[spi])
            goto next_iter;
        exclude[final] &= finalmask;
        /* better solution is not possible unless it include both 1 and n */
        if (exclude[final] & (1 << finalbit))
            goto next_iter;

        uint maxrest = 0;
        for (uint i = spi + 1; i <= final; ++i) {
            chunk_t v = exclude[i];
            if (i == final)
                v |= ~finalmask;
            maxrest += bits_avail(v);
        }

        if (bitcount + maxrest <= max)
            goto next_iter;

        /* advance for next element */
        ++spi;
        nextsp->chunk_off = 0;
        nextsp->bitcount = bitcount;
        continue;

      derecurse:
        --spi;
        sp = &stack[spi];
      next_iter:
        ++sp->chunk_off;
        continue;
    }
    /* no solution: f(n) = f(n-1) */
    disp_result(0);
    f3[n] = max;
    ++n;
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
    init_rev();
    init_lookup_sub();
    init_chunks();
    t0 = utime();
    findmax();
    return 0;
}
