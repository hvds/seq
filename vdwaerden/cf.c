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

static_assert(CHUNK_BYTES == 1, "required for MAX_CHUNK_SET");
#define MAX_CHUNK_SET 4

chunk_t *allchunks;
size_t allchunks_size;
size_t allchunks_i;
typedef struct sized_extent_s {
    size_t start;
    size_t end;
} sized_extent_t;
typedef struct extent_s {
    uint maxbits;
    sized_extent_t range[MAX_CHUNK_SET + 1];
} extent_t;

/* An extent for each mask/reverse mask combination. */
extent_t extents[1 << (2 * CHUNK_BITS)];
/* An extent with bit 0 set, with zero mask and reverse mask, for each width
 * (with msb set), and one for full width. */
extent_t start_extents[CHUNK_BITS + 1];
/* An extent for each mask/reverse mask combination, for each width,
 * with msb set. */
extent_t tail_extents[CHUNK_BITS * (1 << (2 * CHUNK_BITS))];

chunk_t self_prod[1 << CHUNK_BITS];
chunk_t self_prod_reverse[1 << CHUNK_BITS];
chunk_t cross_prod[3 * (1 << (CHUNK_BITS * 2))];

typedef struct stack_s {
    size_t off;
    size_t end;
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

static inline void push_chunk(chunk_t val) {
    if (allchunks_i >= allchunks_size) {
        size_t newsize = (allchunks_size) ? (allchunks_size * 3 / 2) : 1024;
        allchunks = realloc(allchunks, newsize * sizeof(chunk_t));
        allchunks_size = newsize;
    }
    allchunks[allchunks_i++] = val;
}

static inline void fill_extent(extent_t *ep) {
    uint maxbits = ep->maxbits;
    if (maxbits == 0)
        return;
    size_t start = ep->range[0].start;
    size_t end = ep->range[0].end;
    /* lexically sorted, 0 will be at the end */
    ep->range[1].start = start;
    ep->range[1].end = (allchunks[end - 1] == 0) ? end - 1 : end;

    for (uint bitc = 2; bitc <= maxbits; ++bitc) {
        size_t next = allchunks_i;
        for (size_t off = start; off < end; ++off) {
            chunk_t val = allchunks[off];
            if (__builtin_popcount((uint32_t)val) < bitc)
                continue;
            push_chunk(val);
        }
        ep->range[bitc].start = start = next;
        ep->range[bitc].end = end = allchunks_i;
    }
}

/* For middle chunks, we will specify the mask (2^8), the reverse
 * mask (2^8) and the required minimum number of bits (0 .. 4)
 * for a total of 5 x 2^16 entries.
 * For the first chunk, the mask and reverse mask are always zero, but
 * we will specify the required minimum number of bits (0 .. 4) and
 * require that all entries set the lsb.
 * For the last chunk, we will specify the mask, reverse mask and minbits,
 * and require that all entries set the (floating) msb, for a total of
 * 5 x 2^19 entries.
 */
void init_chunks(void) {
    for (uint i = 0; i <= CHUNK_MASK; ++i) {
        self_prod[i] = self_product((chunk_t)i);
        for (uint j = 0; j <= CHUNK_MASK; ++j) {
            chunk_t *cp = &cross_prod[3 * ((i << CHUNK_BITS) | j)];
            cross_product(cp, (chunk_t)i, (chunk_t)j);
        }
    }
    for (uint i = 0; i <= CHUNK_MASK; ++i)
        self_prod_reverse[i] = revchunk(self_prod[revchunk((chunk_t)i)]);

    uint maxbits0 = 0;
    size_t start0 = allchunks_i;
    allsub_iter_t *aip = init_iter(CHUNK_BITS);
    while (1) {
        chunk_t val = (chunk_t)do_iter(aip);
        push_chunk(val);
        uint bits = __builtin_popcount((uint32_t)val);
        if (maxbits0 < bits)
            maxbits0 = bits;
        if (val == 0)
            break;
    }
    done_iter(aip);
    size_t end0 = allchunks_i;

    /* sort once, then fill all the derived sets in sorted order */
    qsort(&allchunks[0], end0 - start0, sizeof(chunk_t), cmp_lexical);

    {
        extent_t *ep = &extents[0];
        ep->maxbits = maxbits0;
        ep->range[0].start = start0;
        ep->range[0].end = end0;
        fill_extent(ep);
    }

    for (uint mask = 0; mask <= CHUNK_MASK; ++mask) {
        for (uint revmask = 0; revmask <= CHUNK_MASK; ++revmask) {
            if (mask == 0 && revmask == 0)
                continue;
            extent_t *ep = &extents[(mask << CHUNK_BITS) | revmask];
            uint maxbits = 0;
            size_t start = allchunks_i;
            for (size_t off = start0; off < end0; ++off) {
                chunk_t val = allchunks[off];
                if (val & (chunk_t)mask)
                    continue;
                if (self_prod_reverse[val] & (chunk_t)revmask)
                    continue;
                push_chunk(val);
                uint bits = __builtin_popcount((uint32_t)val);
                if (maxbits < bits)
                    maxbits = bits;
            }
            ep->maxbits = maxbits;
            ep->range[0].start = start;
            ep->range[0].end = allchunks_i;
            fill_extent(ep);
        }
    }

    {
        extent_t *ep = &start_extents[CHUNK_BITS];
        uint maxbits = 0;
        size_t start = allchunks_i;
        for (size_t off = start0; off < end0; ++off) {
            chunk_t val = allchunks[off];
            if ((val & 1) == 0)
                continue;
            push_chunk(val);
            uint bits = __builtin_popcount((uint32_t)val);
            if (maxbits < bits)
                maxbits = bits;
        }
        ep->maxbits = maxbits;
        ep->range[0].start = start;
        ep->range[0].end = allchunks_i;
        fill_extent(ep);
    }

    for (uint width = 0; width < CHUNK_BITS; ++width) {
        chunk_t expect = (1 << width);
        chunk_t matchmask = 1 + ~expect;
        extent_t *ep = &start_extents[width];
        uint maxbits = 0;
        size_t start = allchunks_i;
        for (size_t off = start0; off < end0; ++off) {
            chunk_t val = allchunks[off];
            if ((val & 1) == 0)
                continue;
            if ((val & matchmask) != expect)
                continue;
            push_chunk(val);
            uint bits = __builtin_popcount((uint32_t)val);
            if (maxbits < bits)
                maxbits = bits;
        }
        ep->maxbits = maxbits;
        ep->range[0].start = start;
        ep->range[0].end = allchunks_i;
        fill_extent(ep);
    }

    for (uint width = 0; width < CHUNK_BITS; ++width) {
        chunk_t expect = (1 << width);
        chunk_t matchmask = 1 + ~expect;
        for (uint mask = 0; mask <= CHUNK_MASK; ++mask) {
            for (uint revmask = 0; revmask <= CHUNK_MASK; ++revmask) {
                extent_t *sep = &extents[(mask << CHUNK_BITS) | revmask];
                size_t sstart = sep->range[0].start;
                size_t send = sep->range[0].end;
                extent_t *ep = &tail_extents[
                    (((width << CHUNK_BITS) | mask) << CHUNK_BITS) | revmask
                ];
                uint maxbits = 0;
                size_t start = allchunks_i;
                for (size_t off = sstart; off < send; ++off) {
                    chunk_t val = allchunks[off];
                    if ((val & matchmask) != expect)
                        continue;
                    push_chunk(val);
                    uint bits = __builtin_popcount((uint32_t)val);
                    if (maxbits < bits)
                        maxbits = bits;
                }
                ep->maxbits = maxbits;
                ep->range[0].start = start;
                ep->range[0].end = allchunks_i;
                fill_extent(ep);
            }
        }
    }
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

static inline uint min_needed(uint bitcount, chunk_t *exclude, chunk_t mask) {
    uint final = (n - 1) / CHUNK_BITS;
    uint finalbit = (n - 1) % CHUNK_BITS;
    uint bits_needed = max + 1 - bitcount;
    uint max_avail = bits_needed;
    uint known_avail = 0;
    for (uint i = spi + 1; i <= final; ++i) {
        signed int theoretical = f3[n - i * CHUNK_BITS];
        if (max_avail > theoretical + known_avail)
            max_avail = theoretical + known_avail;
        chunk_t mask = exclude[i];
        extent_t *aep = (i == final)
            ? &tail_extents[(((finalbit << CHUNK_BITS) | mask) << CHUNK_BITS) | 0]
            : &extents[(mask << CHUNK_BITS) | 0];
        known_avail += aep->maxbits;
    }
    if (max_avail > known_avail)
        max_avail = known_avail;
    return (bits_needed > max_avail) ? bits_needed - max_avail : 0;
}

void reset_data(uint n) {
    uint final = (n - 1) / CHUNK_BITS;
    bzero(&stack[0].exclude, sizeof(stack[0].exclude));
    spi = 0;
    uint need = min_needed(0, stack[0].exclude, (chunk_t)0);
    extent_t *ep = &start_extents[
        (n <= CHUNK_BITS) ? (n - 1) : CHUNK_BITS
    ];
    if (need > ep->maxbits)
        need = ep->maxbits;
    stack[0].off = ep->range[need].start;
    stack[0].end = ep->range[need].end;
    stack[0].bitcount = 0;
}

char bufstack[4096];
char *show_stack(uint sp) {
    uint off = 0;
    for (uint si = 0; si <= sp; ++si) {
        stack_t *sp = &stack[si];
        chunk_t val = allchunks[sp->off];
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
        chunk_t val = allchunks[sp->off];
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
        size_t off = sp->off;

        if (off >= sp->end)
            goto derecurse;
        chunk_t val = allchunks[off];
#ifdef DEBUG
        if (debug_n && n == debug_n)
            fprintf(stderr, "> <%s> [%s] {%s}\n",
                    show_short(spi), show_stack(spi), show_excl(spi));
#endif

        uint bitcount = sp->bitcount + __builtin_popcount((uint32_t)val);
        if (spi == final) {
            if (bitcount <= max)
                goto next_iter;
            /* we have a solution */
            disp_result(1);
            ++max;
            f3[n] = max;
            ++n;
            goto TRY_NEXT;
        }

        if (bitcount + f3[n - (spi + 1) * CHUNK_BITS] <= max)
            goto next_iter;

        /* ok to continue with this subset */

        /* block the third element of any new 2-element APs */
        stack_t *nextsp = &stack[spi + 1];
        chunk_t *exclude = nextsp->exclude;
/* FIXME:
 * copy only from spi+1 to final;
 * stop loop when low bit would exceed final;
 * apply only up to final
 */
        memcpy(exclude, &sp->exclude, sizeof(sp->exclude));
        exclude[spi + 1] |= self_prod[val];
        for (uint rel = 1; rel <= spi; ++rel) {
            stack_t *relsp = &stack[spi - rel];
            chunk_t relval = allchunks[relsp->off];
            chunk_t *mask = &cross_prod[3 * ((relval << CHUNK_BITS) | val)];
            for (uint eoff = 0; eoff < 3; ++eoff)
                exclude[spi + rel - 1 + eoff] |= mask[eoff];
        }
        exclude[final] &= finalmask;

        /* advance for next element */
        ++spi;
        extent_t *ep;
        uint mask = (uint)exclude[spi];
        if (spi == final) {
            ep = &tail_extents[(((finalbit << CHUNK_BITS) | mask) << CHUNK_BITS) | val];
        } else {
            ep = &extents[(mask << CHUNK_BITS) | val];
        }
        uint need = min_needed(bitcount, exclude, mask);
        if (need > ep->maxbits)
            goto derecurse;
        nextsp->off = ep->range[need].start;
        nextsp->end = ep->range[need].end;
        nextsp->bitcount = bitcount;
        continue;

      derecurse:
        --spi;
        sp = &stack[spi];
      next_iter:
        ++sp->off;
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
