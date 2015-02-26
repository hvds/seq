#ifndef MYGMP_H
#define MYGMP_H

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <gmp.h>

typedef unsigned long ulong;
typedef int bool;

#ifdef IN_MYGMP_C
#define INLINABLE inline
#else
#define INLINABLE extern inline
#endif

#ifdef DEBUG_GMP_LEAK
#include <stdio.h>
#include <stdarg.h>
#define DEBUG_PRINT(q, string, format) { \
    va_list ap; \
    va_start(ap, format); \
    fprintf(stderr, "%s %p: ", string, q); \
    vfprintf(stderr, format, ap); \
    fprintf(stderr, "\n"); \
    va_end(ap); \
}
#else /* DEBUG_GMP_LEAK */
#define DEBUG_PRINT(q, string, format)
#endif /* DEBUG_GMP_LEAK */

INLINABLE void QINIT(mpq_t* q, char* format, ...) {
    mpq_init(*q);
    DEBUG_PRINT(q, "qinit", format);
}

INLINABLE void QCLEAR(mpq_t* q, char* format, ...) {
    mpq_clear(*q);
    DEBUG_PRINT(q, "qclear", format);
}

INLINABLE void ZINIT(mpz_t* z, char* format, ...) {
    mpz_init(*z);
    DEBUG_PRINT(z, "zinit", format);
}

INLINABLE void ZCLEAR(mpz_t* z, char* format, ...) {
    mpz_clear(*z);
    DEBUG_PRINT(z, "zclear", format);
}

INLINABLE void mpq_add_ui(mpq_t q, unsigned long u) {
    mpz_addmul_ui(mpq_numref(q), mpq_denref(q), u);
}

typedef int pack_t; /* assume mp_limb_t and mp_size_t are multiples of this */
#define LIMB_MULT (sizeof(mp_limb_t) / sizeof(pack_t))

typedef struct {
    mp_size_t numsize;
    mp_size_t densize;
    mp_limb_t limbs[0];
} mpq_pack_t;
#define PACKSIZE (sizeof(mpq_pack_t) / sizeof(pack_t))

INLINABLE size_t mpq_packsize(mpq_t q) {
    return PACKSIZE + (mpz_size(mpq_numref(q)) + mpz_size(mpq_denref(q))) * LIMB_MULT;
}

INLINABLE size_t packed_size(mpq_pack_t *p) {
    return PACKSIZE + (p->numsize + p->densize) * LIMB_MULT;
}

INLINABLE void store_packed(mpq_pack_t *p, mpq_t q) {
    p->numsize = mpz_size(mpq_numref(q));
    p->densize = mpz_size(mpq_denref(q));
    mpz_export((void*)&p->limbs[0],
            NULL, 1, sizeof(mp_limb_t), 0, 0, mpq_numref(q));
    mpz_export((void*)&p->limbs[p->numsize],
            NULL, 1, sizeof(mp_limb_t), 0, 0, mpq_denref(q));
}

INLINABLE void fetch_packed(mpq_pack_t *p, mpq_t q) {
    mpz_import(mpq_numref(q), (size_t)p->numsize,
            1, sizeof(mp_limb_t), 0, 0, (void*)&p->limbs[0]);
    mpz_import(mpq_denref(q), (size_t)p->densize,
            1, sizeof(mp_limb_t), 0, 0, (void*)&(p->limbs[p->numsize]));
}

INLINABLE ulong mpz_strict_get_ui(mpz_t z) {
    if (mpz_sgn(z) < 0 || mpz_cmp_ui(z, ULONG_MAX) > 0) {
        gmp_fprintf(stderr, "strict ulong overflow converting %Zd\n", z);
        exit(1);
    }
    return mpz_get_ui(z);
}

INLINABLE ulong mpz_bitsize(mpz_t z) {
    size_t s = mpz_sizeinbase(z, 2);
    return (ulong)s;
}

INLINABLE bool mpq_simpler(mpq_t q1, mpq_t q2) {
    int i = mpz_cmp(mpq_numref(q1), mpq_numref(q2));
    if (i == 0) {
        return (mpz_cmp(mpq_denref(q1), mpq_numref(q2))) < 0 ? 1 : 0;
    } else {
        return (i < 0) ? 1 : 0;
    }
}

#endif /* MYGMP_H */
