#ifndef MYGMP_H
#define MYGMP_H

#include <gmp.h>

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

extern inline void QINIT(mpq_t* q, char* format, ...) {
	mpq_init(*q);
	DEBUG_PRINT(q, "qinit", format);
}

extern inline void QCLEAR(mpq_t* q, char* format, ...) {
	mpq_clear(*q);
	DEBUG_PRINT(q, "qclear", format);
}

extern inline void ZINIT(mpz_t* z, char* format, ...) {
	mpz_init(*z);
	DEBUG_PRINT(z, "zinit", format);
}

extern inline void ZCLEAR(mpz_t* z, char* format, ...) {
	mpz_clear(*z);
	DEBUG_PRINT(z, "zclear", format);
}

#define MPX_MAXLIMBS 8
typedef mp_limb_t* mpx_t;
typedef void (mpx_add_func)(mpx_t xd, mpx_t xs1, mpx_t xs2);
typedef int (mpx_cmp_func)(mpx_t xs1, mpx_t xs2);
typedef struct mpx_s_support {
	int size;
	mpx_add_func* adder;
	mpx_cmp_func* cmper;
} mpx_support;

extern mpx_support* mpx_support_n(int n);
extern mpx_support* mpx_support_z(mpz_t z);
extern void mpx_set_ui(mpx_t x, int size, unsigned int ui);
extern void mpx_set(mpx_t xd, int sized, mpx_t xs, int sizes);
extern void mpx_set_z(mpx_t x, int size, mpz_t z);
extern void mpz_set_x(mpz_t z, mpx_t x, int size);

#endif /* MYGMP_H */
