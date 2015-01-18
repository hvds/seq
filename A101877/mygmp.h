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

#endif /* MYGMP_H */
