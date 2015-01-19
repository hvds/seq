#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mygmp.h"

inline void QINIT(mpq_t* q, char* format, ...) {
	mpq_init(*q);
	DEBUG_PRINT(q, "qinit", format);
}

inline void QCLEAR(mpq_t* q, char* format, ...) {
	mpq_clear(*q);
	DEBUG_PRINT(q, "qclear", format);
}

inline void ZINIT(mpz_t* z, char* format, ...) {
	mpz_init(*z);
	DEBUG_PRINT(z, "zinit", format);
}

inline void ZCLEAR(mpz_t* z, char* format, ...) {
	mpz_clear(*z);
	DEBUG_PRINT(z, "zclear", format);
}

/*
 * Set a fixed-size integer from an unsigned int.
 * Inputs:
 *   mpx_t x: allocated fixed-size integer to write to
 *   int size: size of x (in mp_limb_t)
 *   unsigned int ui: the value to set
 * Returns:
 *   Nothing.
 */
void mpx_set_ui(mpx_t x, int size, unsigned int ui) {
	x[0] = ui;
	if (size > 1)
		memset(&x[1], 0, (size - 1) * sizeof(mp_limb_t));
}

/*
 * Set a fixed-size integer from an mpz_t.
 * Inputs:
 *   mpx_t x: allocated fixed-size integer to write to
 *   int size: size of x (in mp_limb_t)
 *   mpz_t z: the value to set
 * Returns:
 *   Nothing.
 */
void mpx_set_z(mpx_t x, int size, mpz_t z) {
	int zsize = z->_mp_size;
	if (zsize)
		memcpy(&x[0], &z->_mp_d[0], zsize * sizeof(mp_limb_t));
	if (zsize < size)
		memset(&x[zsize], 0, (size - zsize) * sizeof(mp_limb_t));
}

/*
 * Set a fixed-size integer from another fixed-size integer
 * Inputs:
 *   mpx_t xd: fixed-size integer to write to
 *   int sized: size of fixed-size integer to write to
 *   mpx_t xs: fixed-size integer to read from
 *   int sizes: size of fixed-size integer to read from
 * Returns:
 *   Nothing.
 */
void mpx_set(mpx_t xd, int sized, mpx_t xs, int sizes) {
	int head = (sized > sizes) ? sizes : sized;
	int tail = sized - sizes;
	if (head)
		memcpy(&xd[0], &xs[0], head * sizeof(mp_limb_t));
	if (tail)
		memset(&xd[head], 0, tail * sizeof(mp_limb_t));
}

/*
 * Set a GMP integer from a fixed-size integer
 * Inputs:
 *   mpz_t z: the initialized mpz_t to set
 *   mpx_t x: the fixed-size integer to set from
 *   int size: size of x (in mp_limb_t)
 * Returns:
 *   Nothing.
 */
void mpz_set_x(mpz_t z, mpx_t x, int size) {
	while (size && x[size - 1] == 0)
		--size;
	if (z->_mp_alloc < size)
		mpz_realloc2(z, size * GMP_LIMB_BITS);
	if (size)
		__GMPN_COPY(&z->_mp_d[0], &x[0], size);
	z->_mp_size = size;
}

/*
 * Fixed-size add routines
 */

#define ADD_BODY(text) \
	asm( \
		text \
		: /* no output operands */ \
		: "r" (xd), "r" (xs1), "r" (xs2) \
		: "%rax", "memory" \
	)

#define ADD_FIRST() \
	"movq 0(%1), %%rax; \n\t" \
	"addq 0(%2), %%rax; \n\t" \
	"movq %%rax, 0(%0); \n\t"

#define ADD_NEXT(offset) \
	"movq " #offset "(%1), %%rax; \n\t" \
	"adcq " #offset "(%2), %%rax; \n\t" \
	"movq %%rax, " #offset "(%0); \n\t"

void mpx_add_1(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST()
	);
}

void mpx_add_2(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8)
	);
}

void mpx_add_3(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16)
	);
}

void mpx_add_4(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16) ADD_NEXT(24)
	);
}

void mpx_add_5(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16) ADD_NEXT(24)
		ADD_NEXT(32)
	);
}

void mpx_add_6(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16) ADD_NEXT(24)
		ADD_NEXT(32) ADD_NEXT(40)
	);
}

void mpx_add_7(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16) ADD_NEXT(24)
		ADD_NEXT(32) ADD_NEXT(40) ADD_NEXT(48)
	);
}

void mpx_add_8(mpx_t xd, mpx_t xs1, mpx_t xs2) {
	ADD_BODY(
		ADD_FIRST() ADD_NEXT(8) ADD_NEXT(16) ADD_NEXT(24)
		ADD_NEXT(32) ADD_NEXT(40) ADD_NEXT(48) ADD_NEXT(56)
	);
}

/*
 * Fixed-size compare routines
 */
#define CMP_BODY(text) \
	int result; \
	asm( \
		text \
		: "=&a" (result) \
		: "c" (xs1), "d" (xs2) \
	); \
	return result

#define CMP_CHUNK(size, offset) \
	"movq " #offset "(%1), %%rax; \n\t" \
	"cmpq " #offset "(%2), %%rax; \n\t" \
	"jc fixed_cmp_" #size "_lt; \n\t" \
	"jne fixed_cmp_" #size "_gt; \n\t"

#define CMP_TAIL(size) \
	"movl $0, %0; \n\t" \
	"jmp fixed_cmp_" #size "_return; \n\t" \
  "fixed_cmp_" #size "_lt: \n\t" \
	"movl $-1, %0; \n\t" \
	"jmp fixed_cmp_" #size "_return; \n\t" \
  "fixed_cmp_" #size "_gt: \n\t" \
	"movl $1, %0; \n\t" \
  "fixed_cmp_" #size "_return: \n\t"

int mpx_cmp_1(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(1, 0)
		CMP_TAIL(1)
	);
}

int mpx_cmp_2(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(2, 8) CMP_CHUNK(2, 0)
		CMP_TAIL(2)
	);
}

int mpx_cmp_3(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(3, 16) CMP_CHUNK(3, 8) CMP_CHUNK(3, 0)
		CMP_TAIL(3)
	);
}

int mpx_cmp_4(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(4, 24) CMP_CHUNK(4, 16) CMP_CHUNK(4, 8) CMP_CHUNK(4, 0)
		CMP_TAIL(4)
	);
}

int mpx_cmp_5(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(5, 32)
		CMP_CHUNK(5, 24) CMP_CHUNK(5, 16) CMP_CHUNK(5, 8) CMP_CHUNK(5, 0)
		CMP_TAIL(5)
	);
}

int mpx_cmp_6(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(6, 40) CMP_CHUNK(6, 32)
		CMP_CHUNK(6, 24) CMP_CHUNK(6, 16) CMP_CHUNK(6, 8) CMP_CHUNK(6, 0)
		CMP_TAIL(6)
	);
}

int mpx_cmp_7(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(7, 48) CMP_CHUNK(7, 40) CMP_CHUNK(7, 32)
		CMP_CHUNK(7, 24) CMP_CHUNK(7, 16) CMP_CHUNK(7, 8) CMP_CHUNK(7, 0)
		CMP_TAIL(7)
	);
}

int mpx_cmp_8(mpx_t xs1, mpx_t xs2) {
	CMP_BODY(
		CMP_CHUNK(8, 56) CMP_CHUNK(8, 48) CMP_CHUNK(8, 40) CMP_CHUNK(8, 32)
		CMP_CHUNK(8, 24) CMP_CHUNK(8, 16) CMP_CHUNK(8, 8) CMP_CHUNK(8, 0)
		CMP_TAIL(8)
	);
}

mpx_support xsupport[MPX_MAXLIMBS] = {
	{ 1, &mpx_add_1, &mpx_cmp_1 },
	{ 2, &mpx_add_2, &mpx_cmp_2 },
	{ 3, &mpx_add_3, &mpx_cmp_3 },
	{ 4, &mpx_add_4, &mpx_cmp_4 },
	{ 5, &mpx_add_5, &mpx_cmp_5 },
	{ 6, &mpx_add_6, &mpx_cmp_6 },
	{ 7, &mpx_add_7, &mpx_cmp_7 },
	{ 8, &mpx_add_8, &mpx_cmp_8 }
};

mpx_support* mpx_support_n(int n) {
	if (n <= 0 || n > MPX_MAXLIMBS) {
		fprintf(stderr, "Error: no mpx support available for size %d\n", n);
		exit(1);
	}
	return &xsupport[n - 1];
}

mpx_support* mpx_support_z(mpz_t z) {
	return mpx_support_n(z->_mp_size);
}
