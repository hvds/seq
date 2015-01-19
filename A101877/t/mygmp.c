#include <stdio.h>
#include <stdlib.h>

#include "mygmp.h"

uint g_fail = 0;
uint g_test = 0;
mpz_t z1, z2, z3;
mp_limb_t limbs[100];

#define GUARD 0xa5a5a5a5a5a5a5a5

void guard_mpx_set_z(mpx_t x, int size, mpz_t z) {
	x[-1] = GUARD;
	x[size] = GUARD;
	mpx_set_z(x, size, z);
	if (x[-1] != GUARD || x[size] != GUARD) {
		printf("%d error: mpx_set_z() wrote outside area (x[-1]=%#lx, x[%d]=%#lx)\n",
				g_test, x[-1], size, x[size]);
		++g_fail;
	}
	++g_test;
}

void guard_mpz_set_x(mpz_t z, mpx_t x, int size) {
	x[-1] = GUARD;
	x[size] = GUARD;
	mpz_set_x(z, x, size);
	if (x[-1] != GUARD || x[size] != GUARD) {
		printf("%d error: mpz_set_x() wrote outside area (x[-1]=%#lx, x[%d]=%#lx)\n",
				g_test, x[-1], size, x[size]);
		++g_fail;
	}
	++g_test;
}

void guard_mpx_add(mpx_support* xsup, mpx_t xd, mpx_t xs1, mpx_t xs2) {
	int size = xsup->size;
	xd[-1] = GUARD;
	xd[size] = GUARD;
	xsup->adder(xd, xs1, xs2);
	if (xd[-1] != GUARD || xd[size] != GUARD) {
		printf("Error: mpx adder wrote outside area (x[-1]=%#lx, x[%d]=%#lx)\n",
				xd[-1], size, xd[size]);
		++g_fail;
	}
	++g_test;
}

void test_zero(mpx_t x, int size, char* legend) {
	int i;
	for (i = 0; i < size; ++i) {
		if (x[i] != 0) {
			printf("Error: %s[%d] (size %d) is %#lx, expected zero\n",
					legend, i, size, x[i]);
			++g_fail;
		}
		++g_test;
	}
}

void test_mpx(int size) {
	mpx_support* xsup;
	mpx_t x1, x2, x3;
	int i, result;

	mpz_set_ui(z1, ((uint)1) << (size * 2 + 1));
	for (i = 1; i < size; ++i) {
		mpz_mul_ui(z1, z1, ((uint)1) << 31);
		mpz_mul_ui(z1, z1, ((uint)1) << 31);
    }
	xsup = mpx_support_z(z1);
	if (!xsup || xsup->size != size) {
		if (!xsup)
			gmp_printf("Error: mpx_support_z(size=%d) returned NULL\n", size);
		else
			gmp_printf("Error: mpx_support_z(size=%d) returned xsup[size=%d]\n",
					size, xsup->size);
		++g_fail;
	}
	++g_test;

	x1 = &limbs[1];
	x2 = &limbs[11];
	x3 = &limbs[21];
	mpz_set_ui(z2, 0);
	guard_mpx_set_z(x1, size, z2);
	test_zero(x1, size, "x1");
	guard_mpx_set_z(x2, size, z2);
	test_zero(x2, size, "x2");
	guard_mpx_add(xsup, x3, x1, x2);
	test_zero(x3, size, "x1+x2=x3");
	if ((result = xsup->cmper(x1, x3)) != 0) {
		gmp_printf("Error: 0 cmp 0+0 (size %d) gave %d\n", size, result);
		++g_fail;
	}
	++g_test;
	guard_mpz_set_x(z3, x3, size);
	if (mpz_cmp_ui(z3, 0) != 0) {
		gmp_printf("Error: 0 + 0 (size %d) gave %Zd\n", size, z3);
		++g_fail;
	}
	++g_test;
}

int main(int argc, char** argv) {
	int i;
	ZINIT(&z1, "test z1");
	ZINIT(&z2, "test z2");
	ZINIT(&z3, "test z3");

	for (i = 1; i <= 8; ++i)
		test_mpx(i);

	ZCLEAR(&z3, "test z3");
	ZCLEAR(&z2, "test z2");
	ZCLEAR(&z1, "test z1");
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}

