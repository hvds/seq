#include <stdio.h>
#include <stdlib.h>
#include "walker.h"
#include "pp.h"
#include "mygmp.h"

int g_fail = 0;
int g_test = 0;

double timing(void) { return 0; }

void test_empty(whp wh, wrhp wrh) {
	if (wrh) {
		mpz_t z;
		ZINIT(&z, "test_empty temp");
		mpz_set_x(z, wr_discard(wh, wrh), WP(wh)->numsize);
		gmp_printf("Error: expected end of walker iterator, got <%Zd %d %d>\n",
				z, WRP(wh, wrh)->invsum, wr_vec(wh, wrh)[0]);
		ZCLEAR(&z, "test_empty temp");
		++g_fail;
	} else {
		printf("Ok: end of walker iterator\n");
	}
	++g_test;
}

void test_wr(whp wh, wrhp wrh, int discard, int invsum, int vec0) {
	mpz_t z;
	if (!wrh) {
		printf("Error: expected walker iteration <%d %d %d>, got empty\n",
			discard, invsum, vec0);
		++g_fail;
	} else {
		ZINIT(&z, "test_wr temp");
		mpz_set_x(z, wr_discard(wh, wrh), WP(wh)->numsize);
		if (mpz_cmp_ui(z, discard) != 0
				|| WRP(wh, wrh)->invsum != invsum
				|| wr_vec(wh, wrh)[0] != vec0) {
			gmp_printf("Error: expected walk_result <%d %d %d>, got <%Zd %d %d>\n",
				discard, invsum, vec0, z, WRP(wh, wrh)->invsum, wr_vec(wh, wrh)[0]);
			ZCLEAR(&z, "test_wr temp");
			++g_fail;
		} else {
			printf("Ok: walker iteration <%d %d %d>\n", discard, invsum, vec0);
		}
	}
	++g_test;
}

void set_value(pp_pp* pp, int index, int val, int inv) {
	mpx_set_ui(ppv_mpx(pp_value_i(pp, index)), pp->valnumsize, val);
	pp_value_i(pp, index)->inv = inv;
}

int main(int argc, char** argv) {
	int i, j;
	pp_pp pp;
	whp wh, wh2;
	wrhp wrh;
	mpz_t limit;
	mpx_support* xsup;

	ZINIT(&limit, "test limit");
	setup_walker();
	++g_test;
	pp.p = 3;
	pp.pp = 27;
	pp.valsize = 3;
	pp.valnumsize = 1;
	pp.value = calloc(10, pp_valsize_n(3));
	xsup = mpx_support_n(1);
	pp.adder = xsup->adder;
	pp.cmper = xsup->cmper;
	set_value(&pp, 0, 10, 1);
	set_value(&pp, 1, 7, 1);
	set_value(&pp, 2, 5, 2);
	mpz_set_ui(limit, 22);

	wh = new_walker(&pp, limit, -1);
	++g_test;
	test_wr(wh, walker_findnext(wh), 0, 0, 0);
	test_wr(wh, walker_findnext(wh), 5, 2, 4);
	test_wr(wh, walker_findnext(wh), 7, 1, 2);
	test_wr(wh, walker_findnext(wh), 10, 1, 1);
	mpz_set_ui(limit, 13);
	wh2 = new_walker(&pp, limit, 0);
	test_wr(wh2, walker_findnext(wh2), 0, 0, 0);
	test_wr(wh2, walker_findnext(wh2), 12, 0, 6);
	test_empty(wh2, walker_findnext(wh2));
	delete_walker(wh2);
	test_wr(wh, walker_findnext(wh), 12, 0, 6);
	test_wr(wh, walker_findnext(wh), 15, 0, 5);
	test_wr(wh, walker_findnext(wh), 17, 2, 3);
	test_wr(wh, walker_findnext(wh), 22, 1, 7);
	test_empty(wh, walker_findnext(wh));
	delete_walker(wh);
	teardown_walker();
	free(pp.value);
	ZCLEAR(&limit, "test limit");
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
