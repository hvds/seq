#include <stdio.h>
#include <stdlib.h>
#include "walker.h"
#include "pp.h"
#include "mygmp.h"

int g_fail = 0;
int g_test = 0;

void test_empty(walker* w, walk_result* wr) {
	if (wr) {
		mpz_t z;
		ZINIT(&z, "test_empty temp");
		mpz_set_x(z, DISCARD(w, wr), w->numsize);
		gmp_printf("Error: expected end of walker iterator, got <%Zd %d %d>\n",
				z, wr->invsum, VEC(w, wr)[0]);
		ZCLEAR(&z, "test_empty temp");
		++g_fail;
	} else {
		printf("Ok: end of walker iterator\n");
	}
	++g_test;
}

void test_wr(walker* w, walk_result* wr, int discard, int invsum, int vec0) {
	mpz_t z;
	if (!wr) {
		printf("Error: expected walker iteration <%d %d %d>, got empty\n",
			discard, invsum, vec0);
		++g_fail;
	} else {
		ZINIT(&z, "test_wr temp");
		mpz_set_x(z, DISCARD(w, wr), w->numsize);
		if (mpz_cmp_ui(z, discard) != 0
				|| wr->invsum != invsum
				|| VEC(w, wr)[0] != vec0) {
			gmp_printf("Error: expected walk_result <%d %d %d>, got <%Zd %d %d>\n",
				discard, invsum, vec0, z, wr->invsum, VEC(w, wr)[0]);
			ZCLEAR(&z, "test_wr temp");
			++g_fail;
		} else {
			printf("Ok: walker iteration <%d %d %d>\n", discard, invsum, vec0);
		}
	}
	++g_test;
}

void set_value(pp_pp* pp, int index, int val, int inv) {
	mpx_set_ui(MPX(VALUE_I(pp, index)), pp->valnumsize, val);
	VALUE_I(pp, index)->inv = inv;
}

int main(int argc, char** argv) {
	int i, j;
	pp_pp pp;
	walker *w, *w2;
	walk_result* wr;
	mpz_t limit;
	mpx_support* xsup;

	ZINIT(&limit, "test limit");
	setup_walker();
	++g_test;
	pp.p = 3;
	pp.pp = 27;
	pp.valsize = 3;
	pp.valnumsize = 1;
	pp.value = calloc(10, VALSIZE_N(3));
	xsup = mpx_support_n(1);
	pp.adder = xsup->adder;
	pp.cmper = xsup->cmper;
	set_value(&pp, 0, 10, 1);
	set_value(&pp, 1, 7, 1);
	set_value(&pp, 2, 5, 2);
	mpz_set_ui(limit, 22);
	w = new_walker(&pp, limit, -1);
	++g_test;
	test_wr(w, walker_next(w), 0, 0, 0);
	test_wr(w, walker_next(w), 5, 2, 4);
	test_wr(w, walker_next(w), 7, 1, 2);
	test_wr(w, walker_next(w), 10, 1, 1);
	mpz_set_ui(limit, 13);
	w2 = new_walker(&pp, limit, 0);
	test_wr(w2, walker_findnext(w2), 0, 0, 0);
	test_wr(w2, walker_findnext(w2), 12, 0, 6);
	test_empty(w2, walker_findnext(w2));
	delete_walker(w2);
	test_wr(w, walker_next(w), 12, 0, 6);
	test_wr(w, walker_next(w), 15, 0, 5);
	test_wr(w, walker_next(w), 17, 2, 3);
	test_wr(w, walker_next(w), 22, 1, 7);
	test_empty(w, walker_next(w));
	delete_walker(w);
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
