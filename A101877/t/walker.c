#include <stdio.h>
#include <stdlib.h>
#include "walker.h"
#include "pp.h"
#include "mygmp.h"

int g_fail = 0;
int g_test = 0;
mpz_t previous;

void test_empty(walk_result* wr) {
	if (wr) {
		gmp_printf("Error: expected end of walker iterator, got <%Zd %d %d>\n",
				wr->discard, wr->invsum, wr->vec[0]);
		++g_fail;
	} else {
		printf("Ok: end of walker iterator\n");
	}
	++g_test;
}

void test_wr(walk_result* wr, int discard, int invsum, int vec0) {
	if (!wr) {
		printf("Error: expected walker iteration <%d %d %d>, got empty\n",
			discard, invsum, vec0);
		++g_fail;
	} else if (mpz_cmp_ui(wr->discard, discard) != 0
			|| wr->invsum != invsum
			|| wr->vec[0] != vec0) {
		gmp_printf("Error: expected walker iteration <%d %d %d>, got <%Zd %d %d>\n",
			discard, invsum, vec0, wr->discard, wr->invsum, wr->vec[0]);
		++g_fail;
	} else {
		printf("Ok: walker iteration <%d %d %d>\n", discard, invsum, vec0);
	}
	++g_test;
}

void set_value(pp_value* v, int val, int inv) {
	ZINIT(&v->value, "test value");
	mpz_set_ui(v->value, val);
	v->inv = inv;
}

int main(int argc, char** argv) {
	int i, j;
	pp_pp pp;
	walker *w, *w2;
	walk_result* wr;
	mpz_t limit;

	ZINIT(&limit, "test limit");
	ZINIT(&previous, "test previous");
	setup_walker();
	++g_test;
	pp.p = 3;
	pp.pp = 27;
	pp.valsize = 3;
	pp.value = calloc(10, sizeof(pp_value));
	set_value(&pp.value[0], 10, 1);
	set_value(&pp.value[1], 7, 1);
	set_value(&pp.value[2], 5, 2);
	mpz_set_ui(limit, 22);
	w = new_walker(&pp, limit, -1);
	++g_test;
	test_wr(walker_next(w), 0, 0, 0);
	test_wr(walker_next(w), 5, 2, 4);
	test_wr(walker_next(w), 7, 1, 2);
	test_wr(walker_next(w), 10, 1, 1);
	mpz_set_ui(limit, 13);
	w2 = new_walker(&pp, limit, 0);
	mpz_set_si(previous, -1);
	test_wr(walker_find(w2, previous), 0, 0, 0);
	test_wr(walker_find(w2, previous), 12, 0, 6);
	test_empty(walker_find(w2, previous));
	delete_walker(w2);
	test_wr(walker_next(w), 12, 0, 6);
	test_wr(walker_next(w), 15, 0, 5);
	test_wr(walker_next(w), 17, 2, 3);
	test_wr(walker_next(w), 22, 1, 7);
	test_empty(walker_next(w));
	delete_walker(w);
	teardown_walker();
	for (i = 0; i < 3; ++i)
		ZCLEAR(&pp.value[i].value, "test value");
	free(pp.value);
	ZCLEAR(&previous, "test previous");
	ZCLEAR(&limit, "test limit");
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
