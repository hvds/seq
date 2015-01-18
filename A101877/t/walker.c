#include "walker.h"
#include "pp.h"
#include <stdio.h>
#include <stdlib.h>

int g_fail = 0;
int g_test = 0;

void test_empty(walk_result* wr) {
	if (wr) {
		printf("Error: expected end of walker iterator, got <%d %d %d>\n",
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
	} else if (wr->discard != discard
			|| wr->invsum != invsum
			|| wr->vec[0] != vec0) {
		printf("Error: expected walker iteration <%d %d %d>, got <%d %d %d>\n",
			discard, invsum, vec0, wr->discard, wr->invsum, wr->vec[0]);
		++g_fail;
	} else {
		printf("Ok: walker iteration <%d %d %d>\n", discard, invsum, vec0);
	}
	++g_test;
}

int main(int argc, char** argv) {
	int i, j;
	pp_pp pp;
	walker *w, *w2;
	walk_result* wr;

	setup_walker();
	++g_test;
	pp.p = 3;
	pp.pp = 27;
	pp.valsize = 3;
	pp.value = calloc(10, sizeof(pp_value));
	pp.value[0].value = 10;
	pp.value[0].inv = 1;
	pp.value[1].value = 7;
	pp.value[1].inv = 1;
	pp.value[2].value = 5;
	pp.value[2].inv = 2;
	w = new_walker(&pp, 0, -1);
	++g_test;
	test_wr(walker_next(w), 0, 0, 0);
	test_wr(walker_next(w), 5, 2, 4);
	test_wr(walker_next(w), 7, 1, 2);
	test_wr(walker_next(w), 10, 1, 1);
	w2 = new_walker(&pp, 13, 0);
	test_wr(walker_find(w2), 0, 0, 0);
	test_wr(walker_find(w2), 12, 0, 6);
	test_empty(walker_find(w2));
	delete_walker(w2);
	test_wr(walker_next(w), 12, 0, 6);
	test_wr(walker_next(w), 15, 0, 5);
	test_wr(walker_next(w), 17, 2, 3);
	test_wr(walker_next(w), 22, 1, 7);
	test_empty(walker_next(w));
	delete_walker(w);
	teardown_walker();
	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
