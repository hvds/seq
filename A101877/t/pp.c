#include "pp.h"
#include <stdio.h>
#include <stdlib.h>

int g_fail = 0;
int g_test = 0;

double timing(void) { return 0; }

void dump_pp(int n) {
	int i, j;
	pp_pp* pp;
	mpz_t zv;

	ZINIT(&zv, "dump_pp zv");
	printf("pp list(%d): ", pplistsize);
	for (i = 0; i < pplistsize; ++i) {
		printf("%p(%d), ", pplist[i], pplist[i]->pp);
	}
	printf("\n");
	for (i = 1; i <= n; ++i) {
		pp = &pppp[i];
		if (pp->p) {
			printf("pppp[%d] = { p = %d; pp = %d; depend = %d; valsize = %d; valmax = %d; value = {",
					i, pp->p, pp->pp, pp->depend, pp->valsize, pp->valmax);
			for (j = 0; j < pp->valsize; ++j) {
				pp_value* v = VALUE_I(pp, j);
				mpz_set_x(zv, MPX(v), pp->valnumsize);
				gmp_printf("{ value = %Zd; parent = %d; inv = %d }",
						zv, v->parent, v->inv);
			}
			gmp_printf("}; min_discard = %Zd; invtotal = %d; total = %Zd; denominator = %Zd; invdenom = %d }\n",
					pp->min_discard, pp->invtotal, pp->total, pp->denominator, pp->invdenom);
		}
	}
	printf("\n");
	ZCLEAR(&zv, "dump_pp zv");
}

int main(int argc, char** argv) {
	int i, j;
	pp_pp* pp;

	setup_pp(24);
	dump_pp(24);
	pp_study(3);
	dump_pp(24);

	if (g_fail) {
		printf("FAIL: failed %u of %u tests.\n", g_fail, g_test);
	} else {
		printf("PASS: passed %u tests.\n", g_test);
	}
	return 0;
}
