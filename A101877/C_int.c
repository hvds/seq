#include <stdio.h>
#include <stdlib.h>
#include "pp.h"

int main(int argc, char** argv) {
	int n = 0, k = 0;
	int i, j;

	if (argc > 1)
		n = atoi(argv[1]);
	if (argc > 2)
		k = atoi(argv[2]);
	for (i = n ? n : 1; n ? (i <= n) : 1; ++i) {
		for (j = k ? k : 1; k ? (j <= k) : 1; ++j) {
			init_pp(j);
			pp_study();
			pp_find(i);
			/* close_pp(); */
		}
	}
	return 0;
}
