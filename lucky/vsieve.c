#include <stdlib.h>
#include <stdio.h>
typedef unsigned long long int ullong;

//========================================================================
// This function generates a b-file for the lucky numbers (A000959)
// Run "lucky <elems>" to print elements 1 through <elems>
// It uses a "virtual sieve" and requires O(<elems>) space
// It runs as fast of faster than an explicit sieving algorithm
int main(int argc, char** argv) {
	ullong elems = argc <= 1 ? 1 : (ullong)atol(argv[1]);
	ullong* lucky = (ullong*)malloc(elems * sizeof(ullong));
	ullong g, k, i, n;

	if (lucky == NULL) {
		fprintf(stderr, "Too many elements requested\n");
		return(1);
	}

	lucky[0] = 1;
	lucky[1] = 3;
	if (elems >= 1) printf("1\n");
	if (elems >= 2) printf("3\n");

	// g is the largest index with lucky[g] <= n+1
	for (n = 2, g = 0; n < elems; ++n) {

		// Update g to largest index with lucky[g] <= n+1
		if (lucky[g + 1] <= n + 1) ++g;

		// Now we are going to trace the position k of the nth
		// lucky number backwards through the sieving process.
		// k is the nth lucky number, so it is at position n
		// after all the sieves.
		k = n;

		// If lucky[i] > n+1, the sieve on lucky[i] does not alter
		// the position of the nth lucky number, that is, does not
		// alter k. So we need to run backwards through the sieves
		// for which lucky[i] <= n+1. The last such sieve is the
		// sieve for lucky[g], by definition of g.

		// So, we run backwards through the sieves for lucky[g]
		// down to the sieve for lucky[1] = 3.
		for (i = g; i > 0; i--) {
			// Here k is the position of the nth lucky number
			// after the sieve on lucky[i]. Adjust the position
			// prior to the sieve on lucky[i].
			k = k * lucky[i] / (lucky[i] - 1);
		}

		// Here k is the position of the nth lucky number prior to
		// sieve on 3, that is, after the sieve on 2. Adjust the
		// position prior to the sieve on 2.
		k = 2 * k;

		// Here k is the position of the nth lucky number prior to
		// the sieve on 2, that is, within the natural numbers
		// (1, 2, 3, ...) indexed on 0. So the nth lucky number is
		lucky[n] = k + 1;

		// Adjust n for 1-indexing and print our new value
		printf("%u\n", k + 1);
	}

	free(lucky);
	return 0;
}
