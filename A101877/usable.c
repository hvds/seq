#include <stdlib.h>
#include "usable.h"
#include "prime.h"

int* usable = (int*)NULL;
int usable_size = 0;
int usable_limit;

void usable_grow(int size) {
	usable_size = size;
	usable = realloc((void*)usable, usable_size * sizeof(int));
}

void setup_usable(int k) {
	int hnum = 0, hden = 1, rnum, rden, shared;
	int index = 0;
	usable_grow(10);
	while (hden < k) {
		if (index >= usable_size)
			usable_grow(usable_size * 3 / 2);
		usable[index++] = hnum;

		rnum = 1;
		rden = index;
		shared = gcd(hden, rden);
		hden /= shared;
		rden /= shared;
		hnum = hnum * rden + rnum * hden;
		hden = hden * rden * shared;
	}
	usable[index] = hnum;
	usable_limit = index;
}

void teardown_usable(void) {
	free(usable);
	usable = (int*)NULL;
	usable_size = 0;
}

int is_usable(int p, int pp, int k) {
	int count = k / pp;
	if (count > usable_limit)
		return 1;
	if (usable[count] >= p)
		return 1;
	return 0;
}
