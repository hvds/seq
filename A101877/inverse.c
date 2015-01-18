/*
 * Support for inverse (mod p).
 */

#include "inverse.h"

/*
 * Extended Euclidean algorithm
 * Input:
 *   unsigned int n, 1 <= n < 2^31
 *   unsigned int m, 1 <= m < 2^31
 *   gcd(n, m) == 1
 * Returns:
 *   unsigned int r, 0 <= r < m, nr == 1 (mod m)
 * Time:
 *   O(log^2(m))
 */
uint inveuclid(uint n, uint m) {
	uint a = m, b = n, q, nexta;
	int x = 1, lastx = 0, nextx;
	while (b != 0) {
		q = a / b;
		nexta = b;
		b = a % b;
		a = nexta;
		nextx = lastx - q * x;
		lastx = x;
		x = nextx;
	}
#ifdef DEBUG
	if ((lastx < -(int)m) || (lastx >= (int)m)) {
		fprintf(stderr, "Value out of range: %u (mod %u) => %d\n", n, m, lastx);
	}
#endif
	return (lastx < 0) ? lastx + m : lastx;
}

/*
 * Table inverse
 * Input:
 *   unsigned int p, 3 <= p < 2^31, p an odd prime
 *   unsigned int t[p], pointer to a mutable array of p unsigned ints
 * Returns:
 *   Nothing, fills in t[k] == k^-1 (mod p), 1 <= k <= p-1
 * Time:
 *   ?
 */
void invtable(uint p, uint* t) {
	uint n = 1;
	uint r = 1;
	uint gap = 1;
	memset((void*)t, 0, p * sizeof(uint));
	while (1) {
		t[n] = r;
		t[r] = n;
		t[p - n] = p - r;
		t[p - r] = p - n;
		n <<= 1;
		if (n > p) n -= p;
		r = ((r & 1) ? (r + p) : r) >> 1;
		if (t[n] == 0) continue;
		while (gap < p && t[gap]) ++gap;
		if (gap == p) break;
		n = gap;
		r = inveuclid(n, p);
	}
}

uint** inverse;
uint inverse_size = 0;

/*
 * Initialize for use of inverse caching (see invfast())
 */
void setup_inverse(void) {
	inverse_size = 0;
	inverse = (uint**)NULL;
}

/*
 * Free memory used for inverse caching (see invfast())
 */
void teardown_inverse(void) {
	int i;
	for (i = 0; i < inverse_size; ++i) {
		if (inverse[i]) free(inverse[i]);
	}
	free(inverse);
}

/*
 * Set up cache of inverses (mod p)
 * Input:
 *   unsigned int p, 2 <= p < 2^31, p prime
 *   previous initialization call to setup_inverse()
 *   eventual teardown call to teardown_inverse()
 * Returns:
 *   Nothing, but caches all inverses (mod p) for later calls to invfast()
 * Time:
 *   O(invtable(p)) + O(p)
 */
void inverse_table(uint p) {
	uint size = p + 1;
	uint effective_size = p + ((p <= 2) ? 1 : 0);
	if (inverse_size < size) {
		inverse = realloc(inverse, size * sizeof(uint*));
		memset(&inverse[inverse_size], 0,
				(size - inverse_size) * sizeof(uint*));
		inverse_size = size;
	}
	inverse[p] = calloc(effective_size, sizeof(uint));
	invtable(p, inverse[p]);
}

/*
 * Fast modular inverse
 * Input:
 *   unsigned int n, 0 <= n < p
 *   unsigned int p, 2 <= p < 2^31, p prime
 *   previous call to inverse_table(p)
 * Returns:
 *   unsigned int r, 0 <= r < p, nr == 1 (mod p)
 * Time:
 *   1
 * Caveats:
 *   For speed, no checking is done: constraints above must be satisfied.
 *   All calls are expected to use the inline definition in "inverse.h".
 */
uint invfast(uint n, uint p) {
	return inverse[p][n];
}
