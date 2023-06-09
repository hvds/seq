/*
  A169858(n) = min(m: forall d in {1..n}: exists k in {0..log_10(m)}:
      d | floor(m / 10^k))
  A177834(n) = min(m: forall d in {1..n}: exists j, k in {0..log_10(m)}: j < k &
      exists s: s = floor(m / 10^k) % 10^{k-j} & s != 0 & d | s)

That is to say:
- A169858(n) is the least value m such that all integers from 1..n divide some
prefix of m (written in base 10).
- A177834(n) is the least value m such that all integers from 1..n divide some
non-zero substring of m (written in base 10).

For the purposes of calculation:
1) we need only check divisibility of d in {ceil(n/2)..n};
2) we cache a bitvector of divisibility of the relevant substrings for each
prefix of a candidate m;
3) after testing a candidate m and finding it is not a solution, we choose
the largest divisor d_0 < n that does not divide the prefix of m, and pick
the next candidate m' = m + d_0 - (m % d_0), noting that this cannot skip
past a candidate for which d_0 divides a prefix.

*/

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef NO_TIMING
# include <sys/times.h>
# include <unistd.h>
#endif

typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct s_prefix {
	uint divisors;	/* location of bit vector of unsatisfied divisors */
	mpz_t value;	/* actual value of this prefix */
	uint digit;		/* new digit compared to previous prefix (== value % 10) */
} t_prefix;

uint n;				/* Number of divisors to test for */
uint mindiv;		/* ceil(n/2), the minimum divisor that needs testing */
mpz_t m;			/* Current candidate number */
uint digits;		/* number of digits in m */
t_prefix* prefix;	/* list of prefixes of m */
uint vecsize;		/* size of the divisor vectors */
uchar* vectors;		/* space for the divisor vectors */
#ifndef NO_TIMING
long clock_tick;	/* used to show timings */
#endif

/* Return pointer to divisors bitvector for prefix <i> */
static inline uchar* VEC(uint i) {
	return vectors + prefix[i].divisors;
}
/* Set bit number <d> in the bitvector for prefix <i> */
static inline void VECSET(uint i, uint d) {
	uint bit = d & 7;
	uint offset = d >> 3;
	VEC(i)[offset] |= 1 << bit;
}
/* Clear bit number <d> in the bitvector for prefix <i> */
static inline void VECCLEAR(uint i, uint d) {
	uint bit = d & 7;
	uint offset = d >> 3;
	VEC(i)[offset] &= ~(1 << bit);
}
/* Test bit number <d> in the bitvector for prefix <i>, return TRUE if set */
static inline int VECTEST(uint i, uint d) {
	uint bit = d & 7;
	uint offset = d >> 3;
	return (VEC(i)[offset] & (1 << bit)) ? 1 : 0;
}

/*
 * Set up the bitvector for prefix <i> to have bits set only for those
 * divisors that do not divide this prefix or any of its prefixes.
 * Requires 1 <= i <= digits.
 */
void test_divisors(uint i) {
	uint d;

	/* Copy my immediate prefix's bitvector */
	memcpy(VEC(i), VEC(i - 1), vecsize);

	/* Test divisibility for those divisors not already satisfied */
	for (d = mindiv; d <= n; ++d) {
		if (!VECTEST(i, d)) continue;
		if (mpz_divisible_ui_p(prefix[i].value, d))
			VECCLEAR(i, d);
	}
}

/*
 * Initialize data structures and helper variables:
 * - mindiv = ceil(n/2)
 * - vecsize = number of bytes in vector needed to reach bit n
 * - digits = number of digits in m
 * - prefix is an array of digits+1 t_prefix structures
 * - vectors is an array of digits+1 bitvectors, each of vecsize bytes
 * When the initialization is complete, we know m is a solution iff the
 * bitvector for prefix[digits] has all its bits clear.
 */
void init_mn(void) {
	uchar* s;
	uint i, d;

	mindiv = 1 + n / 2;
	vecsize = 1 + n / 8;
	s = (uchar*)mpz_get_str((char*)NULL, 10, m);
	digits = strlen(s);
	prefix = malloc((digits + 1) * sizeof(t_prefix));
	vectors = malloc(vecsize * (digits + 1) * sizeof(uchar));
	for (i = 0; i <= digits; ++i) {
		prefix[i].divisors = i * vecsize;
		mpz_init(prefix[i].value);
		prefix[i].digit = i ? (s[i - 1] - '0') : 0;
		if (i == 0) {
			mpz_set_ui(prefix[i].value, 0);
			memset(VEC(i), 0, vecsize);
			for (d = mindiv; d <= n; ++d)
				VECSET(i, d);
		} else {
			mpz_mul_ui(prefix[i].value, prefix[i - 1].value, 10);
			mpz_add_ui(prefix[i].value, prefix[i].value, prefix[i].digit);
			test_divisors(i);
		}
	}
}

/*
 * Reinitialize data structures and helper variables when n changes:
 * - mindiv and vecsize are recalculated
 * - if vecsize has changed, the space is reallocated, with associated fixups
 * - the bitvectors are recalculated to incorporate the new divisors
 */
void reinit_n(void) {
	uint oldsize, i;
	char* oldvec;

	mindiv = 1 + n / 2;
	oldsize = vecsize;
	vecsize = 1 + n / 8;
	if (vecsize != oldsize) {
		vectors = realloc(vectors, vecsize * (digits + 1) * sizeof(uchar));
		for (i = 0; i <= digits; ++i)
			prefix[i].divisors = i * vecsize;
	}
	VECSET(0, n);
	for (i = 1; i <= digits; ++i)
		test_divisors(i);
}

/*
 * Test if m is a solution for a(n): it is a solution if all bits in the
 * bitvector for prefix[digits] are clear, signalling that all divisors
 * are satisfied.
 */
int test_mn(void) {
	uint d;
#ifdef A175516
	/* A175516(n): look only at even m */
	if (prefix[digits].digit & 1)
		return 0;
#endif
	for (d = mindiv; d <= n; ++d) {
		if (VECTEST(digits, d))
			return 0;
	}
	return 1;
}

/*
 * Decide how far to skip forward for the next candidate m.
 * If some divisor d does not divide any of the prefixes of the current m,
 * we cannot have a solution any earlier than the next multiple of d.
 * We choose the largest unsatisfied divisor, and skip to its next multiple.
 *
 * In some cases one could skip further by choosing the maximum gap to the
 * first multiple of *all* the divisors not satisfied by our prefixes, but
 * the cost of calculating that would probably exceed the benefit.
 */
uint skip_size(void) {
	uint d;

	/* Find the largest divisor not satisfied by our immediate prefix */
	for (d = n; 1; --d) {
		if (VECTEST(digits - 1, d) != 0)
			break;
	}
	/* if some prefix x would be a multiple of d, then we will find x0 here,
	   so we don't need to check the prefixes separately */
	return d - mpz_fdiv_ui(m, d);
}

/*
 * Pick a new candidate m, and update data structures.
 */
void skip_to_next_m(void) {
	uint skip = skip_size();
	uint i = digits + 1, j;

	mpz_add_ui(m, m, skip);
	while (skip) {
		--i;
		if (i == 0) {
			/* We've passed a power of 10, so the new m has more digits than
			 * the old. Resize the data structures and fix up. This happens
			 * rarely, so it doesn't matter if it is a little slow.
			 */
			++digits;
			prefix = realloc(prefix, (digits + 1) * sizeof(t_prefix));
			vectors = realloc(vectors, vecsize * (digits + 1) * sizeof(uchar));
			prefix[digits].divisors = digits * vecsize;
			mpz_init(prefix[digits].value);
			for (j = digits; j > 0; --j) {
				prefix[j].digit = prefix[j - 1].digit;
				mpz_set(prefix[j].value, prefix[j - 1].value);
			}
			++i;
		}
		mpz_add_ui(prefix[i].value, prefix[i].value, skip);
		skip += prefix[i].digit;
		prefix[i].digit = skip % 10;
		skip /= 10;
	}
	/* Re-test divisors only for the values that have changed */
	for (; i <= digits; ++i)
		test_divisors(i);
}

#ifdef NO_TIMING
double timing(void) {
	return 0;
}
#else
/* Return the user CPU time in seconds. */
double timing(void) {
	struct tms ttd;
	times(&ttd);
	return ((double)ttd.tms_utime) / clock_tick;
}
#endif

/*
 * Starting from the current <m>, search for a solution to a(n) and print it.
 */
void search_n(void) {
	while (1) {
		if (test_mn()) {
			gmp_printf("a(%u) = %Zu [%.2f]\n", n, m, timing());
			return;
		}
		skip_to_next_m();
	}
}

int main(int argc, char** argv) {
	if (argc < 2 || argc > 3) {
		fprintf(
			stderr, "Usage: %s <n> [<min>]\n"
			"Search for A169858(n) starting at min (default min=1)\n",
			argv[0]
		);
		return 0;
	}
	n = atoi(argv[1]);
	if (argc > 2) {
		mpz_init_set_str(m, argv[2], 10);
	} else {
		mpz_init_set_ui(m, 1);
	}
#ifndef NO_TIMING
	clock_tick = sysconf(_SC_CLK_TCK);
#endif
	init_mn();
	while (1) {
		search_n();
		++n;
		reinit_n();
	}
	return 0;
}
