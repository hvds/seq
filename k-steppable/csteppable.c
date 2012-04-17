#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

/*

How many of the b^n strings of length n, base b can be partitioned in at least
one way into k non-empty substrings such that the substrings appear in lexical
order?

Input is the 3 arguments (b, k, n); output is (b, k, n, good, bad, time)
where 'good' is the number of strings that can be partitioned in at least
one way as described above; 'bad' is the number that cannot be so
partitioned (and so good + bad = b^n); 'time' is the CPU time taken for
the calculation.

If b^n >= 2^(8 * sizeof(unsigned long)), the code will need to be recompiled
with a larger type, or converted to use bigints.

*/

#ifndef ulong
#define ulong unsigned long
#endif

#define DEBUG 0

/* Context structure holds all the information we need; this was created
   on the assumption that some optimization would involve recursive calls
   to steppable(), but the code doesn't currently do that.
*/
typedef struct ctx_s {
	ulong b;
	ulong k;
	ulong n;
	ulong good;
	ulong bad;
	ulong lim;
	char s[0];
} ctx_t;

long CLK_TCK;

double timer(void) {
	struct tms tbuf;
	times(&tbuf);
	return (double)tbuf.tms_utime / CLK_TCK;
}

/* Return b^n; a mild attempt is made to detect overflow.
*/
ulong upow(ulong b, ulong n) {
	ulong p = 1, q, i = 0;
	while (i < n) {
		q = p * b;
		if (q < p) {
			fprintf(stderr, "overflow: cannot count up to b^n\n");
			exit(-1);
		}
		++i;
		p = q;
	}
	return p;
}

/* Advance the string in ctx->s to the next base-b value, ignoring the
   last <offset-1> characters.

   Returns the offset of the most significant changed character, or 0
   if the advancement wrapped the whole string back to zero.
*/
ulong advance(ctx_t *ctx, ulong offset) {
	int i;
	ulong lim = ctx->b;
	for (i = ctx->n - offset; i >= 0; --i) {
		if (++ctx->s[i] != lim) {
			return ctx->n - i;
		}
		ctx->s[i] = 0;
	}
	return 0;
}

/*
  In the string <s>, attempt to find a further <steps> substrings, such that
  the next substring starts at <offset> or greater, and is lexically greater
  than the preceding substring starting at <prev>.

  If successful, returns the offset at which success was achieved (ie, the
  number of characters used for the full <k> substrings), else returns 0.
*/
ulong try_match_r(char* s, ulong prev, ulong offset, ulong length, ulong steps) {
	ulong thislen, result;

	if (steps < 1) return offset;
	for (; offset < length; ++offset) {
		ulong prevlen = offset - prev;
		ulong maxthislen = length - offset;
		int equal = 1;
		for (thislen = 1; thislen <= maxthislen; ++thislen) {
			if (equal) {
				if (thislen > prevlen) {
					/* <prev> is now strictly a prefix of <this> */
					equal = 0;
				} else {
					signed char c = s[prev + thislen - 1] - s[offset + thislen - 1];
					if (c < (signed char)0) {
						/* <prev> is strictly less than <this> */
						equal = 0;
					} else if (c > (signed char)0) {
						/* <prev> is strictly greater than <this> */
						break;
					} else {
						/* <prev> is same as <this> to this length */
						continue;
					}
				}
			}

			result = try_match_r(s, offset, offset + thislen, length, steps - 1);
			if (result) return result;
		}
	}
	return 0;
}

/*
  Returns a non-zero value if the string ctx->s can be partitioned into ctx->k
  ordered substrings, else 0.

  The non-zero value is actually the number of characters used for the
  substrings, but we don't currently take advantage of this.
*/
ulong try_match(ctx_t *ctx) {
	return try_match_r(ctx->s, 0, 1, ctx->n, ctx->k - 1);
}

void disp (char* s, ulong n, char* msg) {
	ulong i;
	printf("%s: ", msg);
	for (i = 0; i < n; ++i) {
		printf("%c", '0' + s[i]);
	}
	printf("\n");
}

/*
  Calculates the number of base-b length-n substrings that can and cannot
  be partitioned in at least on way into k substrings that are lexically
  ordered, and returns a filled-in context object with the details.

  It is the caller's responsibility to free the context object.
*/
ctx_t *steppable(ulong b, ulong k, ulong n) {
	int zeroes, match_to;
	ctx_t *ctx = (ctx_t *)malloc(sizeof(ctx_t) + n);

	ctx->b = b;
	ctx->k = k;
	ctx->n = n;
	ctx->bad = 0;
	ctx->lim = upow(b, n);	/* verify b^m fits in sizeof(ulong) */

	/* short-circuit for n <= 1, so we can later assert match_to > 0
	   for every match
	*/
	if (n <= 1) {
		ctx->good = ctx->lim;
		return ctx;
	}

	memset(ctx->s, 0, n);
	/* Except for the initial string 0^n, if we find that the string with
	   <z> trailing zeroes is good, all the strings with that prefix must
	   be good.
	   There would probably be a small further optimization to be gained
	   by taking advantage of the return value of try_match.
	*/
	zeroes = 1;
	while (zeroes) {
		if (!try_match(ctx)) {
			if (DEBUG) disp(ctx->s, n, "bad");
			++ctx->bad;
			zeroes = 1;
		}
		zeroes = advance(ctx, zeroes);
	}
	ctx->good = ctx->lim - ctx->bad;
	return ctx;
}

int main(int argc, char** argv) {
	ulong b, k, n;
	double t0;
	ctx_t *ctx;
	if (argc != 4) {
		fprintf(stderr, "Usage: %s <b> <k> <n>\n", argv[0]);
		return -1;
	}
	CLK_TCK = sysconf(_SC_CLK_TCK);
	b = strtoul(argv[1], NULL, 10);
	k = strtoul(argv[2], NULL, 10);
	n = strtoul(argv[3], NULL, 10);
	t0 = timer();
	ctx = steppable(b, k, n);
	printf("%lu %lu %lu %lu %lu (%.2f)\n",
			b, k, n, ctx->good, ctx->bad, timer() - t0);
	free(ctx);
	return 0;
}

