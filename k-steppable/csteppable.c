#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#ifndef ulong
#define ulong unsigned long
#endif

#define DEBUG 0

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
*/
ulong try_match_r(char* s, ulong prev, ulong offset, ulong length, ulong steps);
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

