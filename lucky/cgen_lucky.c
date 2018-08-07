#include <stdlib.h>
#include <stdio.h>

typedef unsigned int uint;
typedef struct s_lucky {
	uint lucky;	/* the lucky number */
	uint count; /* remaining to pass through until next discard */
} t_lucky;

t_lucky* plucky;	/* array of lucky numbers and skip counts */
uint lucky_size;	/* current size of the array */
uint lucky_count;	/* index at which to store next new lucky number */
uint lucky_index;	/* index of first lucky number greater than current count */
uint target;		/* find lucky numbers up to this limit */

/* resize the plucky[] array to make room for more numbers */
void resize_plucky(void) {
	uint new_size = lucky_size * 3 / 2;
	plucky = realloc(plucky, new_size * sizeof(t_lucky));
	lucky_size = new_size;
}

/* initialize globals */
void init(void) {
	plucky = (t_lucky*)NULL;
	lucky_size = 1024;
	resize_plucky();
	lucky_count = 1;
	lucky_index = 1;
	plucky[0].lucky = 1;
	printf("%u\n", (uint)1);
}

/* find and return the next lucky number, or 0 if we've reached the target */
uint next_lucky(uint current) {
	uint i;
	int skip;

	while (current < target) {
		current += 2;
		skip = 0;
		/* we only maintain and check counts for those lucky numbers that
		 * have had at least one discard, ie indices 1 .. lucky_index
		 */
		for (i = 1; i < lucky_index; ++i) {
			if (--plucky[i].count == 0) {
				/* reset the count, and discard this candidate */
				plucky[i].count = plucky[i].lucky;
				skip = 1;
				break;
			}
		}
		if (skip) continue;	/* has been discarded */

		/* assume this is a lucky number */
		if (lucky_size == lucky_count) resize_plucky();
		plucky[lucky_count].lucky = current;
		plucky[lucky_count].count = current;
		++lucky_count;

		/* having made sure that plucky[lucky_index] exists, check if this
		 * is its first discard: if so, discard it.
		 */
		if (plucky[lucky_index].lucky == lucky_count) {
			++lucky_index;
			--lucky_count;
			continue;
		}

		/* we have survived all bullets - lucky us */
		return current;
	}
	/* we've exceeded the target */
	return 0;
}

int main(int argc, char** argv) {
	uint current = 1;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <target>\n", argv[0]);
		exit(1);
	}
	target = (uint)atol(argv[1]);
	init();
	while (current = next_lucky(current)) {
		printf("%u\n", current);
	}
	exit(0);
}

/* timings:
  sudo /usr/bin/time taskset -c 0 ./cgen_lucky <n> >/dev/null
1e6: 0.52s
5e6: 8.42s
1e7: 28.84s
2e7: 99.45s
4e7: 346.30s

This suggests a power law: t ~= cn^a, with c=8.55e-12, a=1.79, predicting:
1e6: 0.47s
1e7: 28.96s
1e8: 1785.72s (0.5h)
1e9: 110106.80s (30.6h)

For n=1e9, I expect memory usage to be about 500MB.

Later: I was over-optimistic. With target 1e9, it reached 3e8 after 17h.
Memory usage at that point was 133MB, so that is still within the range.

*/
