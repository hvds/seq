#include "part.h"
#include "vec.h"
#include "set.h"
#include "symmetries.h"
#include "pieces.h"
#include "clock.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>

#define REPORT_MASK ((1 << 24) - 1)
#define STR_EVALUATE(x) #x
#define STRINGIFY(x) STR_EVALUATE(x)

counter sym_result = 0;
counter all_result = 0;

typedef struct step_s {
	int parent;		/* step index of last preceding piece of different size */
	int piece_size;	/* size of this piece */
	int remain;		/* number of bits set in freevec */

	set_t set;		/* the combined set of pieces so far */
	vec_t freevec;	/* bits still free */
	vec_t shape;	/* a rotated/reflected piece to insert */
	sym_set_t* ss;	/* symmetries that map the set so far to itself */
} step_t;
step_t steps[NODES + 1];
set_t* solution;
sym_set_t* sym_filters[NODES];

void setup_steps(void) {
	uint i;
	step_t* step;
	sym_set_t* full_ss;

	solution = (set_t*)malloc(sizeof(set_t));
	for (i = 0; i < NODES; ++i)
		steps[i].ss = symset_new();

	step = &(steps[0]);
	vec_zero(&(step->freevec));
	vec_not(&(step->freevec));
	step->remain = NODES;
	step->parent = 0;
	step->piece_size = NODES + 1;
	vec_zero(&(step->shape));
	set_init(&(step->set), &(step->shape));

	full_ss = steps[0].ss;
	for (i = 0; i < sym_count; ++i)
		full_ss->index[i] = i;
	full_ss->count = sym_count;
}

void teardown_steps(void) {
	uint i;
	for (i = 0; i < NODES; ++i)
		symset_delete(steps[i].ss);
	free(solution);
}

/*
  Given a set I<set>, returns C<TRUE> if it is in canonical form (ie if in
  its current state it sorts lexically before any symmetry of itself), else
  C<FALSE>. The set of symmetries I<ssfrom> are the only ones considered.
  If canonical, I<ssto> is also set to the subset of the symmetries in
  I<ssfrom> that map I<set> to itself.

  Note we must compare pieces in strict order (rather than just comparing
  the full string lexically), to ensure every canonical partial has a
  canonical predecessor.
*/
int canonicalize(
	set_t* set, uint piece_size, uint piece_count,
	sym_set_t* ssfrom, sym_set_t* ssto
) {
	uint i, sym_i, ssend;
	uint seen, mapping[NODES], source, dest;
	sym_t* sym;

	ssend = ssfrom->count;
	ssto->count = 0;

	if (piece_count == 1) {
		/* the original faster code is safe for a single piece */

		/* Make a copy with 0 mapped to NODES, so we sort the way we want to */
		for (i = 0; i < NODES; ++i)
			solution->p[i] = set->p[i] ? set->p[i] : NODES;

		for (sym_i = 0; sym_i < ssend; ++sym_i) {
			sym = sym_map(ssfrom->index[sym_i]);
			seen = 0;
			memset(mapping, 0, sizeof(mapping));
			for (i = 0; i < NODES; ++i) {
				source = set->p[sym->map[i]];
				dest = source
					? mapping[source - 1]
						? mapping[source - 1]
						: (mapping[source - 1] = ++seen)
					: NODES;
				/* if this map sorts earlier, the original is not canonical */
				if (dest < solution->p[i])
					return 0;
				/* if this map sorts later, it's fine, but not an identity */
				if (dest > solution->p[i])
					goto NOT_IDENTITY1;
			}
			/* if we fall through it's an identity, so preserve it */
			ssto->index[ssto->count++] = ssfrom->index[sym_i];
		  NOT_IDENTITY1:
			;
		}
	} else {
		uint pivot, broken_index;
		int broken_value;
		uint seen_index[NODES], piece_index[NODES];
		uint cmp, next;

		memset(seen_index, 0, sizeof(seen_index));
		/* Make a copy with 0 mapped to NODES, so we sort the way we want to */
		for (i = 0; i < NODES; ++i) {
			uint offset = set->p[i];
			if (!offset)
				continue;
			--offset;
			piece_index[offset * piece_size + seen_index[offset]++] = i;
		}

		ssend = ssfrom->count;
		ssto->count = 0;
		for (sym_i = 0; sym_i < ssend; ++sym_i) {
			sym = sym_map(ssfrom->index[sym_i]);
			seen = 0;
			memset(mapping, 0, sizeof(mapping));

			pivot = 0;
			broken_value = 0;
			broken_index = piece_count;
			next = piece_index[0];
			memset(seen_index, 0, piece_count * sizeof(uint));
			for (i = 0; i < NODES; ++i) {
				if (i > next)
					goto NOT_IDENTITYn;
				source = set->p[sym->map[i]];
				if (!source)
					continue;
				if (!mapping[source - 1])
					mapping[source - 1] = ++seen;
				dest = mapping[source - 1] - 1;
				if (dest >= broken_index)
					continue;
				cmp = piece_index[dest * piece_size + seen_index[dest]++];
				if (cmp != i) {
					broken_value = (cmp > i) ? -1 : 1;
					broken_index = dest;
					if (pivot < broken_index)
						continue;
					if (broken_value < 0)
						return 0;
					goto NOT_IDENTITYn;
				}

				if (dest > pivot)
					continue;
				if (seen_index[dest] < piece_size) {
					next = piece_index[dest * piece_size + seen_index[dest]];
					continue;
				}
				/* we have exactly matched the pivot piece */
				while (++pivot < piece_count) {
					if (pivot == broken_index) {
						if (broken_value < 0) {
							return 0;
						}
						goto NOT_IDENTITYn;
					}
					if (seen_index[pivot] < piece_size)
						break;
				}
				if (pivot >= piece_count)
					break;
				next = piece_index[pivot * piece_size + seen_index[pivot]];
			}
			/* if we fall through it's an identity, so preserve it */
			ssto->index[ssto->count++] = ssfrom->index[sym_i];
		  NOT_IDENTITYn:
			;
		}
	}
	/* if we get through the gauntlet, we are canonical */
	return 1;
}

void check_solution(uint level) {
	step_t* step = &(steps[level]);
	set_t combined;
	uint i;
	uint syms = sym_count / step->ss->count;

	++sym_result;
	all_result += syms;

#ifdef QUIET
	/* build the printable set only if needed */
	if (! (sym_result & REPORT_MASK)) {
		set_zero(&combined);
		for (i = level; i > 0; i = steps[i].parent)
			set_merge(&(steps[i].set), &combined, steps[i].parent);
		fprintf(stderr, "%llu: (%u) ", sym_result, syms);
		fprint_set(stderr, &combined);
		fprintf(stderr, " (%.2f)\n", GTIME);
	}
#else
	set_zero(&combined);
	for (i = level; i > 0; i = steps[i].parent)
		set_merge(&(steps[i].set), &combined, steps[i].parent);

	fprintf(stdout, "(%u) ", syms);
	fprint_set(stdout, &combined);
	printf("\n");
	if (! (sym_result & REPORT_MASK)) {
		fprintf(stderr, "%llu: (%u) ", sym_result, syms);
		fprint_set(stderr, &combined);
		fprintf(stderr, " (%.2f)\n", GTIME);
	}
#endif
}

inline int insert_piece(uint level, vec_t* piece) {
	step_t* step = &(steps[level]);
	step_t* prev = &(steps[level - 1]);
	uint parent, piece_count;
	if (step->piece_size == prev->piece_size) {
		parent = prev->parent;
		piece_count = level - parent;
		set_append(&(step->set), &(prev->set), piece, piece_count);
	} else {
		piece_count = 1;
		parent = level - 1;
		set_init(&(step->set), piece);
	}
	step->parent = parent;
	return canonicalize(&(step->set), step->piece_size, piece_count,
			steps[parent].ss, step->ss);
}

void try_recurse(uint prev_level) {
	uint level = prev_level + 1;
	step_t* prev = &(steps[prev_level]);
	step_t* step = &(steps[level]);
	uint max_size = prev->remain < prev->piece_size
			? prev->remain : prev->piece_size;
	uint size;

	/* save as a solution the contents of prev_step, assuming all remaining
	 * free locations are filled with pieces of size 1.
	 */
	check_solution(prev_level);

#ifdef REVERSE
	for (size = max_size; size >= 2; --size) {
#else
	for (size = 2; size <= max_size; ++size) {
#endif
		uint pi, end_index;
		step->piece_size = size;
		step->remain = prev->remain - size;
		if (step->remain == 0) {
			/* no point iterating over the pieces, only one can fit */
			if (is_connected(&(prev->freevec)))
				if (insert_piece(level, &(prev->freevec)))
					check_solution(level);
			continue;
		}
		end_index = piece_array[size + 1];
		for (pi = piece_array[size]; pi < end_index; ++pi) {
			piece_t* this_piece = pieces_vec(pi);
			sym_set_t* ss = SSO(this_piece->sso);
			uint sym_c = ss->count;
			uint sym_i;
			for (sym_i = 0; sym_i < sym_c; ++sym_i) {
				apply_map2(
					sym_map(ss->index[sym_i]), &(this_piece->v), &(step->shape)
				);
				if (! vec_contains(&(prev->freevec), &(step->shape)))
					continue;
				vec_xor3(&(prev->freevec), &(step->shape), &(step->freevec));
				if (!insert_piece(level, &(step->shape)))
					continue;
				try_recurse(level);
			}
		}
		if (prev_level == 0) {
			fprintf(stderr, "solutions %u: %llu/%llu (%.2f)\n",
					size, sym_result, all_result, GTIME);
		}
	}
}

void import_pieces(char* filename) {
	struct {
		uint section;
		uint size;
	} section_header;
	int expect = (1 << 0) | (1 << 1);
	int seen = 0;
	FILE* f = fopen(filename, "r");
	int read_ok;

	if (!f) {
		fprintf(stderr, "%s: %s (%d)\n", filename, strerror(errno), errno);
		exit(-1);
	}
	while (1) {
		read_ok = fread(&section_header, sizeof(section_header), 1, f);
		if (read_ok < 1) {
			fprintf(stderr, "%s: %s (%d)\n", filename, strerror(errno), errno);
			fclose(f);
			exit(-1);
		}
		switch (section_header.section) {
		  case 0:
			if (seen & (1 << 0)) {
				fprintf(stderr, "%s: repeated section 0\n", filename);
				fclose(f);
				exit(-1);
			}
			load_sym_sets(f, section_header.size);
			seen = seen | (1 << 0);
			break;
		  case 1:
			if (seen & (1 << 1)) {
				fprintf(stderr, "%s: repeated section 1\n", filename);
				fclose(f);
				exit(-1);
			}
			load_pieces(f, section_header.size);
			seen = seen | (1 << 1);
			break;
		  case 2:
			if (seen == expect) {
				fclose(f);
				return;
			}
			fprintf(stderr,
				"%s: expected sections %d before footer, got sections %d\n",
				filename, expect, seen
			);
			fclose(f);
			exit(-1);
		  default:
			fprintf(stderr, "%s: unexpected section %#08x\n",
					filename, section_header.section);
			fclose(f);
			exit(-1);
		}
	}
}
 
void teardown(void) {
	teardown_steps();
	teardown_pieces();
	teardown_sym_set();
	teardown_symmetries();
	teardown_vec();
	teardown_clock();
}

void setup(void) {
	char* piece_file = "results/p" STRINGIFY(NBASE);

	setup_clock();
	setup_vec();
	setup_symmetries();
	setup_sym_set();
	setup_pieces();
	import_pieces(piece_file);
	setup_steps();
}

int main(int argc, char** argv) {
	uint first;
	double t1, t2;
	counter prev;

	setup();

	t1 = TIMETHIS({
		try_recurse(0);
	});
	printf("%u: Total %llu/%llu (%.2f)\n", NBASE, sym_result, all_result, t1);
	teardown();
	return 0;
}

