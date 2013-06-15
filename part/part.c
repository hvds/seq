#include "part.h"
#include "vec.h"
#include "set.h"
#include "symmetries.h"
#include "pieces.h"
#include "clock.h"
#include <assert.h>

#define REPORT_MASK ((1 << 18) - 1)

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
	uint pivot, broken_index;
	int broken_value;
	uint seen_index[NODES], piece_index[NODES];
	uint cmp;

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
		memset(seen_index, 0, piece_count * sizeof(uint));
		for (i = 0; i < NODES; ++i) {
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
				goto NOT_IDENTITY;
			}

			if (dest > pivot)
				continue;
			if (seen_index[dest] < piece_size)
				continue;
			/* we have exactly matched the pivot piece */
			while (++pivot < piece_count) {
				if (pivot == broken_index) {
					if (broken_value < 0) {
						return 0;
					}
					goto NOT_IDENTITY;
				}
				if (seen_index[pivot] < piece_size)
					break;
			}
			if (pivot >= piece_count)
				break;
		}
		/* if we fall through it's an identity, so preserve it */
		ssto->index[ssto->count++] = ssfrom->index[sym_i];
	  NOT_IDENTITY:
		;
	}
	/* if we get through the gauntlet, we are canonical */
	assert(ssto->count > 0);
	assert(0 == (sym_count % ssto->count));
	return 1;
}

void check_solution(uint level) {
	step_t* step = &(steps[level]);
	set_t combined;
	uint i;
	uint syms = sym_count / step->ss->count;
	assert(0 == (sym_count % step->ss->count));

	++sym_result;
	all_result += syms;

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
#ifdef NOPREP
		if (prev_level == 0) {
			prep_pieces(size);
		}
#endif
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
			sym_set_t* ss = this_piece->ss;
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
 
void teardown(void) {
	teardown_steps();
	teardown_pieces();
	teardown_sym_set();
	teardown_symmetries();
	teardown_vec();
	teardown_clock();
}

void setup(void) {
	setup_clock();
	setup_vec();
	setup_symmetries();
	setup_sym_set();
	setup_pieces();
	setup_steps();
}

int main(int argc, char** argv) {
	uint first;
	double t1, t2;
	counter prev;

	setup();

#ifndef NOPREP
	t1 = TIMETHIS({
		for (first = 1; first <= NODES; ++first) {
			t2 = TIMETHIS({
				prep_pieces(first);
			});
			fprintf(stderr, "pieces %u: %u [%u symsets] (%.2f)\n",
					first, piece_array[first + 1] - piece_array[first],
					sym_set_count, t2);
		}
	});
	fprintf(stderr, "Total pieces: %u (%.2f)\n", pieces_used, t1);
	reset_clock();
#endif

#ifndef PIECE_ONLY
	t1 = TIMETHIS({
		try_recurse(0);
	});
	printf("%u: Total %llu/%llu (%.2f)\n", NBASE, sym_result, all_result, t1);
#endif
	teardown();
	return 0;
}

