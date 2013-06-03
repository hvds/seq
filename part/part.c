#include "part.h"
#include "vec.h"
#include "set.h"
#include "symmetries.h"
#include "clock.h"

#define REPORT_MASK ((1 << 20) - 1)

#define BLOCK_START
#define BLOCK_END

BLOCK_START
	/* pieces functions */
	uint pieces_size;
	uint piece_array_size;
	uint pieces_used;
	uint piece_array_used;
	uint* piece_array;
	vec_t* pieces;
	vec_t* canonical_v;

	inline vec_t* pieces_vec(uint index) {
		return pieces + index;
	}
	inline vec_t* pieces_for(uint size) {
		return pieces_vec(piece_array[size]);
	}

	/*
	  Given a vector I<v> of I<NODES> bits, returns the index of the symmetry
	  that transforms it to canonical form. The canonicalized vector will be
	  in canonical_v immediately after the call.
	*/
	uint canonical_piece(vec_t* v) {
		uint i, j, best_i;
		uint first, sc_start, sc_end;
		uint bits = 0;
		uint contiguous = 0;
		vec_t w;
		sym_t* sym;
		/* The (sorted) symmetries evenly map each bit to the zero bit */
		uint sym_block = sym_count >> NBASE;

		/* count the bits, and the number contiguous at start */
		for (i = 0; i < NODES; ++i) {
			if (vec_testbit(v, i)) {
				++bits;
				if (i == contiguous)
					++contiguous;
			}
		}

		w.v[0] = 0;
		best_i = 0;
		vec_copy(v, canonical_v);

		if (bits == 0)
			return 0;

		for (first = 0; first < NODES; ++first) {
			/* if this bit is not set, we can skip the whole block */
			if (!vec_testbit(v, first))
				continue;
			sc_start = sym_block * first;
			sc_end = sc_start + sym_block;
			for (i = sc_start; i < sc_end; ++i) {
				sym = sym_map(i);
				for (j = 1; j < contiguous; ++j)
					if (!vec_testbit(v, sym->map[j]))
						goto SYM_FAIL;
				apply_map2(sym, v, &w);
				if (vec_cmp(&w, canonical_v) > 0) {
					best_i = i;
					vec_copy(&w, canonical_v);
					for (j = contiguous; j < bits; ++j) {
						if (vec_testbit(&w, j)) {
							++contiguous;
						} else {
							break;
						}
					}
				}
			  SYM_FAIL:
				;
			}
		}
		return best_i;
	}

	void setup_pieces(void) {
		canonical_v = (vec_t*)malloc(sizeof(vec_t));

		piece_array_size = NODES + 2;
		piece_array_used = 1;
		piece_array = (uint*)malloc(sizeof(uint) * piece_array_size);
		piece_array[1] = 0;
		piece_array[2] = 1;

		pieces_size = 100;
		pieces_used = 1;
		pieces = (vec_t*)calloc(pieces_size, sizeof(vec_t));

		vec_setbit(pieces_vec(0), 0);
		canonical_piece(pieces_vec(0));
		vec_copy(canonical_v, pieces_vec(0));
	}

	void teardown_pieces(void) {
		free(pieces);
		free(piece_array);
		free(canonical_v);
	}

	void prep_pieces(uint size) {
		uint smalli, small_lim;
		vec_t scratch, *smallv;
		uint new;
		vech_tree* seen;
		double t;

		if (size <= piece_array_used)
			return;
		if (size > piece_array_used + 1)
			prep_pieces(size - 1);
		small_lim = pieces_used;
		seen = vech_new();
		for (smalli = piece_array[size - 1]; smalli < small_lim; ++smalli) {
			smallv = pieces_vec(smalli);
			for (new = 0; new < NODES; ++new) {
				if (vec_testbit(smallv, new))
					continue;
				vec_and3(smallv, connect_vec(new), &scratch);
				if (vec_empty(&scratch))
					continue;
				vec_copy(smallv, &scratch);
				vec_setbit(&scratch, new);
				canonical_piece(&scratch);
				if (vech_seen(seen, canonical_v) == VECH_EXISTS)
					continue;
				if (pieces_used >= pieces_size) {
					pieces_size *= 1.5;
					pieces = (vec_t*)realloc(pieces, pieces_size * sizeof(vec_t));
					smallv = pieces_vec(smalli);
				}
				vec_copy(canonical_v, pieces_vec(pieces_used));
				++pieces_used;
			}
		}
		vech_delete(seen);
		piece_array[size + 1] = pieces_used;
		piece_array_used = size;
	}
BLOCK_END

counter sym_result;

typedef struct step_s {
	set_t set;		/* the combined set of pieces so far */
	vec_t freevec;	/* bits still free */
	vec_t shape;	/* a rotated/reflected piece to insert */
} step_t;
step_t steps[NODES];
set_t* solution;
seth_tree** solutions_seen;

void setup_stack(void) {
	uint i;
	solution = (set_t*)malloc(sizeof(set_t));
	solutions_seen = (seth_tree**)malloc(NODES * sizeof(seth_tree*));
	for (i = 0; i < NODES; ++i)
		solutions_seen[i] = seth_new();
}

void teardown_stack(void) {
	uint i;
	for (i = 0; i < NODES; ++i)
		seth_delete(solutions_seen[i]);
	free(solutions_seen);
	free(solution);
}

/*
  Given a set_t I<s>, returns the index of the symmetry that
  transforms it to canonical form. This is the symmetry mapping under
  which it sorts lexically first.
  The canonical set is in I<solution> immediately after calling.
*/
uint canonical_set(set_t* set, uint pieces) {
	uint i, sym_i, best_sym;
	uint seen, mapping[NODES + 1], source, dest;
	sym_t* sym;

	/* the set may not start in canonical form, so a simple set_copy() here
	 * is not enough.
	 */
	sym_i = 0;
	sym = sym_map(sym_i);
	seen = 0;
	memset(mapping, 0, sizeof(mapping));
	for (i = 0; i < NODES; ++i) {
		source = set->p[sym->map[i]];
		if (source && !mapping[source])
			mapping[source] = ++seen;
		solution->p[i] = mapping[source];
	}
	best_sym = 0;

	for (sym_i = 1; sym_i < sym_count; ++sym_i) {
		sym = sym_map(sym_i);

		/* while applying the map, at each step we either know it is worse
		 * (and abort) or that it is better (in which case we go to the
		 * accept state) or that it is identical so far. So we need no
		 * scratch space for the result, only for the assignments.
		 */
		seen = 0;
		memset(mapping, 0, sizeof(mapping));

		for (i = 0; i < NODES; ++i) {
			source = set->p[sym->map[i]];
			if (source && !mapping[source])
				mapping[source] = ++seen;
			dest = mapping[source];
			if (dest > solution->p[i])
				goto SYM_FAIL;
			if (dest < solution->p[i])
				goto SYM_ACCEPT;
		}
		/* if we fall through it's identical, treat same as fail */
	  SYM_FAIL:
		continue;

	  SYM_ACCEPT:
		for ( ; i < NODES; ++i) {
			source = set->p[sym->map[i]];
			if (source && !mapping[source]) 
				mapping[source] = ++seen;
			solution->p[i] = mapping[source];
		}
		best_sym = sym_i;
	}
	return best_sym;
}

void check_solution(set_t* set, uint pieces) {
	canonical_set(set, pieces);
	if (seth_seen(solutions_seen[pieces - 1], solution) == SETH_EXISTS)
		return;

	fprint_set(stdout, solution, pieces);
	printf("\n");
	++sym_result;
	if (! (sym_result & REPORT_MASK)) {
		fprintf(stderr, "%llu: ", sym_result);
		fprint_set(stderr, solution, pieces);
		fprintf(stderr, "\n");
	}
}

void try_recurse(
	uint prev_level, uint remain, uint maxpiece, uint prev_index,
	vech_tree* prev_seen
) {
	uint level = prev_level + 1;
	uint size, piece_index, end_index, sym_index;
	step_t* step = &(steps[level]);
	step_t* prev_step = &(steps[prev_level]);
	vec_t *prev_free, *this_free, *this_piece, *this_shape;
	vech_tree* seen;

	if (remain == 0) {
		/* we have a solution in the stack: check for uniqueness, and record */
		check_solution(&(prev_step->set), level);
		return;
	}
	prev_free = &(prev_step->freevec);
	this_free = &(step->freevec);
	this_shape = &(step->shape);
	size = (remain < maxpiece) ? remain : maxpiece;
	for ( ; size > 1; --size) {
		if (size == maxpiece) {
			seen = vech_dup(prev_seen);
			piece_index = prev_index;
		} else {
			prep_pieces(size);
			seen = vech_new();
			piece_index = piece_array[size];
		}
		end_index = piece_array[size + 1];
		for ( ; piece_index < end_index; ++piece_index) {
			this_piece = pieces_vec(piece_index);
			for (sym_index = 0; sym_index < sym_count; ++sym_index) {
				apply_map2(sym_map(sym_index), this_piece, this_shape);
				if (! vec_contains(prev_free, this_shape))
					continue;
				if (vech_seen(seen, this_shape) == VECH_EXISTS)
					continue;
				vec_xor3(prev_free, this_shape, this_free);
				set_append(&(step->set), &(prev_step->set), this_shape, level);
				try_recurse(
					level,			/* recursion depth */
					remain - size,	/* remaining bits free */
					size,			/* max size for pieces */
					piece_index,	/* min index for pieces (if same size) */
					seen			/* heap of seen shapes */
				);
			}
		}
		vech_delete(seen);
	}
	/* handle size=1 separately */
	{
		uint cur = level, i;
		set_t* s = &(step->set);
		set_copy(&(prev_step->set), s);
		for (i = 0; i < NODES; ++i) {
			if (s->p[i] != 0)
				continue;
			s->p[i] = ++cur;
		}
		check_solution(s, cur);
	}
}

void try_first(uint first) {
	uint piece_index, end_index;
	step_t* step = &(steps[0]);
	vec_t *this_piece;
	vech_tree* seen;

	piece_index = piece_array[first];
	end_index = piece_array[first + 1];

	seen = vech_new();
	for ( ; piece_index < end_index; ++piece_index) {
		this_piece = pieces_vec(piece_index);
		vec_copy(this_piece, &(step->shape));
		vec_not2(this_piece, &(step->freevec));
		set_init(&(step->set), &(step->shape));
		vech_seen(seen, this_piece); /* cannot exist */
		try_recurse(
			0,				/* recursion depth */
			NODES - first,	/* remaining bits free */
			first,			/* max size for pieces */
			piece_index,	/* min index for pieces (if same size) */
			seen			/* heap of seen shapes */
		);
	}
	vech_delete(seen);
}

void teardown(void) {
	teardown_stack();
	teardown_pieces();
	teardown_symmetries();
	teardown_vec();
	teardown_clock();
}

void setup(void) {
	setup_clock();
	setup_vec();
	setup_symmetries();
	setup_pieces();
	setup_stack();
}

int main(int argc, char** argv) {
	uint first;
	double t1, t2;
	counter prev;

	setup();

	t1 = TIMETHIS({
		for (first = 1; first <= NODES; ++first) {
			t2 = TIMETHIS({
				prep_pieces(first);
			});
			fprintf(stderr, "pieces %u: %u (%.2f)\n",
					first, piece_array[first + 1] - piece_array[first], t2);
		}
	});
	fprintf(stderr, "Total pieces: %u (%.2f)\n", pieces_used, t1);
#ifndef PIECE_ONLY
	t1 = TIMETHIS({
#ifndef REVERSE
		for (first = 1; first <= NODES; ++first) {
#else
		for (first = NODES; first > 0; --first) {
#endif
			prev = sym_result;
			t2 = TIMETHIS({
				try_first(first);
			});
			fprintf(stderr, "solutions %u: %llu/%llu (%.2f)\n",
					first, sym_result - prev, sym_result, t2);
		}
	});
	printf("%u: Total %llu (%.2f)\n", NBASE, sym_result, t1);
#endif
	teardown();
	return 0;
}

