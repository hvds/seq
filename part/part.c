#include "part.h"
#include "vec.h"
#include "set.h"
#include "symmetries.h"
#include "clock.h"

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
	inline uint piece_count(uint size) {
		return piece_array[size + 1] - piece_array[size];
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
		t = TIMETHIS({
			small_lim = piece_array[size];
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
			piece_array_used = size + 1;
		});
		fprintf(stderr, "P%u: %u (%.2f)\n", size, pieces_used - small_lim, t);
	}
BLOCK_END

counter sym_result;
counter all_result;

vec_t* piecestack;
vec_t* filledstack;
vec_t* solution;
seth_tree** solutions_seen;

inline vec_t* piecestack_vec(uint index) {
	return piecestack + index;
}
inline vec_t* filledstack_vec(uint index) {
	return filledstack + index;
}
inline vec_t* solution_vec(uint index) {
	return solution + index;
}

void setup_stack(void) {
	uint i;
	piecestack = (vec_t*)calloc(NODES, sizeof(vec_t));
	filledstack = (vec_t*)calloc(NODES, sizeof(vec_t));
	solution = (vec_t*)calloc(NODES, sizeof(vec_t));
	solutions_seen = (seth_tree**)malloc(NODES * sizeof(seth_tree*));
	for (i = 0; i < NODES; ++i)
		solutions_seen[i] = seth_new(i + 1);
}

void teardown_stack(void) {
	uint i;
	for (i = 0; i < NODES; ++i)
		seth_delete(solutions_seen[i]);
	free(solutions_seen);
	free(solution);
	free(filledstack);
	free(piecestack);
}

void dump_solution(FILE* stream, vec_t* v, uint size) {
	uint i, j;

	for (i = 0; i < size; ++i) {
		if (i > 0)
			fprintf(stream, " ");
		for (j = 0; j < NODES; ++j)
			fprintf(stream, vec_testbit(v, j) ? "1" : "0");
		++v;
	}
}

/*
  Given a set I<s> of vectors, returns the index of the symmetry that
  transforms it to canonical form. If the first shape in the set is an
  a_0-omino, that means it is the symmetry map under which all of the
  a_0-ominoes in the set sort lexically first (under vec_cmp); remaining
  shapes are used for tie-breaking, next-largest pieces first.
  The canonical set of vectors is in solution[] immediately after calling.
*/
uint canonical_set(vec_t* v, uint pieces) {
	uint i, this_size, cur_size, sym_i, best_sym;
	uint group_start[pieces], group_size[pieces], groups;
	int cmp;
	vec_t *source, *dest, tmpstack[pieces];
	sym_t* sym;

	/* find the size of each piece, to group them by size */
	groups = 0;
	for (i = 0; i < pieces; ++i) {
		vec_t* source = piecestack_vec(i);
		this_size = vec_bitcount(source);
		if (i == 0 || this_size != cur_size) {
			groups = (i == 0) ? 0 : groups + 1;
			group_start[groups] = i;
			group_size[groups] = 1;
			cur_size = this_size;
		} else {
			++group_size[groups];
		}
	}
	++groups;

	memcpy(solution, v, pieces * sizeof(vec_t));
	/* and make sure it's sorted canonically */
	for (i = 0; i < groups; ++i) {
		if (group_size[i] > 1)
			qsort(
				(void*)(solution + group_start[i]),
				group_size[i], sizeof(vec_t), (__compar_fn_t)vec_cmp
			);
	}

	best_sym = 0;
	for (sym_i = 1; sym_i < sym_count; ++sym_i) {
		/* apply the map */
		sym = sym_map(sym_i);
		for (i = 0; i < pieces; ++i) {
			source = piecestack_vec(i);
			dest = &tmpstack[i];
			apply_map2(sym, source, dest);
		}

		/* sort each group of same-sized pieces */
		for (i = 0; i < groups; ++i) {
			if (group_size[i] > 1)
				qsort(
					(void*)(tmpstack + group_start[i]),
					group_size[i], sizeof(vec_t), (__compar_fn_t)vec_cmp
				);
		}

		/* compare against best so far */
		for (i = 0; i < pieces; ++i) {
			cmp = vec_cmp(tmpstack + i, solution + i);
			if (cmp < 0)
				break;
			else if (cmp > 0)
				break;
		}
		if (cmp < 0) {
			memcpy(solution, tmpstack, pieces * sizeof(vec_t));
			best_sym = sym_i;
		}
	}
	return best_sym;
}

void check_solution(uint pieces) {
	canonical_set(piecestack, pieces);
	if (seth_seen(solutions_seen[pieces - 1], (set_t*)solution) == SETH_EXISTS)
		return;

	dump_solution(stdout, solution, pieces);
	printf("\n");
	++sym_result;
	++all_result;
}

void try_recurse(
	uint prev_level, uint remain, uint maxpiece, uint prev_index, vech_tree* prev_seen
) {
	uint level = prev_level + 1;
	uint size, piece_index, end_index, sym_index;
	vec_t *prev_filled, *this_filled, *this_piece, *this_shape;
	vech_tree* seen;

	if (remain == 0) {
		/* we have a solution in the stack: check for uniqueness, and record */
		check_solution(level);
		return;
	}
	prev_filled = filledstack_vec(prev_level);
	this_filled = filledstack_vec(level);
	this_shape = piecestack_vec(level);
	size = (remain < maxpiece) ? remain : maxpiece;
	for ( ; size > 0; --size) {
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
				if (! vec_contains(prev_filled, this_shape))
					continue;
				if (vech_seen(seen, this_shape) == VECH_EXISTS)
					continue;
				vec_xor3(prev_filled, this_shape, this_filled);
				try_recurse(level, remain - size, size, piece_index, seen);
			}
		}
		vech_delete(seen);
	}
}

void try_first(uint first) {
	uint piece_index, end_index;
	vec_t *this_piece, *ps, *fs;
	vech_tree* seen;

	prep_pieces(first);
	piece_index = piece_array[first];
	end_index = piece_array[first + 1];
	ps = piecestack_vec(0);
	fs = filledstack_vec(0);

	seen = vech_new();
	for ( ; piece_index < end_index; ++piece_index) {
		this_piece = pieces_vec(piece_index);
		vec_copy(this_piece, ps);
		vec_not2(this_piece, fs);
		vech_seen(seen, this_piece); /* cannot exist */
		try_recurse(
			0,				/* recursion depth */
			NODES - first,	/* remaining bits in filledstack */
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
	double t;

	setup();

	t = TIMETHIS({
#ifdef PIECE_ONLY
		prep_pieces(NODES);
#else
		for (first = NODES; first > 0; --first) {
			try_first(first);
		}
#endif
	});
	printf("%u: %llu, %llu (%.2f)\n", NBASE, sym_result, all_result, t);
	teardown();
	return 0;
}

