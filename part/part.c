#include "part.h"
#include "vec.h"
#include "set.h"
#include "symmetries.h"
#include "pieces.h"
#include "clock.h"

#define REPORT_MASK ((1 << 16) - 1)

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
	vec_t *prev_free, *this_free, *this_shape;
	piece_t* this_piece;
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
			seen = vech_new();
			piece_index = piece_array[size];
		}
		end_index = piece_array[size + 1];
		for ( ; piece_index < end_index; ++piece_index) {
			this_piece = pieces_vec(piece_index);
			for (sym_index = 0; sym_index < this_piece->ss->count; ++sym_index) {
				apply_map2(sym_map(this_piece->ss->index[sym_index]), &(this_piece->v), this_shape);
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
	piece_t *this_piece;
	vech_tree* seen;

	piece_index = piece_array[first];
	end_index = piece_array[first + 1];

	seen = vech_new();
	for ( ; piece_index < end_index; ++piece_index) {
		this_piece = pieces_vec(piece_index);
		vec_copy(&(this_piece->v), &(step->shape));
		vec_not2(&(this_piece->v), &(step->freevec));
		set_init(&(step->set), &(step->shape));
		vech_seen(seen, &(this_piece->v)); /* cannot exist */
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
			fprintf(stderr, "pieces %u: %u [%u symsets] (%.2f)\n",
					first, piece_array[first + 1] - piece_array[first],
					sym_set_count, t2);
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

