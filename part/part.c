#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/times.h>
#include "avl.h"

typedef unsigned int uint;
typedef unsigned long long counter;
typedef unsigned char vec;
typedef unsigned char mapping;

uint n = 0;
uint nodes;
uint vecsize;
counter sym_result;
counter all_result;
int bit_lookup[256];
uint bit_count[256];
vec transform[256];

vec* connections;
uint sym_count;
mapping* symmetries;
vec* canonical_v;

uint pieces_size;
uint piece_array_size;
uint pieces_used;
uint piece_array_used;
uint* piece_array;
vec* pieces;

vec* fullvec;
vec* piecestack;
vec* filledstack;
vec* solution;
avl_tree** solutions_seen;

int clk_tck;
void init_clk(void) {
	clk_tck = sysconf(_SC_CLK_TCK);
}

double difftime(clock_t t0, clock_t t1) {
	return ((double)t1 - t0) / clk_tck;
}
int curtime(void) {
	struct tms t;
	times(&t);
	return (int) t.tms_utime;
}

inline vec* connect_vec(uint i) {
	return connections + i * vecsize;
}

inline void vec_zero(vec* v) {
	memset(v, 0, vecsize);
}

inline void vec_copy(vec* src, vec* dest) {
	memcpy(dest, src, vecsize);
}

inline void vec_setbit(vec* v, uint i) {
	v[i >> 3] |= 1 << (i & 7);
}

inline void vec_clearbit(vec* v, uint i) {
	v[i >> 3] &= ~(1 << (i & 7));
}

inline uint vec_testbit(vec* v, uint i) {
	return (v[i >> 3] & (1 << (i & 7))) ? 1 : 0;
}

inline uint vec_bitcount(vec* v) {
	uint c = 0, i;
	for (i = 0; i < vecsize; ++i)
		c += bit_count[v[i]];
	return c;
}

#define DO_VEC(state) { \
	uint i; \
	for (i = 0; i < vecsize; ++i) { \
		state; \
	} \
}

inline void vec_or(vec* src, vec* dest)
	DO_VEC(dest[i] |= src[i])
inline void vec_and(vec* src, vec* dest)
	DO_VEC(dest[i] &= src[i])
inline void vec_xor(vec* src, vec* dest)
	DO_VEC(dest[i] ^= src[i])
inline void vec_or3(vec* s1, vec* s2, vec* dest)
	DO_VEC(dest[i] = s1[i] | s2[i])
inline void vec_and3(vec* s1, vec* s2, vec* dest)
	DO_VEC(dest[i] = s1[i] & s2[i])
inline void vec_xor3(vec* s1, vec* s2, vec* dest)
	DO_VEC(dest[i] = s1[i] ^ s2[i])

inline int vec_empty(vec* v) {
	uint i;
	for (i = 0; i < vecsize; ++i) {
		if (v[i])
			return 0;
	}
	return 1;
}

inline int vec_cmp(vec* s1, vec* s2) {
	uint i;
	signed int c;
	for (i = 0; i < vecsize; ++i) {
		c = transform[s1[i]] - transform[s2[i]];
		if (c > 0)
			return 1;
		if (c < 0)
			return -1;
	}
	return 0;
}

int vec_comparator(vec* s1, vec* s2, uint size) {
	return vec_cmp(s1, s2);
}

int set_comparator(vec* s1, vec* s2, uint size) {
	uint i;
	signed int c;
	for (i = 0; i < size; ++i) {
		c = transform[s1[i]] - transform[s2[i]];
		if (c > 0)
			return 1;
		if (c < 0)
			return -1;
	}
	return 0;
}

inline mapping* sym_map(uint i) {
	return symmetries + i * nodes;
}

inline void apply_map2(mapping* m, vec* src, vec* dest) {
	uint i, map;
	for (i = 0; i < nodes; ++i) {
		if (vec_testbit(src, m[i]))
			vec_setbit(dest, i);
		else
			vec_clearbit(dest, i);
	}
}

int first_bit(vec* v) {
	uint i;
	for (i = 0; i < vecsize; ++i) {
		if (v[i])
			return bit_lookup[v[i]] + (i << 3);
	}
	return -1;
}

int next_bit(vec* v, int first) {
	uint i = (first + 1) >> 3;
	vec c;
	if (i > vecsize)
		return -1;
	c = v[i] & ~((1 << ((first + 1) & 7)) - 1);
	if (c)
		return bit_lookup[c] + (i << 3);
	for (++i; i < vecsize; ++i) {
		if (v[i])
			return bit_lookup[v[i]] + (i << 3);
	}
	return -1;
}

inline vec* piecestack_vec(uint index) {
	return piecestack + index * vecsize;
}
inline vec* filledstack_vec(uint index) {
	return filledstack + index * vecsize;
}
inline vec* solution_vec(uint index) {
	return solution + index * vecsize;
}
inline vec* pieces_vec(uint index) {
	return pieces + index * vecsize;
}
inline vec* pieces_for(uint size) {
	return pieces_vec(piece_array[size]);
}
inline uint piece_count(uint size) {
	return piece_array[size + 1] - piece_array[size];
}

void dump_solution(FILE* stream, vec* v, uint size) {
	uint i, j;

	for (i = 0; i < size; ++i) {
		if (i > 0)
			fprintf(stream, " ");
		for (j = 0; j < nodes; ++j)
			fprintf(stream, vec_testbit(v, j) ? "1" : "0");
		v += vecsize;
	}
}

/*
  Given a vector I<v> of I<nodes> bits, returns the index of the symmetry
  that transforms it to canonical form. The canonicalized vector will be
  in canonical_v immediately after the call.
*/
uint canonical_piece(vec* v) {
	uint i, best_i;
	vec w[vecsize];
	mapping* m;

	w[0] = 0;
	best_i = 0;
	vec_copy(v, canonical_v);
	for (i = 1; i < sym_count; ++i) {
		m = sym_map(i);
		apply_map2(m, v, w);
		if (vec_cmp(w, canonical_v) > 0) {
			best_i = i;
			vec_copy(w, canonical_v);
		}
	}
	return best_i;
}

/*
  Given a set I<s> of vectors, returns the index of the symmetry that
  transforms it to canonical form. If the first shape in the set is an
  a_0-omino, that means it is the symmetry map under which all of the
  a_0-ominoes in the set sort lexically first (under vec_cmp); remaining
  shapes are used for tie-breaking, next-largest pieces first.
  The canonical set of vectors is in solution[] immediately after calling.
*/
uint canonical_set(vec* v, uint pieces) {
	uint i, this_size, cur_size, sym_i, best_sym;
	uint group_start[pieces], group_size[pieces], groups;
	int cmp;
	vec *source, *dest, tmpstack[pieces * vecsize];
	mapping* map;

	/* find the size of each piece, to group them by size */
	groups = 0;
	for (i = 0; i < pieces; ++i) {
		vec* source = piecestack_vec(i);
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

	memcpy(solution, v, pieces * vecsize);
	/* and make sure it's sorted canonically */
	for (i = 0; i < groups; ++i) {
		if (group_size[i] > 1)
			qsort(
				(void*)(solution + group_start[i] * vecsize),
				group_size[i], vecsize, (__compar_fn_t)vec_cmp
			);
	}

	best_sym = 0;
	for (sym_i = 1; sym_i < sym_count; ++sym_i) {
		/* apply the map */
		map = sym_map(sym_i);
		for (i = 0; i < pieces; ++i) {
			source = piecestack_vec(i);
			dest = tmpstack + i * vecsize;
			apply_map2(map, source, dest);
		}

		/* sort each group of same-sized pieces */
		for (i = 0; i < groups; ++i) {
			if (group_size[i] > 1)
				qsort(
					(void*)(tmpstack + group_start[i] * vecsize),
					group_size[i], vecsize, (__compar_fn_t)vec_cmp
				);
		}

		/* compare against best so far */
		for (i = 0; i < pieces; ++i) {
			cmp = vec_cmp(tmpstack + i * vecsize, solution + i * vecsize);
			if (cmp < 0)
				break;
			else if (cmp > 0)
				break;
		}
		if (cmp < 0) {
			memcpy(solution, tmpstack, pieces * vecsize);
			best_sym = sym_i;
		}
	}
	return best_sym;
}

void init_pieces(void) {
	canonical_v = (vec*)malloc(vecsize);

	piece_array_size = nodes + 2;
	piece_array_used = 1;
	piece_array = (uint*)malloc(sizeof(uint) * piece_array_size);
	piece_array[1] = 0;
	piece_array[2] = 1;

	pieces_size = 100;
	pieces_used = 1;
	pieces = (vec*)calloc(pieces_size, vecsize);

	vec_setbit(pieces_vec(0), 0);
	canonical_piece(pieces_vec(0));
	vec_copy(canonical_v, pieces_vec(0));
}

void prep_pieces(uint size) {
	uint smalli, small_lim;
	vec scratch[vecsize], *smallv;
	uint new;
	avl_tree* seen;
clock_t t0, t1;

	if (size <= piece_array_used)
		return;
	if (size > piece_array_used + 1)
		prep_pieces(size - 1);
t0 = curtime();
	small_lim = piece_array[size];
	seen = avl_new(vecsize);
	for (smalli = piece_array[size - 1]; smalli < small_lim; ++smalli) {
		smallv = pieces_vec(smalli);
		for (new = 0; new < nodes; ++new) {
			if (vec_testbit(smallv, new))
				continue;
			vec_and3(smallv, connect_vec(new), scratch);
			if (vec_empty(scratch))
				continue;
			vec_copy(smallv, scratch);
			vec_setbit(scratch, new);
			canonical_piece(scratch);
			if (avl_seen(seen, canonical_v) == AVL_EXISTS)
				continue;
			if (pieces_used >= pieces_size) {
				pieces_size *= 1.5;
				pieces = (vec*)realloc(pieces, pieces_size * vecsize);
				smallv = pieces_vec(smalli);
			}
			vec_copy(canonical_v, pieces_vec(pieces_used));
			++pieces_used;
		}
	}
	avl_delete(seen);
	piece_array[size + 1] = pieces_used;
	piece_array_used = size + 1;
t1 = curtime();
	fprintf(stderr, "P%u: %u (%.2f)\n", size, pieces_used - small_lim, difftime(t0, t1));
}

/*
  Given a vector I<v> of I<nodes> bits, returns C<TRUE> if it represents
  a single connected piece, else C<FALSE>.
  If no bits in the vector are set, behaviour is undetermined.
*/
int connected(vec* v) {
	int i = first_bit(v);
	vec unallocated[vecsize];
	vec current[vecsize];
	vec next[vecsize];

	vec_copy(v, unallocated);
	vec_clearbit(unallocated, i);
	vec_zero(current);
	vec_setbit(current, i);

	while (!vec_empty(current)) {
		vec_zero(next);
		for (i = first_bit(current); i >= 0; i = next_bit(current, i)) {
			vec* w = connect_vec(i);
			vec_or(w, next);
		}
		vec_and3(unallocated, next, current);
		vec_xor(current, unallocated);
	}
	return vec_empty(unallocated) ? 1 : 0;
}

void init_connections(void) {
	uint i, j;
	connections = (vec*)calloc(nodes, vecsize);
	for (i = 0; i < nodes; ++i) {
		vec* v = connect_vec(i);
		for (j = 0; j < n; ++j) {
			uint k = i ^ (1 << j);
			vec_setbit(v, k);
		}
	}
}

void init_bit_lookup(void) {
	uint i;
	bit_lookup[0] = -1;
	bit_count[0] = 0;
	for (i = 1; i < 256; ++i) {
		if (i & 1) {
			bit_lookup[i] = 0;
			bit_count[i] = 1 + bit_count[i >> 1];
		} else {
			bit_lookup[i] = 1 + bit_lookup[i >> 1];
			bit_count[i] = bit_count[i >> 1];
		}
	}
}

void next_perm(mapping* p, uint size) {
	uint last = size - 1;
	uint j, k;
	mapping temp;
	int i = (int)last - 1;

	while (i >= 0 && p[i] > p[i+1])
		--i;
	if (i < 0)
		return;
	for (j = i + 1, k = last; j < k; ++j, --k) {
		temp = p[j];
		p[j] = p[k];
		p[k] = temp;
	}
	for (j = i + 1; p[j] < p[i]; ++j)
		;
	temp = p[i];
	p[i] = p[j];
	p[j] = temp;
}

void init_transform(void) {
	uint i, j;
	vec value;
	for (i = 0; i < 256; ++i) {
		value = 0;
		for (j = 0; j < 8; ++j) {
			if (i & (1 << j))
				value |= 1 << (7 - j);
		}
		transform[i] = value;
	}
}

void init_symmetries(void) {
	uint fac = 1;
	uint i, j, k;
	mapping *m, *m2, value;
	mapping perm[n];

	for (i = 2; i <= n; ++i)
		fac *= i;
	sym_count = nodes * fac;
	symmetries = (mapping*)calloc(sym_count, nodes * sizeof(mapping));

	for (i = 0; i < n; ++i) {
		perm[i] = (mapping)i;
	}
	for (i = 0; i < fac; ++i, next_perm(perm, n)) {
		m = sym_map(i);
		for (j = 0; j < nodes; ++j) {
			value = 0;
			for (k = 0; k < n; ++k) {
				if (j & (1 << perm[k]))
					value |= 1 << k;
			}
			m[j] = value;
		}
	}

	m = sym_map(0);
	for (i = 1; i < nodes; ++i) {
		m2 = sym_map(i * fac);
		for (j = 0; j < nodes * fac; ++j)
			m2[j] = m[j] ^ i;
	}
}

void init_stack(void) {
	uint i;
	piecestack = (vec*)calloc(nodes, vecsize);
	filledstack = (vec*)calloc(nodes, vecsize);
	solution = (vec*)calloc(nodes, vecsize);

	solutions_seen = (avl_tree**)malloc(nodes * sizeof(avl_tree*));
	for (i = 0; i < nodes; ++i) {
		solutions_seen[i] = avl_new(vecsize * (i + 1));
	}

	fullvec = (vec*)calloc(1, vecsize);
	memset(fullvec, -1, vecsize);
	if (nodes < 8) {
		fullvec[0] = (1 << nodes) - 1;
	}
}

void check_solution(uint pieces) {
	canonical_set(piecestack, pieces);
	if (avl_seen(solutions_seen[pieces - 1], solution) == AVL_EXISTS)
		return;

	dump_solution(stdout, solution, pieces);
	printf("\n");
	++sym_result;
	++all_result;
}

void try_recurse(
	uint prev_level, uint remain, uint maxpiece, uint prev_index, avl_tree* prev_seen
) {
	uint level = prev_level + 1;
	uint size, piece_index, end_index, sym_index;
	vec *prev_filled, *this_filled, *this_piece, *this_shape;
	vec scratch[vecsize];
	avl_tree* seen;

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
			seen = avl_dup(prev_seen);
			piece_index = prev_index;
		} else {
			prep_pieces(size);
			seen = avl_new(vecsize);
			piece_index = piece_array[size];
		}
		end_index = piece_array[size + 1];
		for ( ; piece_index < end_index; ++piece_index) {
			this_piece = pieces_vec(piece_index);
			for (sym_index = 0; sym_index < sym_count; ++sym_index) {
				apply_map2(sym_map(sym_index), this_piece, this_shape);
				vec_and3(this_shape, prev_filled, scratch);
				if (vec_cmp(this_shape, scratch) != 0)
					continue;
				if (avl_seen(seen, this_shape) == AVL_EXISTS)
					continue;
				vec_xor3(prev_filled, this_shape, this_filled);
				try_recurse(level, remain - size, size, piece_index, seen);
			}
		}
		avl_delete(seen);
	}
}

void try_first(uint first) {
	uint piece_index, end_index;
	vec *this_piece, *ps, *fs;
	avl_tree* seen;

	prep_pieces(first);
	piece_index = piece_array[first];
	end_index = piece_array[first + 1];
	ps = piecestack_vec(0);
	fs = filledstack_vec(0);

	seen = avl_new(vecsize);
	for ( ; piece_index < end_index; ++piece_index) {
		this_piece = pieces_vec(piece_index);
		vec_copy(this_piece, ps);
		vec_xor3(this_piece, fullvec, fs);
		avl_seen(seen, this_piece); /* cannot exist */
		try_recurse(
			0,				/* recursion depth */
			nodes - first,	/* remaining bits in filledstack */
			first,			/* max size for pieces */
			piece_index,	/* min index for pieces (if same size) */
			seen			/* heap of seen shapes */
		);
	}
	avl_delete(seen);
}

void teardown(void) {
	uint i;
	free(piece_array);
	free(pieces);
	free(canonical_v);
	free(connections);
	free(symmetries);
	free(piecestack);
	free(filledstack);
	free(solution);
	free(fullvec);
	for (i = 0; i < nodes; ++i) {
		avl_delete(solutions_seen[i]);
	}
	free(solutions_seen);
}

int main(int argc, char** argv) {
	uint first;
	clock_t t0 = curtime(), t1;

	if (argc > 1)
		n = atoi(argv[1]);
	if (argc != 2 || n < 1 || n > 31) {
		fprintf(stderr, "Usage: %s <dimension>\n", argv[0]);
		return -1;
	}
	nodes = 1 << n;
	vecsize = nodes >> 3;
	if (vecsize < 1) vecsize = 1;

	init_clk();
	init_connections();
	init_bit_lookup();
	init_transform();
	init_symmetries();
	init_pieces();
	init_stack();

#if 1
	for (first = nodes; first > 0; --first) {
		try_first(first);
	}
#else
	prep_pieces(nodes);
#endif
	t1 = curtime();
	printf("%u: %llu, %llu (%.2f)\n", n, sym_result, all_result, difftime(t0, t1));
	teardown();
	return 0;
}

