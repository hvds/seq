#define BUILD_PIECES

#include "part.h"
#include "clock.h"
#include "vec.h"
#define IS_SYM_SET_C
#include "sym_set.h"
#define IS_PIECES_C
#include "pieces.h"
#include "symmetries.h"

uint* sym_set_lookup;
char* sym_set_content;

uint pieces_size;
uint pieces_used;
uint piece_array_used;
uint piece_array[NODES + 2];
uint piece_count[NODES];
piece_t* pieces;
uint counts_size;
uint* counts;
vec_t* canonical_v;

char* sym_set_arena;
uint sym_set_size;
uint sym_set_used;
uint sym_set_count;

typedef unsigned int ssh;
typedef struct ssh_s {
	ssh left;
	ssh right;
	int balance;
	uint sym_set_offset;
} ssh_t;
typedef struct ssh_tree_s {
	ssh_t* ssharena;
	ssh ssha_size;
	ssh ssha_used;
	ssh ssha_root;
} ssh_tree_t;
ssh_tree_t* set_tree;

#define SSH_OK ((uint)1)
#define SSH_GROW ((uint)2)

inline ssh_t* SSH(ssh_tree_t* t, ssh index) {
	return (ssh_t*)(t->ssharena + (index - 1));
}
inline sym_set_t* SYMSET(uint offset) {
	return (sym_set_t*)(sym_set_arena + offset);
}
inline sym_set_t* SSHV(ssh_t* elem) {
	return SYMSET(elem->sym_set_offset);
}
inline signed int sym_set_cmp(sym_set_t* s1, sym_set_t* s2) {
	uint i, u1, u2;
	if (s1->count != s2->count)
		return s1->count > s2->count ? 1 : -1;
	for (i = 0; i < s1->count; ++i) {
		u1 = s1->index[i];
		u2 = s2->index[i];
		if (u1 != u2)
			return u1 > u2 ? 1 : -1;
	}
	return 0;
}

void setup_sym_set(void) {
	sym_set_count = 0;
	sym_set_used = 0;
	sym_set_size = 1024;
	sym_set_arena = (char*)malloc(sym_set_size);
	/* new tree */
	set_tree = (ssh_tree_t*)malloc(sizeof(ssh_tree_t));
	set_tree->ssha_size = 100;
	set_tree->ssha_used = 0;
	set_tree->ssha_root = 0;
	set_tree->ssharena = (ssh_t*)malloc(set_tree->ssha_size * sizeof(ssh_t));
}

void teardown_sym_set(void) {
	uint i;
	free(set_tree->ssharena);
	free(set_tree);
	free(sym_set_arena);
}

uint save_set(sym_set_t* set) {
	uint size = sizeof(sym_set_t) + set->count * sizeof(uint);
	uint offset = sym_set_used;
	++sym_set_count;
	sym_set_used += size;
	if (sym_set_used > sym_set_size) {
		sym_set_size = sym_set_size * 3 / 2;
		if (sym_set_used > sym_set_size) {
			sym_set_size = sym_set_used;
		}
		sym_set_arena = (char*)realloc(sym_set_arena, sym_set_size);
	}
	memcpy(SYMSET(offset), set, size);
	return offset;
}

void resize_ssh(ssh_tree_t* t, uint size) {
	if (size <= t->ssha_size)
		return;
	t->ssha_size = t->ssha_size * 3 / 2;
	if (size > t->ssha_size)
		t->ssha_size = size;
	t->ssharena = (ssh_t*)realloc(t->ssharena, t->ssha_size * sizeof(ssh_t));
}

uint ssh_alloc(ssh_tree_t* t, uint offset) {
	ssh ai = ++t->ssha_used;	/* 1-based; 0 means not present */
	ssh_t* ap = SSH(t, ai);
	ap->left = 0;
	ap->right = 0;
	ap->balance = 0;
	ap->sym_set_offset = offset;
	return ai;
}

void ssh_rotateleft(ssh_tree_t* t, ssh* root) {
	ssh a = *root;
	ssh b = SSH(t, a)->right;
	*root = b;
	SSH(t, a)->right = SSH(t, b)->left;
	SSH(t, b)->left = a;
}

void ssh_rotateright(ssh_tree_t* t, ssh* root) {
	ssh a = *root;
	ssh b = SSH(t, a)->left;
	*root = b;
	SSH(t, a)->left = SSH(t, b)->right;
	SSH(t, b)->right = a;
}

void ssh_rebalance(ssh_tree_t* t, ssh base) {
	ssh_t* bp = SSH(t, base);
	switch (bp->balance) {
	  case 0:
		SSH(t, bp->left)->balance = 0;
		SSH(t, bp->right)->balance = 0;
		break;
	  case -1:
		SSH(t, bp->left)->balance = 0;
		SSH(t, bp->right)->balance = 1;
		break;
	  case 1:
		SSH(t, bp->left)->balance = -1;
		SSH(t, bp->right)->balance = 0;
		break;
   }
   bp->balance = 0;
}

/* Find or insert an element into a SSH tree.
 * Returns the offset of the matching element if it was found, else
 * returns either SSH_OK (if the element was inserted without growing)
 * or SSH_GROW (if the element was inserted and the subtree grew).
 */
uint sshx_seen(ssh_tree_t* t, ssh* rootp, uint offset) {
	ssh root = *rootp;	/* note: invalidated by any rotate */
	ssh_t* rnp;
	uint result;

	if (!root) {
		*rootp = ssh_alloc(t, offset);
		return SSH_GROW;
	}
	rnp = SSH(t, root);

	switch (sym_set_cmp(SSHV(rnp), SYMSET(offset))) {
	  case 0:
		return rnp->sym_set_offset;
	  case 1:
		/* insert into the left subtree */
		if (!rnp->left) {
			rnp->left = ssh_alloc(t, offset);
			if (rnp->balance--)
				return SSH_OK;
			return SSH_GROW;
		}
		result = sshx_seen(t, &(rnp->left), offset);
		if (result != SSH_GROW)
			return result;
		switch (rnp->balance--) {
		  case 1:
			return SSH_OK;
		  case 0:
			return SSH_GROW;
		}
		if (SSH(t, rnp->left)->balance < 0) {
			ssh_rotateright(t, rootp);
			rnp = SSH(t, *rootp);
			rnp->balance = 0;
			SSH(t, rnp->right)->balance = 0;
		} else {
			ssh_rotateleft(t, &(rnp->left));
			ssh_rotateright(t, rootp);
			ssh_rebalance(t, *rootp);
		}
		return SSH_OK;
	  case -1:
		/* insert into the right subtree */
		if (!rnp->right) {
			rnp->right = ssh_alloc(t, offset);
			if (rnp->balance++)
				return SSH_OK;
			return SSH_GROW;
		}
		result = sshx_seen(t, &(rnp->right), offset);
		if (result != SSH_GROW) 
			return result;
		switch (rnp->balance++) {
		  case -1:
			return SSH_OK;
		  case 0:
			return SSH_GROW;
		}
		if (SSH(t, rnp->right)->balance > 0) {
			ssh_rotateleft(t, rootp);
			rnp = SSH(t, *rootp);
			rnp->balance = 0;
			SSH(t, rnp->left)->balance = 0;
		} else {
			ssh_rotateright(t, &(rnp->right));
			ssh_rotateleft(t, rootp);
			ssh_rebalance(t, *rootp);
		}
		return SSH_OK;
	}
}

uint ssh_seen(ssh_tree_t* tree, sym_set_t* set) {
	uint offset, result;

	/* do any realloc now, so the tree won't move under our feet */
	resize_ssh(tree, tree->ssha_used + 1);

	/* allocate a copy of the set */
	offset = save_set(set);
	result = sshx_seen(tree, &(tree->ssha_root), offset);
	if (result == SSH_GROW || result == SSH_OK) {
		/* keep the allocated copy */
		return offset;
	} else {
		/* discard the allocated copy, return offset of the copy */
		sym_set_used = offset;
		--sym_set_count;
		return result;
	}
}

sym_set_t* symset_new(void) {
	return (sym_set_t*)malloc(sizeof(sym_set_t) + sym_count * sizeof(int));
}
sym_set_t* symset_resize(sym_set_t* ss) {
	return (sym_set_t*)realloc(ss, sizeof(sym_set_t) + ss->count * sizeof(int));
}
void symset_delete(sym_set_t* ss) {
	free(ss);
}

uint distinct_symmetries_for_piece(vec_t* v) {
	vech_tree* seen = vech_new();
	sym_set_t* ss = symset_new();
	uint sso_real;
	vec_t scratch;
	uint i;

	ss->count = 0;
	for (i = 0; i < sym_count; ++i) {
		apply_map2(sym_map(i), v, &scratch);
		if (vech_seen(seen, &scratch) == VECH_EXISTS)
			continue;
		ss->index[ss->count++] = i;
	}
	vech_delete(seen);
	ss = symset_resize(ss);
	sso_real = ssh_seen(set_tree, ss);
	symset_delete(ss);
	return sso_real;
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

	piece_array_used = 1;
	piece_array[1] = 0;
	piece_array[2] = 1;

	pieces_size = 100;
	pieces_used = 1;
	pieces = (piece_t*)calloc(pieces_size, sizeof(piece_t));

	vec_setbit(&(pieces_vec(0)->v), 0);
	canonical_piece(&(pieces_vec(0)->v));
	vec_copy(canonical_v, &(pieces_vec(0)->v));
	pieces_vec(0)->sso = distinct_symmetries_for_piece(canonical_v);
}

void teardown_pieces(void) {
	free(pieces);
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
		smallv = &(pieces_vec(smalli)->v);
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
				pieces = (piece_t*)realloc(pieces, pieces_size * sizeof(piece_t));
				smallv = &(pieces_vec(smalli)->v);
			}
			vec_copy(canonical_v, &(pieces_vec(pieces_used)->v));
			pieces_vec(pieces_used)->sso = distinct_symmetries_for_piece(canonical_v);
			++pieces_used;
		}
	}
	vech_delete(seen);
	piece_array[size + 1] = pieces_used;
	piece_array_used = size;
}

void prep_counts(void) {
	uint i, j;
	uint* base;

	counts_size = 0;
	for (i = 1; i <= NODES; ++i) {
		piece_count[i - 1] = piece_array[i + 1] - piece_array[i];
		counts_size += i * piece_count[i - 1];
	}
	counts = (uint*)calloc(counts_size, sizeof(uint));
	/* init count(piece, 1) = 1 for all pieces */
	base = counts;
	for (i = 1; i <= NODES; ++i) {
		for (j = 0; j < piece_count[i - 1]; ++j) {
			base[j * i] = 1;
		}
		base += piece_count[i - 1] * i;
	}
}

/*
  Writes out the data for pieces, a sequence of sections of the form:
	typedef struct section_s {
		uint section;
		uint size;
		char data[size];
	};

  The sections are:
	section 0: sym_sets, data is of form:
		uint count;
		sym_set_t symsets[count];
	section 1: pieces, data is of form:
		uint count_of_size_n_plus_1[NODES];
		// symset offsets match the section 0 data
		piece_t pieces[total count];
	section 2: counts, data is of form:
		// entries in same order as for pieces above
		counts for size 1 pieces (1 * 4 bytes each)
		counts for size 2 pieces (2 * 4 bytes each)
		...
		counts for size NODES pieces (NODES * 4 bytes each)
	section 3: end marker (no data)
*/

void write_pieces(void) {
	uint section, size, i;

	section = 0;
	size = sizeof(uint) + sym_set_used * sizeof(char);
	fwrite(&section, sizeof(uint), 1, stdout);
	fwrite(&size, sizeof(uint), 1, stdout);
	fwrite(&sym_set_count, sizeof(uint), 1, stdout);
	fwrite(sym_set_arena, sizeof(char), sym_set_used, stdout);

	section = 1;
	size = NODES * sizeof(uint) + pieces_used * sizeof(piece_t);
	fwrite(&section, sizeof(uint), 1, stdout);
	fwrite(&size, sizeof(uint), 1, stdout);
	fwrite(piece_count, sizeof(uint), NODES, stdout);
	fwrite(pieces, sizeof(piece_t), pieces_used, stdout);

	section = 2;
	size = counts_size * sizeof(uint);
	fwrite(&section, sizeof(uint), 1, stdout);
	fwrite(&size, sizeof(uint), 1, stdout);
	fwrite(counts, sizeof(uint), counts_size, stdout);

	section = 3;
	size = 0;
	fwrite(&section, sizeof(uint), 1, stdout);
	fwrite(&size, sizeof(uint), 1, stdout);
}

void teardown(void) {
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
}

int main(int argc, char** argv) {
	unsigned int i;
	double t0;

	setup();
	for (i = 1; i <= NODES; ++i) {
		uint pc = pieces_used;
		uint sc = sym_set_count;
		t0 = TIMETHIS({
			prep_pieces(i);
		});
		fprintf(
			stderr, "pieces %u: %u/%u [%u/%u symsets] (%.2f)\n",
			i, pieces_used - pc, pieces_used,
			sym_set_count - sc, sym_set_count, t0
		);
	}
	prep_counts();
	write_pieces();
	teardown();
	return 0;
}
