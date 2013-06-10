#ifndef _SETH_H
#define _SETH_H 1

#include <stdio.h>
#include <assert.h>
#include "part.h"
#include "vec.h"

typedef unsigned int seth;

typedef struct set_s {
	unsigned char p[NODES];
} set_t;

typedef enum {
	SETH_GROW = -1,		/* insert ok, tree depth has grown */
	SETH_OK = 0,		/* insert ok, tree depth has not grown */
	SETH_EXISTS = 1		/* no insert, element already exists */
} seth_insert_t;

typedef struct seth_s {
	seth left;
	seth right;
	int balance;
	uint ref;
} seth_t;

typedef struct seth_tree_s {
	seth_t* setharena;
	seth sha_size;
	seth sha_used;
	seth sha_root;
	set_t* setarena;
	uint sa_size;
	uint sa_used;
} seth_tree;

seth_tree* seth_new(void);
void seth_delete(seth_tree* tree);
void seth_reset(seth_tree* tree);
seth_tree* seth_dup(seth_tree* source);
seth_insert_t seth_seen(seth_tree* tree, set_t* v);

#ifndef IS_SET_C
#define SET_INLINE extern inline
#else
#define SET_INLINE inline
#endif

SET_INLINE void set_copy(set_t* source, set_t* dest) {
	memcpy(dest, source, sizeof(set_t));
}

SET_INLINE void set_init(set_t* set, vec_t* v) {
	uint i;
	for (i = 0; i < NODES; ++i)
		set->p[i] = vec_testbit(v, i);
}

SET_INLINE void set_append(set_t* dest, set_t* source, vec_t* v, uint size) {
	uint i;
	set_copy(source, dest);
	for (i = 0; i < NODES; ++i)
		if (vec_testbit(v, i)) {
			assert(dest->p[i] == 0);
			dest->p[i] = size + 1;
		}
}

void fprint_set(FILE* stream, set_t* set, uint size);

#endif /* seth.h */
