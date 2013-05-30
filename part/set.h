#ifndef _SETH_H
#define _SETH_H 1

#include "part.h"
#include "vec.h"

typedef unsigned int seth;
typedef int(seth_cmp_t)(vec_t* a, vec_t* b, uint size);

typedef struct set_s {
	vec_t s[NODES];
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
	uchar* setarena;
	uint sa_size;
	uint sa_used;
	uint element_size;
} seth_tree;

seth_tree* seth_new(uint element_size);
void seth_delete(seth_tree* tree);
seth_tree* seth_dup(seth_tree* source);
seth_insert_t seth_seen(seth_tree* tree, set_t* v);

#endif /* seth.h */
