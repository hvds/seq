#ifndef _AVL_H
#define _AVL_H 1

typedef unsigned int avl;
typedef unsigned char uchar;
typedef int(avl_cmp_t)(uchar* a, uchar* b, uint size);

typedef enum {
	AVL_GROW = -1,	/* insert ok, tree depth has grown */
	AVL_OK = 0,		/* insert ok, tree depth has not grown */
	AVL_EXISTS = 1	/* no insert, element already exists */
} avl_insert_t;

typedef struct avl_s {
	avl left;
	avl right;
	int balance;
	uint ref;
} avl_t;

typedef struct avl_tree_s {
	avl_t* avlarena;
	avl aa_size;
	avl aa_used;
	avl aa_root;
	uchar* vecarena;
	uint va_size;
	uint va_used;
	uint element_size;
} avl_tree;

avl_tree* avl_new(uint element_size);
void avl_delete(avl_tree* tree);
avl_tree* avl_dup(avl_tree* source);
avl_insert_t avl_seen(avl_tree* tree, uchar* v);

#endif /* avl.h */
