#ifndef _AVL_H
#define _AVL_H 1

typedef unsigned int avl;
typedef unsigned char uchar;
typedef int(avl_cmp_t)(uchar* a, uchar* b, uint size);

typedef struct avl_s {
	avl left;
	avl right;
	int balance;
	uint ref;
	uint size;
} avl_t;

typedef struct avl_tree_s {
	avl_t* avlarena;
	avl aa_size;
	avl aa_used;
	avl aa_root;
	uchar* vecarena;
	uint va_size;
	uint va_used;
	avl_cmp_t* comparator;
} avl_tree;

avl_tree* avl_new(avl_cmp_t* comparator);
void avl_delete(avl_tree* tree);
avl_tree* avl_dup(avl_tree* source);
int avl_seen(avl_tree* tree, uchar* v, uint size);

#endif /* avl.h */