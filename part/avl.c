#include <stdlib.h>
#include <string.h>
#include "avl.h"

typedef unsigned int uint;

avl_tree* avl_new(avl_cmp_t* comparator) {
	avl_tree* t = (avl_tree*)malloc(sizeof(avl_tree));

	t->aa_size = 100;
	t->aa_used = 0;
	t->aa_root = 0;
	t->avlarena = (avl_t*)malloc(t->aa_size * sizeof(avl_t));
	t->va_size = 1024;
	t->va_used = 0;
	t->vecarena = (uchar*)malloc(t->va_size);
	t->comparator = comparator;
	return t;
}

void avl_delete(avl_tree* t) {
	free(t->avlarena);
	free(t->vecarena);
	free(t);
}

avl_tree* avl_dup(avl_tree* source) {
	avl_tree* dest = (avl_tree*)malloc(sizeof(avl_tree));
	memcpy(dest, source, sizeof(avl_tree));
	dest->avlarena = (avl_t*)malloc(dest->aa_size * sizeof(avl_t));
	dest->vecarena = (uchar*)malloc(dest->va_size);
	memcpy(dest->avlarena, source->avlarena, dest->aa_used * sizeof(avl_t));
	memcpy(dest->vecarena, source->vecarena, dest->va_used);
	return dest;
}

void resize_avl(avl_tree* t, uint size) {
	if (size <= t->aa_size)
		return;
	t->aa_size = t->aa_size * 3 / 2;
	if (size > t->aa_size)
		t->aa_size = size;
	t->avlarena = (avl_t*)realloc(t->avlarena, t->aa_size * sizeof(avl_t));
}

inline avl_t* AVL(avl_tree* t, avl index) {
	return (avl_t*)(t->avlarena + (index - 1));
}
inline uchar* VEC(avl_tree* t, uint index) {
	return (uchar*)(t->vecarena + index);
}
inline uchar* AVEC(avl_tree* t, avl index) {
	return VEC(t, AVL(t, index)->ref);
}

int avl_cmp(avl_tree* t, uchar* v1, uint size1, uchar* v2, uint size2) {
	if (size1 < size2)
		return -1;
	if (size1 > size2)
		return 1;
	return (*(t->comparator))(v1, v2, size1);
}

uint vec_alloc(avl_tree* t, uchar* v, uint size) {
	uint vi = t->va_used;
	t->va_used += size;
	if (t->va_used > t->va_size) {
		t->va_size = t->va_size * 3 / 2;
		if (t->va_used > t->va_size)
			t->va_size = t->va_used;
		t->vecarena = (uchar*)realloc(t->vecarena, t->va_size);
	}
	memcpy(VEC(t, vi), v, size);
	return vi;
}

avl avl_alloc(avl_tree* t, uint ref, uint size) {
	avl ai = ++t->aa_used;	/* 1-based; 0 means not present */
	avl_t* ap;
	resize_avl(t, t->aa_used);
	ap = AVL(t, ai);
	ap->left = 0;
	ap->right = 0;
	ap->balance = 0;
	ap->ref = ref;
	ap->size = size;
	return ai;
}

void avl_rotateleft(avl_tree* t, avl* root) {
	avl a = *root;
	avl b = AVL(t, a)->right;
	*root = b;
	AVL(t, a)->right = AVL(t, b)->left;
	AVL(t, b)->left = a;
}

void avl_rotateright(avl_tree* t, avl* root) {
	avl a = *root;
	avl b = AVL(t, a)->left;
	*root = b;
	AVL(t, a)->left = AVL(t, b)->right;
	AVL(t, b)->right = a;
}

void avl_rebalance(avl_tree* t, avl base) {
	avl_t* bp = AVL(t, base);
	switch (bp->balance) {
	  case 0:
		AVL(t, bp->left)->balance = 0;
		AVL(t, bp->right)->balance = 0;
		break;
	  case -1:
		AVL(t, bp->left)->balance = 0;
		AVL(t, bp->right)->balance = 1;
		break;
	  case 1:
		AVL(t, bp->left)->balance = -1;
		AVL(t, bp->right)->balance = 0;
		break;
   }
   bp->balance = 0;
}

/* Insert an element into an AVL tree.
 * Returns 1 if the depth of the tree has grown.
 */
int avlx_insert(avl_tree* t, avl* rootp, avl element) {
	avl root = *rootp;	/* note: invalidated by any rotate */
	avl_t* elp = AVL(t, element);
	avl_t* rnp;

	if (!root) {
		*rootp = element;
		return 1;
	}
	rnp = AVL(t, root);

	if (avl_cmp(
		t, VEC(t, rnp->ref), rnp->size, VEC(t, elp->ref), elp->size
	) > 0) {
		/* insert into the left subtree */
		if (!rnp->left) {
			rnp->left = element;
			if (rnp->balance--)
				return 0;
			return 1;
		}
		if (avlx_insert(t, &(rnp->left), element)) {
			switch (rnp->balance--) {
			  case 1:
				return 0;
			  case 0:
				return 1;
			}
			if (AVL(t, rnp->left)->balance < 0) {
				avl_rotateright(t, rootp);
				rnp = AVL(t, *rootp);
				rnp->balance = 0;
				AVL(t, rnp->right)->balance = 0;
			} else {
				avl_rotateleft(t, &(rnp->left));
				avl_rotateright(t, rootp);
				avl_rebalance(t, *rootp);
			}
		}
		return 0;
	} else {
		/* insert into the right subtree */
		if (!rnp->right) {
			rnp->right = element;
			if (rnp->balance++)
				return 0;
			return 1;
		}
		if (avlx_insert(t, &(rnp->right), element)) {
			switch (rnp->balance++) {
			  case -1:
				return 0;
			  case 0:
				return 1;
			}
			if (AVL(t, rnp->right)->balance > 0) {
				avl_rotateleft(t, rootp);
				rnp = AVL(t, *rootp);
				rnp->balance = 0;
				AVL(t, rnp->left)->balance = 0;
			} else {
				avl_rotateright(t, &(rnp->right));
				avl_rotateleft(t, rootp);
				avl_rebalance(t, *rootp);
			}
		}
		return 0;
	}
}

int avlx_exists(avl_tree* t, avl node, uchar* v, uint size) {
	avl_t* np;
	if (!node)
		return 0;
	np = AVL(t, node);
	switch (avl_cmp(t, VEC(t, np->ref), np->size, v, size)) {
	  case 0:
		return 1;
	  case 1:
		return avlx_exists(t, np->left, v, size);
	  case -1:
		return avlx_exists(t, np->right, v, size);
	}
}

void avl_insert(avl_tree* tree, avl element) {
	avlx_insert(tree, &(tree->aa_root), element);
}

int avl_exists(avl_tree* tree, uchar* v, uint size) {
	avl root = tree->aa_root;
	if (!root)
		return 0;
	return avlx_exists(tree, root, v, size);
}

int avl_seen(avl_tree* tree, uchar* v, uint size) {
	avl element;
	uint vi;

	if (avl_exists(tree, v, size))
		return 1;
	vi = vec_alloc(tree, v, size);
	element = avl_alloc(tree, vi, size);
	avl_insert(tree, element);
	return 0;
}

