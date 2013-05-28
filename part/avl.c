#include <stdlib.h>
#include <string.h>
#include "avl.h"

typedef unsigned int uint;

avl_tree* avl_new(uint element_size) {
	avl_tree* t = (avl_tree*)malloc(sizeof(avl_tree));

	t->aa_size = 100;
	t->aa_used = 0;
	t->aa_root = 0;
	t->avlarena = (avl_t*)malloc(t->aa_size * sizeof(avl_t));
	t->va_size = 100;
	t->va_used = 0;
	t->vecarena = (uchar*)malloc(t->va_size * element_size);
	t->element_size = element_size;
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
	dest->vecarena = (uchar*)malloc(dest->va_size * dest->element_size);
	memcpy(dest->avlarena, source->avlarena, dest->aa_used * sizeof(avl_t));
	memcpy(dest->vecarena, source->vecarena, dest->va_used * dest->element_size);
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
	return (uchar*)(t->vecarena + index * t->element_size);
}
inline uchar* AVEC(avl_tree* t, avl index) {
	return VEC(t, AVL(t, index)->ref);
}

int avl_cmp(avl_tree* t, uchar* v1, uchar* v2) {
	return set_comparator(v1, v2, t->element_size);
}

uint vec_alloc(avl_tree* t, uchar* v) {
	uint vi = t->va_used;
	++t->va_used;
	if (t->va_used > t->va_size) {
		t->va_size = t->va_size * 3 / 2;
		if (t->va_used > t->va_size)
			t->va_size = t->va_used;
		t->vecarena = (uchar*)realloc(t->vecarena, t->va_size * t->element_size);
	}
	memcpy(VEC(t, vi), v, t->element_size);
	return vi;
}

avl avl_alloc(avl_tree* t, uchar* v) {
	avl ai = ++t->aa_used;	/* 1-based; 0 means not present */
	avl_t* ap;
	ap = AVL(t, ai);
	ap->left = 0;
	ap->right = 0;
	ap->balance = 0;
	ap->ref = vec_alloc(t, v);
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
avl_insert_t avlx_seen(avl_tree* t, avl* rootp, uchar* v) {
	avl root = *rootp;	/* note: invalidated by any rotate */
	avl_t* rnp;

	if (!root) {
		*rootp = avl_alloc(t, v);
		return AVL_GROW;
	}
	rnp = AVL(t, root);

	switch (avl_cmp(t, VEC(t, rnp->ref), v)) {
	  case 1:
		/* insert into the left subtree */
		if (!rnp->left) {
			rnp->left = avl_alloc(t, v);
			if (rnp->balance--)
				return AVL_OK;
			return AVL_GROW;
		}
		switch (avlx_seen(t, &(rnp->left), v)) {
		  case AVL_GROW:
			switch (rnp->balance--) {
			  case 1:
				return AVL_OK;
			  case 0:
				return AVL_GROW;
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
		  case AVL_OK:
			return AVL_OK;
		  case AVL_EXISTS:
			return AVL_EXISTS;
		}
	  case -1:
		/* insert into the right subtree */
		if (!rnp->right) {
			rnp->right = avl_alloc(t, v);
			if (rnp->balance++)
				return AVL_OK;
			return AVL_GROW;
		}
		switch (avlx_seen(t, &(rnp->right), v)) {
		  case AVL_GROW:
			switch (rnp->balance++) {
			  case -1:
				return AVL_OK;
			  case 0:
				return AVL_GROW;
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
		  case AVL_OK:
			return AVL_OK;
		  case AVL_EXISTS:
			return AVL_EXISTS;
		}
	  case 0:
		return AVL_EXISTS;
	}
}

avl_insert_t avl_seen(avl_tree* tree, uchar* v) {
	/* do any realloc now, so the tree won't move under our feet */
	resize_avl(tree, tree->aa_used + 1);
	switch (avlx_seen(tree, &(tree->aa_root), v)) {
	  case AVL_OK:
	  case AVL_GROW:
		return AVL_OK;
	  case AVL_EXISTS:
		return AVL_EXISTS;
	}
}

