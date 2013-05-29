#include "set.h"

typedef unsigned int uint;

int set_comparator(set_t* s1, set_t* s2, uint size) {
	uint i;
	signed int c;
	uchar* str1 = (uchar*)&(s1[0].s[0].v);
	uchar* str2 = (uchar*)&(s2[0].s[0].v);
	for (i = 0; i < size; ++i) {
		c = transform[str1[i]] - transform[str2[i]];
		if (c > 0)
			return 1;
		if (c < 0)
			return -1;
	}
	return 0;
}

seth_tree* seth_new(uint element_size) {
	seth_tree* t = (seth_tree*)malloc(sizeof(seth_tree));

	t->sha_size = 100;
	t->sha_used = 0;
	t->sha_root = 0;
	t->setharena = (seth_t*)malloc(t->sha_size * sizeof(seth_t));
	t->sa_size = 100;
	t->sa_used = 0;
	t->setarena = (uchar*)malloc(t->sa_size * element_size * VECSIZE);
	t->element_size = element_size;
	return t;
}

void seth_delete(seth_tree* t) {
	free(t->setharena);
	free(t->setarena);
	free(t);
}

seth_tree* seth_dup(seth_tree* source) {
	seth_tree* dest = (seth_tree*)malloc(sizeof(seth_tree));
	memcpy(dest, source, sizeof(seth_tree));
	dest->setharena = (seth_t*)malloc(dest->sha_size * sizeof(seth_t));
	dest->setarena = (uchar*)malloc(dest->sa_size * dest->element_size * VECSIZE);
	memcpy(dest->setharena, source->setharena, dest->sha_used * sizeof(seth_t));
	memcpy(dest->setarena, source->setarena, dest->sa_used * dest->element_size * VECSIZE);
	return dest;
}

void resize_seth(seth_tree* t, uint size) {
	if (size <= t->sha_size)
		return;
	t->sha_size = t->sha_size * 3 / 2;
	if (size > t->sha_size)
		t->sha_size = size;
	t->setharena = (seth_t*)realloc(t->setharena, t->sha_size * sizeof(seth_t));
}

inline seth_t* SETH(seth_tree* t, seth index) {
	return (seth_t*)(t->setharena + (index - 1));
}
inline set_t* SET(seth_tree* t, uint index) {
	return (set_t*)(t->setarena + index * t->element_size * VECSIZE);
}
inline set_t* ASET(seth_tree* t, seth index) {
	return SET(t, SETH(t, index)->ref);
}

uint set_alloc(seth_tree* t, set_t* v) {
	uint vi = t->sa_used;
	++t->sa_used;
	if (t->sa_used > t->sa_size) {
		t->sa_size = t->sa_size * 3 / 2;
		if (t->sa_used > t->sa_size)
			t->sa_size = t->sa_used;
		t->setarena = (uchar*)realloc(t->setarena, t->sa_size * t->element_size * VECSIZE);
	}
	memcpy(SET(t, vi), v, t->element_size * VECSIZE);
	return vi;
}

seth seth_alloc(seth_tree* t, set_t* v) {
	seth ai = ++t->sha_used;	/* 1-based; 0 means not present */
	seth_t* ap;
	ap = SETH(t, ai);
	ap->left = 0;
	ap->right = 0;
	ap->balance = 0;
	ap->ref = set_alloc(t, v);
	return ai;
}

void seth_rotateleft(seth_tree* t, seth* root) {
	seth a = *root;
	seth b = SETH(t, a)->right;
	*root = b;
	SETH(t, a)->right = SETH(t, b)->left;
	SETH(t, b)->left = a;
}

void seth_rotateright(seth_tree* t, seth* root) {
	seth a = *root;
	seth b = SETH(t, a)->left;
	*root = b;
	SETH(t, a)->left = SETH(t, b)->right;
	SETH(t, b)->right = a;
}

void seth_rebalance(seth_tree* t, seth base) {
	seth_t* bp = SETH(t, base);
	switch (bp->balance) {
	  case 0:
		SETH(t, bp->left)->balance = 0;
		SETH(t, bp->right)->balance = 0;
		break;
	  case -1:
		SETH(t, bp->left)->balance = 0;
		SETH(t, bp->right)->balance = 1;
		break;
	  case 1:
		SETH(t, bp->left)->balance = -1;
		SETH(t, bp->right)->balance = 0;
		break;
   }
   bp->balance = 0;
}

/* Insert an element into an SETH tree.
 * Returns 1 if the depth of the tree has grown.
 */
seth_insert_t sethx_seen(seth_tree* t, seth* rootp, set_t* v) {
	seth root = *rootp;	/* note: invalidated by any rotate */
	seth_t* rnp;

	if (!root) {
		*rootp = seth_alloc(t, v);
		return SETH_GROW;
	}
	rnp = SETH(t, root);

	switch (set_comparator(SET(t, rnp->ref), v, t->element_size * VECSIZE)) {
	  case 1:
		/* insert into the left subtree */
		if (!rnp->left) {
			rnp->left = seth_alloc(t, v);
			if (rnp->balance--)
				return SETH_OK;
			return SETH_GROW;
		}
		switch (sethx_seen(t, &(rnp->left), v)) {
		  case SETH_GROW:
			switch (rnp->balance--) {
			  case 1:
				return SETH_OK;
			  case 0:
				return SETH_GROW;
			}
			if (SETH(t, rnp->left)->balance < 0) {
				seth_rotateright(t, rootp);
				rnp = SETH(t, *rootp);
				rnp->balance = 0;
				SETH(t, rnp->right)->balance = 0;
			} else {
				seth_rotateleft(t, &(rnp->left));
				seth_rotateright(t, rootp);
				seth_rebalance(t, *rootp);
			}
		  case SETH_OK:
			return SETH_OK;
		  case SETH_EXISTS:
			return SETH_EXISTS;
		}
	  case -1:
		/* insert into the right subtree */
		if (!rnp->right) {
			rnp->right = seth_alloc(t, v);
			if (rnp->balance++)
				return SETH_OK;
			return SETH_GROW;
		}
		switch (sethx_seen(t, &(rnp->right), v)) {
		  case SETH_GROW:
			switch (rnp->balance++) {
			  case -1:
				return SETH_OK;
			  case 0:
				return SETH_GROW;
			}
			if (SETH(t, rnp->right)->balance > 0) {
				seth_rotateleft(t, rootp);
				rnp = SETH(t, *rootp);
				rnp->balance = 0;
				SETH(t, rnp->left)->balance = 0;
			} else {
				seth_rotateright(t, &(rnp->right));
				seth_rotateleft(t, rootp);
				seth_rebalance(t, *rootp);
			}
		  case SETH_OK:
			return SETH_OK;
		  case SETH_EXISTS:
			return SETH_EXISTS;
		}
	  case 0:
		return SETH_EXISTS;
	}
}

seth_insert_t seth_seen(seth_tree* tree, set_t* v) {
	/* do any realloc now, so the tree won't move under our feet */
	resize_seth(tree, tree->sha_used + 1);
	switch (sethx_seen(tree, &(tree->sha_root), v)) {
	  case SETH_OK:
	  case SETH_GROW:
		return SETH_OK;
	  case SETH_EXISTS:
		return SETH_EXISTS;
	}
}

