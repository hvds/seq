#define IS_VEC_C
#include "vec.h"

typedef unsigned int uint;
uchar transform[256];

void init_transform(void) {
	uint i, j;
	uchar value;
	for (i = 0; i < 256; ++i) {
		value = 0;
		for (j = 0; j < 8; ++j) {
			if (i & (1 << j))
				value |= 1 << (7 - j);
		}
		transform[i] = value;
	}
}

vech_tree* vech_new(void) {
	vech_tree* t = (vech_tree*)malloc(sizeof(vech_tree));

	t->vha_size = 100;
	t->vha_used = 0;
	t->vha_root = 0;
	t->vecharena = (vech_t*)malloc(t->vha_size * sizeof(vech_t));
	t->va_size = 100;
	t->va_used = 0;
	t->vecarena = (vec_t*)malloc(t->va_size * sizeof(vec_t));
	return t;
}

void vech_delete(vech_tree* t) {
	free(t->vecharena);
	free(t->vecarena);
	free(t);
}

vech_tree* vech_dup(vech_tree* source) {
	vech_tree* dest = (vech_tree*)malloc(sizeof(vech_tree));
	memcpy(dest, source, sizeof(vech_tree));
	dest->vecharena = (vech_t*)malloc(dest->vha_size * sizeof(vech_t));
	dest->vecarena = (vec_t*)malloc(dest->va_size * sizeof(vec_t));
	memcpy(dest->vecharena, source->vecharena, dest->vha_used * sizeof(vech_t));
	memcpy(dest->vecarena, source->vecarena, dest->va_used * sizeof(vec_t));
	return dest;
}

void resize_vech(vech_tree* t, uint size) {
	if (size <= t->vha_size)
		return;
	t->vha_size = t->vha_size * 3 / 2;
	if (size > t->vha_size)
		t->vha_size = size;
	t->vecharena = (vech_t*)realloc(t->vecharena, t->vha_size * sizeof(vech_t));
}

inline vech_t* VECH(vech_tree* t, vech index) {
	return (vech_t*)(t->vecharena + (index - 1));
}
inline vec_t* VEC(vech_tree* t, uint index) {
	return (vec_t*)(t->vecarena + index);
}
inline vec_t* AVEC(vech_tree* t, vech index) {
	return VEC(t, VECH(t, index)->ref);
}

uint vec_alloc(vech_tree* t, vec_t* v) {
	uint vi = t->va_used;
	++t->va_used;
	if (t->va_used > t->va_size) {
		t->va_size = t->va_size * 3 / 2;
		if (t->va_used > t->va_size)
			t->va_size = t->va_used;
		t->vecarena = (vec_t*)realloc(t->vecarena, t->va_size * sizeof(vec_t));
	}
	vec_copy(v, VEC(t, vi));
	return vi;
}

vech vech_alloc(vech_tree* t, vec_t* v) {
	vech ai = ++t->vha_used;	/* 1-based; 0 means not present */
	vech_t* ap;
	ap = VECH(t, ai);
	ap->left = 0;
	ap->right = 0;
	ap->balance = 0;
	ap->ref = vec_alloc(t, v);
	return ai;
}

void vech_rotateleft(vech_tree* t, vech* root) {
	vech a = *root;
	vech b = VECH(t, a)->right;
	*root = b;
	VECH(t, a)->right = VECH(t, b)->left;
	VECH(t, b)->left = a;
}

void vech_rotateright(vech_tree* t, vech* root) {
	vech a = *root;
	vech b = VECH(t, a)->left;
	*root = b;
	VECH(t, a)->left = VECH(t, b)->right;
	VECH(t, b)->right = a;
}

void vech_rebalance(vech_tree* t, vech base) {
	vech_t* bp = VECH(t, base);
	switch (bp->balance) {
	  case 0:
		VECH(t, bp->left)->balance = 0;
		VECH(t, bp->right)->balance = 0;
		break;
	  case -1:
		VECH(t, bp->left)->balance = 0;
		VECH(t, bp->right)->balance = 1;
		break;
	  case 1:
		VECH(t, bp->left)->balance = -1;
		VECH(t, bp->right)->balance = 0;
		break;
   }
   bp->balance = 0;
}

/* Insert an element into an VECH tree.
 * Returns VECH_GROW if the depth of the tree has grown.
 * Returns VECH_OK if the element was inserted without growing
 * Returns VECH_EXISTS if the element was already in the tree.
 */
vech_insert_t vechx_seen(vech_tree* t, vech* rootp, vec_t* v) {
	vech root = *rootp;	/* note: invalidated by any rotate */
	vech_t* rnp;

	if (!root) {
		*rootp = vech_alloc(t, v);
		return VECH_GROW;
	}
	rnp = VECH(t, root);

	switch (vec_cmp(VEC(t, rnp->ref), v)) {
	  case 1:
		/* insert into the left subtree */
		if (!rnp->left) {
			rnp->left = vech_alloc(t, v);
			if (rnp->balance--)
				return VECH_OK;
			return VECH_GROW;
		}
		switch (vechx_seen(t, &(rnp->left), v)) {
		  case VECH_GROW:
			switch (rnp->balance--) {
			  case 1:
				return VECH_OK;
			  case 0:
				return VECH_GROW;
			}
			if (VECH(t, rnp->left)->balance < 0) {
				vech_rotateright(t, rootp);
				rnp = VECH(t, *rootp);
				rnp->balance = 0;
				VECH(t, rnp->right)->balance = 0;
			} else {
				vech_rotateleft(t, &(rnp->left));
				vech_rotateright(t, rootp);
				vech_rebalance(t, *rootp);
			}
		  case VECH_OK:
			return VECH_OK;
		  case VECH_EXISTS:
			return VECH_EXISTS;
		}
	  case -1:
		/* insert into the right subtree */
		if (!rnp->right) {
			rnp->right = vech_alloc(t, v);
			if (rnp->balance++)
				return VECH_OK;
			return VECH_GROW;
		}
		switch (vechx_seen(t, &(rnp->right), v)) {
		  case VECH_GROW:
			switch (rnp->balance++) {
			  case -1:
				return VECH_OK;
			  case 0:
				return VECH_GROW;
			}
			if (VECH(t, rnp->right)->balance > 0) {
				vech_rotateleft(t, rootp);
				rnp = VECH(t, *rootp);
				rnp->balance = 0;
				VECH(t, rnp->left)->balance = 0;
			} else {
				vech_rotateright(t, &(rnp->right));
				vech_rotateleft(t, rootp);
				vech_rebalance(t, *rootp);
			}
		  case VECH_OK:
			return VECH_OK;
		  case VECH_EXISTS:
			return VECH_EXISTS;
		}
	  case 0:
		return VECH_EXISTS;
	}
}

vech_insert_t vech_seen(vech_tree* tree, vec_t* v) {
	/* do any realloc now, so the tree won't move under our feet */
	resize_vech(tree, tree->vha_used + 1);
	switch (vechx_seen(tree, &(tree->vha_root), v)) {
	  case VECH_OK:
	  case VECH_GROW:
		return VECH_OK;
	  case VECH_EXISTS:
		return VECH_EXISTS;
	}
}

uint bit_count[256];
void init_bit_count(void) {
	uint i;
	bit_count[0] = 0;
	for (i = 1; i < 256; ++i)
		if (i & 1)
			bit_count[i] = 1 + bit_count[i >> 1];
		else
			bit_count[i] = bit_count[i >> 1];
}

vec_t* connections;
void init_connections(void) {
	uint i, j;
	connections = (vec_t*)calloc(NODES, sizeof(vec_t));
	for (i = 0; i < NODES; ++i) {
		vec_t* v = connect_vec(i);
		for (j = 0; j < NBASE; ++j) {
			uint k = i ^ (1 << j);
			vec_setbit(v, k);
		}
	}
}
