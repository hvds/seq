#ifndef MBH_H
#define MBH_H

#include <stdint.h>

typedef size_t bhsize_t;    /* must be same size as void* */

extern bhsize_t* mbharena;
typedef struct bh_s_heap {
	bhsize_t size;
	void* context;
	void* heap[0];
} bh_heap;
typedef bhsize_t bhp;
#define BHP(h) ((bh_heap*)&mbharena[h])

extern void mbh_init(void);
extern void mbh_final(void);
extern bhp mbh_new(void* context);
extern void mbh_delete(bhp h);
extern void mbh_insert(bhp h, void* v);
extern void* mbh_shift(bhp h);
extern inline bhsize_t mbh_size(bhp h) {
	return BHP(h)->size;
}
extern int mbh_compare(bhp h, bhsize_t left, bhsize_t right);

#define P2I(x) (int)(intptr_t)(x)
#define I2P(x) (void*)(intptr_t)(x)

#endif /* MBH_H */
