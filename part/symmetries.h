#ifndef SYMMETRIES_H
#define SYMMETRIES_H

#include "vec.h"

typedef struct sym_s {
	unsigned int map[NODES];
} sym_t;

extern unsigned int sym_count;
extern sym_t* symmetries;

void setup_symmetries(void);
void teardown_symmetries(void);

#ifndef IS_SYMMETRIES_C
#define SYM_INLINE extern inline
#else
#define SYM_INLINE inline
#endif

SYM_INLINE sym_t* sym_map(unsigned int i) {
	return symmetries + i;
}

SYM_INLINE void apply_map2(sym_t* sym, vec_t* src, vec_t* dest) {
	unsigned int i;
	for (i = 0; i < NODES; ++i) {
		if (vec_testbit(src, sym->map[i]))
			vec_setbit(dest, i); 
		else
			vec_clearbit(dest, i); 
	}
}   

#endif
