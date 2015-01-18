#ifndef INVERSE_H
#define INVERSE_H

/* printf */
#include <stdio.h>

/* malloc/calloc/realloc */
#include <stdlib.h> 

/* memset */
#include <string.h>

typedef unsigned int uint;
extern uint inveuclid(uint n, uint m);
extern void invtable(uint p, uint* t);
extern void setup_inverse(void);
extern void teardown_inverse(void);
extern void inverse_table(uint p);

extern uint** inverse;
extern uint inverse_size;

extern inline uint invfast(uint n, uint p) {
	return inverse[p][n];
}

#endif /* INVERSE_H */
