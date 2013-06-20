#ifndef PIECES_H
#define PIECES_H

#include "part.h"
#include "vec.h"
#include "sym_set.h"

#ifndef IS_PIECES_C
#define PIECE_INLINE extern inline
#else
#define PIECE_INLINE inline
#endif

typedef struct piece_s {
	vec_t v;
	uint sso;	/* handle to the symmetries that yield distinct shapes */
} piece_t;

extern uint piece_array[];
extern piece_t* pieces;

PIECE_INLINE piece_t* pieces_vec(uint index) {
	return pieces + index;
}
PIECE_INLINE piece_t* pieces_for(uint size) {
	return pieces_vec(piece_array[size]);
}

extern void setup_pieces(void);
extern void teardown_pieces(void);

#endif
