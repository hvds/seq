#ifndef SYM_SET_H
#define SYM_SET_H

#include "part.h"
#include "vec.h"

#ifndef IS_SYM_SET_C
#define SYM_SET_INLINE extern inline
#else
#define SYM_SET_INLINE inline
#endif

extern char* sym_set_arena;
typedef struct sym_set_s {
	unsigned int count;
	unsigned int index[0];
} sym_set_t;

sym_set_t* symset_new(void);
sym_set_t* symset_resize(sym_set_t* ss);
void symset_delete(sym_set_t* ss);

SYM_SET_INLINE sym_set_t* SSO(uint offset) {
	return (sym_set_t*)(sym_set_arena + offset);
}

SYM_SET_INLINE uint SSSIZE(sym_set_t* ss) {
	return sizeof(ss->count) + sizeof(uint) * ss->count;
}

#endif
