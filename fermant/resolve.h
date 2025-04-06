#ifndef RESOLVE_H
#define RESOLVE_H

#include "types.h"

/* Struct representing the resolution of paths[pi] against paths[pj] */
typedef struct {
    uint pi;
    uint pj;
} resolve_t;

extern uint nresolve;

extern void init_resolve(void);
extern void done_resolve(void);
extern void split_all(void);

#endif /* RESOLVE_H */
