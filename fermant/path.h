#ifndef PATH_H
#define PATH_H

#include <string.h>

#include "types.h"

/* a vector of nv bits representing the grid lines comprising a path */
typedef uint path_t;

/* a vector of npaths bits representing the paths available in a fragment */
typedef uint pathset_t;

/* a pair of path indices representing the process of resolving which path
 * is shorter in a given context */
typedef struct {
    uint pi;
    uint pj;
} resolve_t;

extern uint npaths;     /* the number of paths */
extern path_t *paths;
extern uint sizepaths;
extern uint nresolve;   /* the number of resolutions */

extern void init_paths(int strategy);
extern void done_paths(void);
extern uint render_path(char *buf, uint buflen, uint pi);
/* takes (const path_t *) */
extern int path_comparator(const void *a, const void *b);
extern pathset_t all_paths(void);
extern uint split_all(uint recover);

__inline path_t path_p(uint pi) {
    return paths[pi];
}

/* Return a new path_t that adds the line segment represented by the
 * specified var to the provided path.
 */
__inline path_t include_line(path_t cur, uint var) {
    return cur | (1 << var);
}

/* Return index of first variable included in a path_t, or undefined if none.
 */
__inline uint path_first(path_t p) {
    return ffs((int)p) - 1;
}

__inline path_t path_xor(uint pi, uint pj) {
    return paths[pi] ^ paths[pj];
}

/* Return the number of variables included in a path_t
 */
__inline uint path_count(path_t p) {
    return __builtin_popcount(p);
}

/* Return index of first path included in a pathset_t, or undefined if none.
 */
__inline path_t pathset_first(pathset_t ps) {
    return ffs((int)ps) - 1;
}

/* Return the number of paths included in a pathset_t
 */
__inline uint pathset_count(pathset_t ps) {
    return __builtin_popcount(ps);
}

#endif /* PATH_H */
