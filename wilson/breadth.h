#ifndef BREADTH_H
#define BREADTH_H

#include "mygmp.h"

/* The two arrays below account for the bulk of the memory use. I don't
 * ever want them to exceed 48GB.
 */
#define HARD_LIMIT (48UL * (1 << 30UL) / sizeof(pack_t))

typedef struct rat_array_s {
    pack_t* space;
    size_t size;    /* malloced size, in sizeof(pack_t) */
    size_t count;   /* space used, in sizeof(pack_t) */
    size_t actual;  /* number of values stored */
} rat_array_t;

extern rat_array_t ra1, ra2;

/* The generations alternate, pick values from one array and pushing new
 * values to the other array.
 */
extern inline rat_array_t *choose_cur(ulong gen) {
    return (gen & 1) ? &ra1 : &ra2;
}

extern inline rat_array_t *choose_next(ulong gen) {
    return (gen & 1) ? &ra2 : &ra1;
}

extern void init_breadth(mpq_t r);
extern void finish_breadth(void);

/* Given the pending queue for generation I<gen> in the C<choose_cur(gen)>
 * array, check that list for solutions, pushing any new values to check
 * onto the C<choose_next(gen)> array.
 * Returns boolean TRUE if a solution was found (in which case the contents
 * of the C<choose_next(gen)> array may be incomplete), else FALSE.
 */
bool breadth_one(ulong gen);

#endif
