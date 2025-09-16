#ifndef ALLSUB_H
#define ALLSUB_H

#include <stddef.h>
#include <stdint.h>

/* for calculating A051013() we will need 64 bits, but as a helper for cf
 * 32 bits is sufficient.
 */
#ifndef MAX_ALLSUB
#   define MAX_ALLSUB 32
#endif
#if MAX_ALLSUB > 32
    typedef uint64_t sets_t;
#else
    typedef uint32_t sets_t;
#endif

extern void init_allsub(void);
extern void find_allsub(uint n);
extern sets_t count_allsub(uint n);

typedef struct allsub_iter_s allsub_iter_t;

extern allsub_iter_t *init_iter(uint n);
extern void done_iter(allsub_iter_t *aip);
extern sets_t do_iter(allsub_iter_t *aip);

typedef struct allsub_s {
    sets_t *sets;
    size_t sets_size;
    size_t sets_count;
} allsub_t;

extern allsub_t allsub[];

#endif
