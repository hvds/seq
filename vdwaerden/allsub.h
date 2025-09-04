#ifndef ALLSUB_H
#define ALLSUB_H

#include <stddef.h>
#include <stdint.h>

extern void init_allsub(void);
extern void find_allsub(uint n);

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

typedef struct allsub_s {
    sets_t *sets;
    size_t sets_size;
    size_t sets_count;
} allsub_t;

extern allsub_t allsub[];

#endif
