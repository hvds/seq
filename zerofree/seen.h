#ifndef SEEN_H

#include <gmp.h>

extern void init_seen(void);
extern void free_seen(void);
extern int seen(mpz_t z);
extern void dump_seen(void);
extern void max_seen(mpz_t z);

#endif
