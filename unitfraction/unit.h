#ifndef UNIT_H
#define UNIT_H

#include <gmp.h>

typedef unsigned char bool;
typedef unsigned int uint;

/* maxdepth is the greatest cardinality we expect to need. */
extern void init_unit(uint maxdepth);
extern void done_unit(void);

/* Returns TRUE if there exists S: |S| < c, sum{1/s_i} = q. */
extern bool better_set(mpq_t q, uint c);
/* Returns TRUE if there exists M: |M| < c, sum{1/m_i} = q. */
extern bool better_multi(mpq_t q, uint c);

/* Returns min(| S |): sum{1/s_i} = q. */
extern uint find_set(mpq_t q);
/* Returns min(| M |): sum{1/m_i} = q. */
extern uint find_multi(mpq_t q);
/* Returns min(| S |): sum{1/s_i^2} = q. */
extern uint find_square_set(mpq_t q, uint min_depth, uint max_depth);

/* handy macros in case we need to debug our GMP memory use */
#define ZINIT(z) mpz_init(z)
#define ZCLEAR(z) mpz_clear(z)
#define QINIT(q) mpq_init(q)
#define QCLEAR(q) mpq_clear(q)

#endif
