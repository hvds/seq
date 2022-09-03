#ifndef ROOTMOD_H
#define ROOTMOD_H 1

#include <gmp.h>
#include "factor.h"

typedef struct s_results {
    uint size;
    uint count;
    mpz_t *r;
} t_results;

extern void init_rootmod(uint levels);
extern void done_rootmod(void);
extern void allrootmod(uint level, mpz_t a, uint k, mpz_t n_factors);
extern t_results *res_array(uint level);
extern void res_copy(uint new_level, uint old_level);
extern void root_extract(uint new_level, uint old_level, uint k, mpz_t n);
extern void root_extend(uint new_level, uint old_level, mpz_t n,
        mpz_t a, uint k, ulong p, uint e, mpz_t px);

#endif
