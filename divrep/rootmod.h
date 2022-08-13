#ifndef ROOTMOD_H
#define ROOTMOD_H 1

#include <gmp.h>
#include "factor.h"

extern void init_rootmod(void);
extern void done_rootmod(void);
extern uint allrootmod(mpz_t **r, mpz_t a, uint k, mpz_t n_factors);

#endif
