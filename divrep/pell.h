#ifndef PELL_H
#define PELL_H 1

#include <gmp.h>
#include "types.h"

typedef unsigned char bool;

extern void init_pell(void);
extern void done_pell(void);
/* solve Ax^2 - By^2 + D = 0, with 0 < x <= limit */
extern void new_pell(mpz_t A, mpz_t B, int D, mpz_t limit);
extern bool next_pell(mpz_t x, mpz_t y);

#endif
