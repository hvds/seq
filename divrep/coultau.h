/* vim: set et ts=2 sw=2 sts=2: */
#ifndef COULTAU_H
#define COULTAU_H

#include <gmp.h>
#include "ptypes.h"

typedef enum {
  FS_INIT = 0,
  FS_TRIAL,
  FS_POWER,
  FS_LARGE,
  FS_TERM,
} fs_state_t;

#define MAX_FACTORS 128

typedef struct factor_state_s {
  fs_state_t state;
  mpz_t n;  /* remaining number to be factored */
  mpz_t f;  /* new factor found */
  int e;    /* new exponent found */
  int ef;   /* exponent multiplier */
  UV tlim;  /* p^2 limit checked by trial division */

  /* used only for trial division phase */
  UV sp;  /* smallprime index */
  UV un;

  /* used only after trial division phase */
  int log;    /* verbose_level */
  int ntofac; /* number of additional factors in tofac_stack[] */
  mpz_t tofac_stack[MAX_FACTORS];
} factor_state;

extern void init_tau(void);
extern int factor_one(factor_state* fs);
extern int is_taux(mpz_t n, uint32_t k, uint32_t x);

#endif
