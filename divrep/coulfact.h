#ifndef COULFACT_H
#define COULFACT_H 1
#include <gmp.h>

typedef struct s_ppow {
    ulong p;
    uint e;
} t_ppow;
typedef struct s_fact {
    uint count;
    uint size;
    t_ppow *ppow;
} t_fact;
typedef struct s_zpow {
    mpz_t p;
    uint e;
} t_zpow;
typedef struct s_zfact {
    uint count;
    uint size;
    t_zpow *ppow;
} t_zfact;

extern void init_fact(t_fact *f);
extern void free_fact(t_fact *f);
extern void add_fact(t_fact *f, t_ppow pp);
extern void reverse_fact(t_fact *f);
extern void init_zfact(t_zfact *f);
extern void free_zfact(t_zfact *f);
extern void add_zfact(t_zfact *f, t_zpow pp);
extern void simple_fact(uint n, t_fact *f);
extern uint simple_tau(t_fact *f);
extern uint simple_valuation(ulong n, ulong p);
extern uint simple_prime_count(ulong n);
extern uint tiny_gcd(uint a, uint b);
extern ulong simple_gcd(ulong a, ulong b);

#endif
