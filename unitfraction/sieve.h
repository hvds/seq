#ifndef SIEVE_H
#define SIEVE_H

typedef unsigned int uint;
typedef unsigned char bool;

extern void extend_sieve(uint limit);
extern uint *prime;
extern uint nprime;
extern uint maxsieve;

#endif
