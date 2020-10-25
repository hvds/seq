#ifndef SYM_H
#define SYM_H

#include "loc.h"

typedef unsigned char bool;

typedef enum {
    xy = 0,
    xY = 1,
    Xy = 2,
    XY = 3,
    yx = 4,
    yX = 5,
    Yx = 6,
    YX = 7
} sym_t;
#define MAXSYM 7

extern int bitcount[256];
typedef struct pack_set_s {
    int count;
    int set[70];
} pack_set_t;
extern pack_set_t pack_set[9];

extern void init_sym(void);
extern void finish_sym(void);
extern bool is_transpose(sym_t s);
extern bool is_reflect(sym_t s);
extern bool sym_dup(sym_t s, sym_t t);
extern bool sym_check(sym_t s, int x, int y, int *vals);
extern bool sym_checkloc(sym_t s, int x, int y, loc_t l);
extern sym_t sym_reflect(int syms, int x, int y, loc_t l);
extern int *sym_lookup(sym_t s);
extern int *sym_transform(sym_t s, int x, int y, int *vals);
extern loc_t sym_transloc(sym_t s, int x, int y, loc_t l);

#endif
