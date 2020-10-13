#ifndef SYM_H
#define SYM_H

#include "group.h"

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

extern void init_sym(void);
extern bool is_transpose(sym_t s);
extern bool sym_check(sym_t s, int x, int y, int *vals);
extern bool sym_checkloc(sym_t s, int x, int y, loc_t l);
extern int *sym_transform(sym_t s, int x, int y, int *vals);
extern loc_t sym_transloc(sym_t s, int x, int y, loc_t l);

#endif
