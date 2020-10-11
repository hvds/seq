#ifndef BOARD_H
#define BOARD_H

#include "group.h"

typedef struct board_s {
    int k;
    int unused;
    int refcount;
    int groups;
    group_t *group[2];
} board_t;

extern int best_k;
extern board_t *best_board;
extern unsigned long board_count;

extern board_t *init_board(int n, int feedback);
extern void finish_board(void);
extern board_t *new_board(int k, int unused, group_t *g0, group_t *g1);
extern void try_board(board_t *b);
extern void print_board(board_t *b);

#endif
