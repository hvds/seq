#include <stdlib.h>
#include <stdio.h>

#include "board.h"

int n;
int feedback;
int best_k;
board_t *best_board;
unsigned long board_count;
board_t *b0;

void ref_board(board_t *b) {
    ++b->refcount;
}

void unref_board(board_t *b) {
    if (--b->refcount > 0)
        return;

    for (int i = 0; i < b->groups; ++i) {
        unref_group(b->group[i]);
    }
    free(b);
}

board_t *init_board(int n0, int feedback0) {
    n = n0;
    feedback = feedback0;
    board_count = 0;

    b0 = new_board(2, n, (group_t *)NULL, (group_t *)NULL);
    best_k = 1;
    best_board = b0;
    ref_board(b0);

    return b0;
}

void finish_board(void) {
    unref_board(best_board);
    unref_board(b0);
}

board_t *new_board(int k, int unused, group_t *g0, group_t *g1) {
    board_t *b = malloc(sizeof(board_t));

    b->k = k;
    b->unused = unused;
    b->refcount = 1;
    b->group[0] = g0;
    b->group[1] = g1;
    b->groups = g0 ? g1 ? 2 : 1 : 0;
    for (int i = 0; i < b->groups; ++i) {
        ref_group(b->group[i]);
    }
    return b;
}

void print_board(board_t *b) {
    printf("b=%lu; k=%d; u=%d; g=%d ->  ",
            board_count, b->k, b->unused, b->groups);
    for (int i = 0; i < b->groups; ++i) {
        if (i)
            printf("  //  ");
        print_group(b->group[i]);
    }
    printf("\n");
}

void recurse(board_t *b, int unused, group_t *g0, group_t *g1) {
    int k = b->k;
    board_t *nb = new_board(k + 1, unused, g0, g1);

    ++board_count;
    if (k <= feedback) {
        print_board(nb);
    }

    if (k > best_k) {
        best_k = k;
        unref_board(best_board);
        best_board = nb;
        ref_board(nb);
    }

    try_board(nb);
    unref_board(nb);
}

void try_board(board_t *b) {
    int k = b->k, unused = b->unused, groups = b->groups;
    group_t *g0 = b->group[0], *g1 = b->group[1];

    /* try making a new group */
    if (unused >= k) {
        int next_unused = unused - k;
        grouplist_t *gl = group_seed(k);
        for (int i = 0; i < gl->count; ++i) {
            if (groups) {
                recurse(b, next_unused, g0, gl->g[i]);
            } else {
                recurse(b, next_unused, gl->g[i], (group_t *)NULL);
            }
        }
        /* seed lists are persistent */
        /* free_grouplist(gl); */
    }

    for (int ig = 0; ig < groups; ++ig) {
        group_t *gi = b->group[ig];
        int *headsi = gi->sum_heads;
        int *chainsi = gi->sum_chains;
        int spare = unused;

        if (spare > 8)
            spare = 8;
        if (spare > k)
            spare = k;

        /* try extending an existing group */
        for (int diff = 0; diff <= spare; ++diff) {
            int rest = k - diff;

            if (rest > gi->maxsum)
                continue;
            for (int chain = headsi[rest]; chain >= 0; chain = chainsi[chain]) {
                int xi = chain / (gi->y + 2) - 1, yi = chain % (gi->y + 2) - 1;
                if (diff == 0) {
                    group_t *gn = group_place(gi, (loc_t){ xi, yi }, k);
                    if (ig) {
                        recurse(b, unused, g0, gn);
                    } else {
                        recurse(b, unused, gn, g1);
                    }
                } else {
                    int next_unused = unused - diff;
                    grouplist_t *gl = group_place_with(
                        gi, (loc_t){ xi, yi }, k, diff
                    );
                    for (int p = 0; p < gl->count; ++p) {
                        if (ig) {
                            recurse(b, next_unused, g0, gl->g[p]);
                        } else {
                            recurse(b, next_unused, gl->g[p], g1);
                        }
                    }
                    free_grouplist(gl);
                }
            }
        }

        /* try by coalesce */
        for (int jg = ig + 1; jg < groups; ++jg) {
            group_t *gj = b->group[jg];
            int *headsj = gj->sum_heads;
            int *chainsj = gj->sum_chains;

            for (int si = 1; si < k && si < gi->maxsum; ++si) {
                int need = k - si;
                int min = (need - unused < 1) ? 1 : (need - unused);

                for (int ci = headsi[si]; ci >= 0; ci = chainsi[ci]) {
                    loc_t li = (loc_t){
                        ci / (gi->y + 2) - 1, ci % (gi->y + 2) - 1
                    };
                    for (int sj = min; sj <= need && sj < gj->maxsum; ++sj) {
                        for (int cj = headsj[sj]; cj >= 0; cj = chainsj[cj]) {
                            loc_t lj = (loc_t){
                                cj / (gj->y + 2) - 1, cj % (gj->y + 2) - 1
                            }; 
                            int use = need - sj;
                            int next_unused = unused - use;
                            grouplist_t *gl = coalesce_group(
                                gi, li, gj, lj, k, use
                            );
                            for (int p = 0; p < gl->count; ++p) {
                                recurse(b, next_unused, gl->g[p], (group_t *)NULL);
                            }
                            free_grouplist(gl);
                        }
                    }
                }
            }
        }
    }
}
