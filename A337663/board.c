#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "board.h"
#include "sym.h"

typedef struct hist_s {
    int type;
    int index;
    int group_index;
    loc_t l1;
    loc_t l2;
} hist_t;

int n;
int freq;
int fcount;
int best_k;
board_t *best_board;
unsigned long board_count;
board_t *b0;

#define MAXHIST 100
#define MAXSTR 16
int in_histc = 0;
int next_hist = 0;
hist_t in_hist[MAXHIST];
hist_t out_hist[MAXHIST];
char out_histstr[MAXHIST * MAXSTR];

/*
 * Increment the refcount of a board.
 */
void ref_board(board_t *b) {
    ++b->refcount;
}

/*
 * Decrement the refcount of a board, freeing it if it reaches zero.
 */
void unref_board(board_t *b) {
    if (--b->refcount > 0)
        return;

    for (int i = 0; i < b->groups; ++i) {
        unref_group(b->group[i]);
    }
    free(b);
}

#define _d1(s) (*s++ - '0')
#define _d2(s) ({                               \
    int tens = _d1(s);                          \
    (tens < 0) ? -_d1(s) : tens * 10 + _d1(s);  \
})
#define _s(s) (void) ({ if (*s == ' ') ++s; })

/*
 * Perform initialization for a run to find A337663(n0), restarting
 * at the point indicated by start_hist if supplied, else from the start.
 */
board_t *init_board(int n0, int freq0, char *start_hist) {
    n = n0;
    freq = freq0;
    fcount = freq0;
    board_count = 0;

    b0 = new_board(2, n, (group_t *)NULL, (group_t *)NULL);
    best_k = 1;
    best_board = b0;
    ref_board(b0);

    in_histc = 0;
    if (start_hist) {
        char *s;
        best_k = strtol(start_hist, &s, 10);
        _s(s);
        in_histc = 2;
        next_hist = 2;
        while (isdigit(*s)) {
            in_hist[in_histc].type = _d1(s);
            switch (in_hist[in_histc].type) {
                case 1:
                    in_hist[in_histc].index = _d2(s);
                    break;
                case 2:
                    in_hist[in_histc].group_index = _d1(s);
                    in_hist[in_histc].l1.x = _d2(s);
                    in_hist[in_histc].l1.y = _d2(s);
                    break;
                case 3:
                    in_hist[in_histc].group_index = _d1(s);
                    in_hist[in_histc].index = _d2(s);
                    in_hist[in_histc].l1.x = _d2(s);
                    in_hist[in_histc].l1.y = _d2(s);
                    break;
                case 4:
                    in_hist[in_histc].group_index = _d1(s);
                    in_hist[in_histc].index = _d2(s);
                    in_hist[in_histc].l1.x = _d2(s);
                    in_hist[in_histc].l1.y = _d2(s);
                    in_hist[in_histc].l2.x = _d2(s);
                    in_hist[in_histc].l2.y = _d2(s);
                    break;
            }
            ++in_histc;
            _s(s);
        }
        printf("reset\n");
    }
    return b0;
}

/*
 * Clean up.
 */
void finish_board(void) {
    unref_board(best_board);
    unref_board(b0);
}

/*
 * Construct and return a new board structure with the specified details.
 * The refcount is initialised to 1.
 */
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

/*
 * Recursive coroutine with try_board(): construct a new board from this
 * board with new details, then call try_board() on it.
 *
 * Progress is tracked here. On return, all possible continuations of
 * the fresh board have been checked.
 */
void recurse(board_t *b, int unused, group_t *g0, group_t *g1) {
    int k = b->k;
    board_t *nb = new_board(k + 1, unused, g0, g1);
    out_histstr[MAXSTR * k] = 0;

    ++board_count;
    if (--fcount == 0) {
        fcount = freq;
        for (int i = k; i >= 2; --i) {
            hist_t h = out_hist[i];
            char *s = &out_histstr[MAXSTR * i];
            if (*s)
                break;
            switch (h.type) {
                case 1:
                    sprintf(s, "%01d%02d",
                            h.type, h.index);
                    break;
                case 2:
                    sprintf(s, "%01d%01d%02d%02d",
                            h.type, h.group_index, h.l1.x, h.l1.y);
                    break;
                case 3:
                    sprintf(s, "%01d%01d%02d%02d%02d",
                            h.type, h.group_index, h.index, h.l1.x, h.l1.y);
                    break;
                case 4:
                    sprintf(s, "%01d%01d%02d%02d%02d%02d%02d",
                            h.type, h.group_index, h.index,
                            h.l1.x, h.l1.y, h.l2.x, h.l2.y);
                    break;
            }
        }
        printf("%d ", best_k);
        for (int i = 2; i <= k; ++i)
            printf("%s ", &out_histstr[MAXSTR * i]);
        print_board(nb);
    }

    if (k >= best_k) {
        printf("%sbest ", (k == best_k) ? "e" : "");
        print_board(nb);
        if (k > best_k) {
            best_k = k;
            unref_board(best_board);
            best_board = nb;
            ref_board(nb);
        }
    }

    try_board(nb);
    unref_board(nb);
}

/*
 * Recursive coroutine with recurse(): try each way to extend this board.
 *
 * On return, best_k will contain the highest k seen, and best_board the
 * first board we saw with that k.
 */
void try_board(board_t *b) {
    int k = b->k, unused = b->unused, groups = b->groups;
    group_t *g0 = b->group[0], *g1 = b->group[1];
    hist_t h, *oh;

    if (k == next_hist) {
        ++next_hist;
        if (k < in_histc)
            h = in_hist[k];
        else {
            next_hist = 0;
            h.type = 0;
        }
    }
    else
        h.type = 0;
    oh = &out_hist[k];

    /* try making a new group */
    if (k <= 4 && unused >= k && h.type <= 1) {
        int next_unused = unused - k;
        grouplist_t *gl = group_seed(k);
        int start_index = (h.type == 1) ? h.index : 0;
        oh->type = 1;
        for (int i = start_index; i < gl->count; ++i) {
            oh->index = i;
            if (groups) {
                recurse(b, next_unused, g0, gl->g[i]);
            } else {
                recurse(b, next_unused, gl->g[i], (group_t *)NULL);
            }
        }
        /* seed lists are persistent */
        /* free_grouplist(gl); */
    }

    int start_group = (h.type > 1) ? h.group_index : 0;
    for (int ig = start_group; ig < groups; ++ig) {
        group_t *gi = b->group[ig];
        int *headsi = gi->sum_heads;
        int *chainsi = gi->sum_chains;
        int spare = unused;

        if (spare > 8)
            spare = 8;
        if (spare > k)
            spare = k;

        oh->group_index = ig;

        /* try extending an existing group */
        if (h.type < 4) {
            bool hist_wait = (h.type > 1 && h.type < 4) ? 1 : 0;
            for (int diff = 0; diff <= spare; ++diff) {
                int rest = k - diff;

                if (rest > gi->maxsum)
                    continue;
                for (int chain = headsi[rest]; chain >= 0; chain = chainsi[chain]) {
                    int xi = chain / (gi->y + 2) - 1, yi = chain % (gi->y + 2) - 1;
                    int p_start = 0;
                    if (hist_wait) {
                        if (xi != h.l1.x || yi != h.l1.y)
                            continue;
                        hist_wait = 0;
                        p_start = h.index;
                    }
                    if (diff == 0) {
                        group_t *gn = group_place(gi, (loc_t){ xi, yi }, k);
                        oh->type = 2;
                        oh->l1.x = xi;
                        oh->l1.y = yi;
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
                        oh->type = 3;
                        oh->l1.x = xi;
                        oh->l1.y = yi;
                        for (int p = p_start; p < gl->count; ++p) {
                            oh->index = p;
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
        }

        /* try by coalesce */
        for (int jg = ig + 1; jg < groups; ++jg) {
            group_t *gj = b->group[jg];
            int *headsj = gj->sum_heads;
            int *chainsj = gj->sum_chains;
            bool iwait = (h.type == 4) ? 1 : 0;

            oh->type = 4;
            for (int si = 1; si < k && si < gi->maxsum; ++si) {
                int need = k - si;
                int min = (need - unused < 1) ? 1 : (need - unused);

                for (int ci = headsi[si]; ci >= 0; ci = chainsi[ci]) {
                    loc_t li = (loc_t){
                        ci / (gi->y + 2) - 1, ci % (gi->y + 2) - 1
                    };
                    int jwait = 0;
                    if (iwait) {
                        if (li.x != h.l1.x || li.y != h.l1.y)
                            continue;
                        iwait = 0;
                        jwait = 1;
                    }
                    oh->l1.x = li.x;
                    oh->l1.y = li.y;
                    for (int sj = min; sj <= need && sj < gj->maxsum; ++sj) {
                        for (int cj = headsj[sj]; cj >= 0; cj = chainsj[cj]) {
                            loc_t lj = (loc_t){
                                cj / (gj->y + 2) - 1, cj % (gj->y + 2) - 1
                            }; 
                            int p_start = 0;
                            if (jwait) {
                                if (lj.x != h.l2.x || lj.y != h.l2.y)
                                    continue;
                                jwait = 0;
                                p_start = h.index;
                            }
                            int use = need - sj;
                            int next_unused = unused - use;
                            grouplist_t *gl = coalesce_group(
                                gi, li, gj, lj, k, use
                            );
                            oh->l2.x = lj.x;
                            oh->l2.y = lj.y;
                            for (int p = p_start; p < gl->count; ++p) {
                                oh->index = p;
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
