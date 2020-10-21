#ifndef GROUP_H
#define GROUP_H

typedef enum {
    AVAIL = 0,
    RES = 1,
    USED = 2,
} avail_t;

typedef struct loc_s {
    int x;
    int y;
} loc_t;

typedef struct group_s {
    int x;
    int y;
    int sym;
    int maxsum;
    int refcount;
    int *vals;
    int **tvals;
    avail_t *avail;
    avail_t **tavail;
    int *sum_heads;
    int *sum_chains;
} group_t;

typedef struct grouplist_t {
    int count;
    group_t *g[0];
} grouplist_t;

extern void init_group(void);
extern void finish_group(void);
extern void print_group(group_t *g);
extern void ref_group(group_t *g);
extern void unref_group(group_t *g);

extern group_t *new_group(int x, int y, int sym, int* vals);
extern group_t *group_place(group_t *g, loc_t loc, int k);

extern void free_grouplist(grouplist_t *gl);
extern grouplist_t *group_seed(int k);
extern grouplist_t *group_place_with(group_t *g, loc_t loc, int k, int use);
extern grouplist_t *coalesce_group(
    group_t *g1, loc_t loc1, group_t *g2, loc_t loc2, int k, int use
);

#endif
