#ifndef WALKER_H
#define WALKER_H

#include "pp.h"
#include "mbh.h"

typedef struct s_walk_result {
	union {
		struct s_walk_result* next;
		int next_discard;
	};
	int discard;
	int invsum;
	int nextbit;
	int vec[0];
} walk_result;

typedef struct s_walker {
	bhp heap;
	struct pp_s_pp* pp;
	int limit;
	int invsum;
	int vecsize;
	char* arena;
	walk_result* arenanext;
	int arenamax;
} walker;

#define MINARENA 16

extern void setup_walker(void);
extern void teardown_walker(void);
extern walker* new_walker(struct pp_s_pp* pp, int limit, int invsum);
extern void w_push_heap(walker* w, int sum, int nextbit);
extern walk_result* walker_next(walker* w);
extern walk_result* walker_find(walker* w);
extern void delete_walker(walker* w);
extern walk_result* wr_clone(walker* w, walk_result* wr);

#endif /* WALKER_H */
