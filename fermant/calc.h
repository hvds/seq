#ifndef CALC_H
#define CALC_H

#include <gmp.h>

#include "types.h"
#include "frag.h"

extern mpq_t totalq;
extern mpq_t *path_total;   /* path_total[npaths] */

extern void init_calc(uint num_paths);
extern void done_calc(void);
extern void integrate_path(uint pi);
extern void report_total();

#endif /* CALC_H */
