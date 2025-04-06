#ifndef CALC_H
#define CALC_H

#include <gmp.h>

#include "types.h"
#include "frag.h"

extern mpq_t totalq;

extern void init_calc(uint num_paths);
extern void done_calc(void);
extern void integrate(fid_t fi);
extern void report_total();

#endif /* CALC_H */
