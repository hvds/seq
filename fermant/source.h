#ifndef SOURCE_H
#define SOURCE_H

#include "frag.h"
#include "path.h"

extern void resolve_open(uint ri);
extern void resolve_close(uint ri);
extern void integrate_open(uint ri, uint pi);
extern void integrate_close(uint ri, uint pi);
extern bool read_frag(fid_t fi);
extern void write_frag(fid_t fi, pathset_t ps);

#endif /* SOURCE_H */
