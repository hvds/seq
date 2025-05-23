#ifndef SOURCE_H
#define SOURCE_H

#include "frag.h"
#include "path.h"

extern bool resolve_open(uint ri);
extern void resolve_close(uint ri);
extern bool integrate_open(uint ri, uint pi);
extern void integrate_close(uint ri, uint pi);
extern bool read_frag(fid_t fi);
extern void write_frag(fid_t fi, pathset_t ps);
extern void resolve_checkpoint(void);
extern void integrate_checkpoint(void);
extern void set_resolve_cp(uint ri, off_t oi, off_t oor, off_t oop1, off_t oop2);
extern void set_integrate_cp(uint ri, uint pi, off_t oi);

extern void mmfrag_init(void);
extern void mmfrag_done(void);
extern uint mmfrag_open(bool resolved, uint ri, uint pi);
extern frag_t *mmfrag_get(uint rec);

#endif /* SOURCE_H */
