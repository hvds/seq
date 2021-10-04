#ifndef SYM_H
#define SYM_H

#include "loc.h"

typedef unsigned char bool;

#define xy 0x01 /* null symmetry */
#define xY 0x02 /* reflect Y */
#define Xy 0x04 /* reflect X */
#define XY 0x08 /* rot-180 */
#define yx 0x10 /* reflect x=y */
#define yX 0x20 /* rot-90 */
#define Yx 0x40 /* rot-270 */
#define YX 0x80 /* reflect x=-y */

/* Bitwise-or of the above symmetry types */
typedef unsigned char sym_t;

extern sym_t sym_check(loclist_t *ll, span_t span, int size);
extern bool sym_best(sym_t s, span_t span, loc_t p);
extern bool sym_best2(sym_t s, span_t span, loc_t p1, loc_t p2);
extern bool sym_axis(sym_t s, span_t span, loc_t p1, loc_t p2);

#endif
