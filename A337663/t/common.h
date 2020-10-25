#include <stdio.h>
#include <stdlib.h>

#include "../group.h"
#include "../sym.h"

#define true 1
#define false 0

extern void init_test(void);
extern void done_testing(void);

extern void ok(char *legend, ...);
extern void fail(char *legend, ...);
extern void fatal(char *legend, ...);
extern void is_bool(bool got, bool expect, char *legend, ...);
extern void is_int(int got, int expect, char *legend, ...);
extern void is_grid(int x, int y, int *got, int *expect, char *legend, ...);
extern void is_loc(loc_t got, loc_t expect, char *legend, ...);

extern void vok(char *legend, va_list argp);
extern void vfail(char *legend, va_list argp);
extern void vfatal(char *legend, va_list argp);
extern void vis_bool(bool got, bool expect, char *legend, va_list argp);
extern void vis_int(int got, int expect, char *legend, va_list argp);
extern void vis_grid(int x, int y, int *got, int *expect, char *legend, va_list argp);
extern void vis_loc(loc_t got, loc_t expect, char *legend, va_list argp);

extern int *parse_vals(int x, int y, char *str);
