#ifndef DIAG_H
#define DIAG_H 1

extern int diag_size;
extern int diag(char *format, ...);
extern void diag_reset(void);
extern void keep_diag(void);
extern int init_diag(void);

#endif
