#ifndef IPART_H
#define IPART_H

#include <stdio.h>

typedef struct ipart_s {
	unsigned int count;
	unsigned int parts[0];
} ipart_t;

void setup_ipart(void);
void teardown_ipart(void);

ipart_t* ipart_first(unsigned int major);
ipart_t* ipart_next(void);
void fprint_ipart(FILE* stream, ipart_t* ipart);

#endif
