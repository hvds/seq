#ifndef PART_H
#define PART_H

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/times.h>

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long long counter;

/* should normally be set as compile time define */
#ifndef NBASE
#define NBASE 4
#endif

#define NODES (1 << NBASE)

#endif
