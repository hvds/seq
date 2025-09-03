#ifndef CF_H
#define CF_H

#include <stdbool.h>

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long ullong;

static inline uint lsb32(uint x) {
    /* assume x != 0 */
    return __builtin_ffs((int)x) - 1;
}

static inline uint lsb64(ullong x) {
    /* assume x != 0 */
    return __builtin_ffsll((long long)x) - 1;
}

#define UINT_HIGHBIT (8 * sizeof(uint) - 1)
static inline uint msb32(uint x) {
    /* assume x != 0 */
    return UINT_HIGHBIT - __builtin_clz((int)x);
}

#define ULLONG_HIGHBIT (8 * sizeof(ullong) - 1)
static inline uint msb64(ullong x) {
    /* assume x != 0 */
    return ULLONG_HIGHBIT - __builtin_clzll((long long)x);
}

#endif
