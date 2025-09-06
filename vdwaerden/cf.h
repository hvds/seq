#ifndef CF_H
#define CF_H

#include <stdbool.h>
#include <stdint.h>

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long ullong;

static inline uint lsb32(uint32_t x) {
    /* assume x != 0 */
    return __builtin_ffs((int)x) - 1;
}

static inline uint lsb64(uint64_t x) {
    /* assume x != 0 */
    return __builtin_ffsll((long long)x) - 1;
}

#define UINT32_HIGHBIT (8 * sizeof(uint32_t) - 1)
static inline uint msb32(uint32_t x) {
    /* assume x != 0 */
    return UINT32_HIGHBIT - __builtin_clz((int)x);
}

#define UINT64_HIGHBIT (8 * sizeof(uint64_t) - 1)
/* FIXME: static assert that sizeof(long long) is 64 bits */
static inline uint msb64(uint64_t x) {
    /* assume x != 0 */
    return UINT64_HIGHBIT - __builtin_clzll((long long)x);
}

#endif
