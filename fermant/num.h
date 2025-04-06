#ifndef NUM_H
#define NUM_H

#include "types.h"

/* FIXME: cache this */
__inline int icomb(int x, int y) {
    if (y == 0)
        return 1;
    if (y < 0 || y > x)
        return 0;
    if (y == 1)
        return x;
    if (y + y > x)
        return icomb(x, x - y);
    return icomb(x - 1, y - 1) + icomb(x - 1, y);
}

__inline uint _ugcd(uint x, uint y) {
    return (x == 0) ? y : (x == 1) ? 1 : _ugcd(y % x, x);
}

__inline uint ugcd(uint x, uint y) {
    return (x > y) ? _ugcd(y, x) : _ugcd(x, y);
}

__inline int igcd(int x, int y) {
    int sign = (y < 0) ? -1 : 1;
    return sign * (int)ugcd(abs(x), abs(y));
}

#endif /* NUM_H */
