#ifndef LOC_H
#define LOC_H

typedef struct {
    int x;
    int y;
} loc_t;

static loc_t loc_diff(loc_t s1, loc_t s2) {
    loc_t d;
    d.x = s2.x - s1.x;
    d.y = s2.y - s1.y;
    return d;
}

static loc_t loc_rot90(loc_t src, loc_t diff) {
    loc_t d;
    d.x = src.x - diff.y;
    d.y = src.y + diff.x;
    return d;
}

static loc_t loc_rot270(loc_t src, loc_t diff) {
    loc_t d;
    d.x = src.x + diff.y;
    d.y = src.y - diff.x;
    return d;
}

static int loc_eq(loc_t s1, loc_t s2) {
    return s1.x == s2.x && s1.y == s2.y;
}

#endif
