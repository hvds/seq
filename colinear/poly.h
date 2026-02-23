#ifndef POLY_H
#define POLY_H

typedef unsigned char uchar;
typedef unsigned int uint;
typedef signed int sint;
typedef unsigned long ulong;
typedef unsigned long long ullong;
typedef unsigned char bool;

/* we use a single byte to store width and height */
#define MAXDIM 255

typedef struct s_point {
    uint x;
    uint y;
} t_point, t_dim;

typedef struct s_shape {
    t_dim d;
    size_t index;
    uchar *points;
    uchar *disallowed;
} t_shape;

typedef struct s_iter {
    t_dim d;
    t_point pos;
    uchar p[0];
} t_iter;

static inline size_t row_size(t_dim d) {
    uint bits = d.y + 2;
    return (bits + 7) / 8;
}

static inline size_t poly_size(t_dim d) {
    return (d.x + 2) * row_size(d);
}

static inline size_t max_poly_size(void) {
    return poly_size((t_dim){ MAXDIM, MAXDIM });
}

static inline size_t pack2_size(t_dim d) {
    return (d.x * d.y + 7) / 8;
}

static inline size_t report_size(t_dim d) {
    /* "X:Y:<bits>\0" */
#if MAXDIM < 1000
    return 3 + 1 + 3 + (d.x + 2) * (d.y + 3) + 1;
#else
#   error "check MAXDIM"
#endif
}

static inline bool poly_test(t_dim d, uchar *py, t_point p) {
    uint bit = p.x * row_size(d) * 8 + p.y;
    return (
        py[bit / 8] & (1 << (bit % 8))
    ) ? 1 : 0;
}

static inline void poly_unmark(t_dim d, uchar *py, t_point p) {
    uint bit = p.x * row_size(d) * 8 + p.y;
    py[bit / 8] &= ~(1 << (bit % 8));
}

static inline void poly_mark(t_dim d, uchar *py, t_point p) {
    uint bit = p.x * row_size(d) * 8 + p.y;
    py[bit / 8] |= 1 << (bit % 8);
}

static inline bool shape_test_point(t_shape *s, t_point p) {
    return poly_test(s->d, s->points, p);
}

static inline bool shape_test_disallowed(t_shape *s, t_point p) {
    return poly_test(s->d, s->disallowed, p);
}

static inline void shape_mark_point(t_shape *s, t_point p) {
    poly_mark(s->d, s->points, p);
}

static inline void shape_mark_disallowed(t_shape *s, t_point p) {
    poly_mark(s->d, s->disallowed, p);
}

static inline void iter_remove(t_iter *i, t_point p) {
    poly_unmark(i->d, &i->p[0], p);
}

static inline void iter_reset(t_iter *i) {
    i->pos = (t_point){ 0, 1 };
}

static inline uint bitcount(size_t size, uchar *poly) {
    uint count = 0;
    for (uint i = 0; i < size; ++i)
        count += __builtin_popcount(poly[i]);
    return count;
}

extern void init_poly(void);
extern void done_poly(void);
extern void shape_free(t_shape *s);
extern t_shape *shape_new(t_dim d);
extern void shape_empty(t_shape *s);
extern void shape_reset(t_shape *s);
extern void report_shape(char *buf, size_t bufsize, t_shape *s);
extern void diag_shape(t_shape *s, char *legend);
extern t_shape *shape_append(t_shape *dest, t_shape *src, t_point p);
extern void shape_neighbours(t_iter *it, t_shape *s);
extern size_t shape_write_canonical(t_shape *s, uchar *p2);
extern t_point shape_iter(t_iter *i);
extern bool test_colinear(t_shape *s, t_point p, uint lim);

#endif
