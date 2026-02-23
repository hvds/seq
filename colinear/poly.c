#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "poly.h"

uchar rev_lookup[256];
t_shape *transed;

void shape_free(t_shape *s) {
    free(s->points);
    free(s->disallowed);
    free(s);
}

t_shape *shape_new(t_dim d) {
    t_shape *s = malloc(sizeof(t_shape));
    size_t size = poly_size(d);
    s->points = malloc(size);
    s->disallowed = malloc(size);
    return s;
}

void init_poly(void) {
    bzero(rev_lookup, 256);
    for (uint i = 0; i < 8; ++i) {
        uint from = 1 << i;
        uint to = 1 << (7 - i);
        for (uint j = 0; j < 256; ++j) {
            if (j & from)
                rev_lookup[j] |= to;
        }
    }
    transed = shape_new((t_dim){ MAXDIM, MAXDIM });
}

void done_poly(void) {
    shape_free(transed);
}

void shape_reset(t_shape *s) {
    bzero(s->points, poly_size(s->d));
    bzero(s->disallowed, poly_size(s->d));
}

void shape_empty(t_shape *s) {
    s->d = (t_dim){ 1, 1 };
    shape_reset(s);
}

bool shiftcpy(uchar *dest, uchar *src, size_t n) {
    bool carry = 0;
    uint off = 0;
    while (n >= 8) {
        uint64_t v = *(uint64_t *)&src[off];
        *(uint64_t *)&dest[off] = (v << 1) | carry;
        carry = (v & 0x8000000000000000ULL) ? 1 : 0;
        off += 8;
        n -= 8;
    }
    if (n >= 4) {
        uint32_t v = *(uint32_t *)&src[off];
        *(uint32_t *)&dest[off] = (v << 1) | carry;
        carry = (v & 0x80000000UL) ? 1 : 0;
        off += 4;
        n -= 4;
    }
    if (n >= 2) {
        uint16_t v = *(uint16_t *)&src[off];
        *(uint16_t *)&dest[off] = (v << 1) | carry;
        carry = (v & 0x8000U) ? 1 : 0;
        off += 2;
        n -= 2;
    }
    if (n >= 1) {
        uint8_t v = *(uint8_t *)&src[off];
        *(uint8_t *)&dest[off] = (v << 1) | carry;
        carry = (v & 0x8000U) ? 1 : 0;
        off += 2;
        n -= 2;
    }
    return carry;
}

void expandcpy(uchar *dest, uchar *src, uint width, uint rows) {
    for (uint i = 0; i < rows; ++i) {
        memcpy(dest, src, width);
        dest[width] = 0;
        src += width;
        dest += width + 1;
    }
}

void explodecpy(uchar *dest, uchar *src, uint width, uint rows) {
    for (uint i = 0; i < rows; ++i) {
        dest[width] = shiftcpy(dest, src, width);
        src += width;
        dest += width + 1;
    }
}

void report_shape(char *buf, size_t bufsize, t_shape *s) {
    t_point p;
    uint off = 0;
    off += snprintf(&buf[off], bufsize - off,
            "%u:%u:", s->d.x, s->d.y);
    for (p.x = 0; p.x <= s->d.x + 1; ++p.x) {
        if (p.x)
            off += snprintf(&buf[off], bufsize - off, "/");
        for (p.y = 0; p.y <= s->d.y + 1; ++p.y) {
            if (shape_test_point(s, p))
                off += snprintf(&buf[off], bufsize - off, "*");
            else if (shape_test_disallowed(s, p))
                off += snprintf(&buf[off], bufsize - off, "o");
            else
                off += snprintf(&buf[off], bufsize - off, ".");
        }
    }
}

void diag_shape(t_shape *s, char *legend) {
    printf("%s [%u x %u]\n", legend, s->d.x, s->d.y);
    t_point p;
    for (p.x = 0; p.x <= s->d.x + 1; ++p.x) {
        printf("  ");
        for (p.y = 0; p.y <= s->d.y + 1; ++p.y) {
            if (shape_test_point(s, p))
                printf("*");
            else if (shape_test_disallowed(s, p))
                printf("o");
            else
                printf(".");
        }
        printf("\n");
    }
}

t_shape *shape_append(t_shape *dest, t_shape *src, t_point p) {
    bool shiftx = 0, shifty = 0;
    if (p.x == 0 || p.x > src->d.x)
        shiftx = 1;
    if (p.y == 0 || p.y > src->d.y)
        shifty = 1;
    dest->d = (t_dim){ src->d.x + shiftx, src->d.y + shifty };

    size_t psize = poly_size(src->d);
    size_t rsize = row_size(src->d);
    bool shiftr = (row_size(dest->d) != rsize);
    if (!shifty || (!shiftr && p.y != 0)) {
        if (!shiftx) {
            /* straight copy */
            memcpy(dest->points, src->points, psize);
            memcpy(dest->disallowed, src->disallowed, psize);
        } else if (p.x != 0) {
            /* straight copy, append extra row */
            memcpy(dest->points, src->points, psize);
            memcpy(dest->disallowed, src->disallowed, psize);
            bzero(&dest->points[psize], rsize);
            bzero(&dest->disallowed[psize], rsize);
        } else {
            /* straight copy, prepend extra row */
            bzero(dest->points, rsize);
            bzero(dest->disallowed, rsize);
            memcpy(&dest->points[rsize], src->points, psize);
            memcpy(&dest->disallowed[rsize], src->disallowed, psize);
        }
    } else if (!shiftr) {
        /* row shifts 1 bit, but same width */
        if (!shiftx) {
            /* just shift */
            shiftcpy(dest->points, src->points, psize);
            shiftcpy(dest->disallowed, src->disallowed, psize);
        } else if (p.x != 0) {
            /* shift, then append row */
            shiftcpy(dest->points, src->points, psize);
            shiftcpy(dest->disallowed, src->disallowed, psize);
            bzero(&dest->points[psize], rsize);
            bzero(&dest->disallowed[psize], rsize);
        } else {
            /* prepend row, then shift */
            bzero(dest->points, rsize);
            bzero(dest->disallowed, rsize);
            shiftcpy(&dest->points[rsize], src->points, psize);
            shiftcpy(&dest->disallowed[rsize], src->disallowed, psize);
        }
    } else if (p.y != 0) {
        /* rows are a byte longer, but no shift */
        if (!shiftx) {
            /* just expand row */
            expandcpy(dest->points, src->points, rsize, src->d.x + 2);
            expandcpy(dest->disallowed, src->disallowed, rsize, src->d.x + 2);
        } else if (p.x != 0) {
            /* expand, then append row */
            expandcpy(dest->points, src->points, rsize, src->d.x + 2);
            expandcpy(dest->disallowed, src->disallowed, rsize, src->d.x + 2);
            bzero(&dest->points[psize], rsize + 1);
            bzero(&dest->disallowed[psize], rsize + 1);
        } else {
            /* prepend row, then expand */
            bzero(dest->points, rsize + 1);
            bzero(dest->disallowed, rsize + 1); 
            expandcpy(&dest->points[rsize], src->points, rsize, src->d.x + 2);
            expandcpy(&dest->disallowed[rsize], src->disallowed, rsize, src->d.x + 2);
        }
    } else {
        /* rows are a byte longer, and shift */
        if (!shiftx) {
            /* just explode row */
            explodecpy(dest->points, src->points, rsize, src->d.x + 2);
            explodecpy(dest->disallowed, src->disallowed, rsize, src->d.x + 2);
        } else if (p.x != 0) {
            /* explode, then append row */
            explodecpy(dest->points, src->points, rsize, src->d.x + 2);
            explodecpy(dest->disallowed, src->disallowed, rsize, src->d.x + 2);
            bzero(&dest->points[psize], rsize + 1);
            bzero(&dest->disallowed[psize], rsize + 1);
        } else {
            /* prepend row, then explode */
            bzero(dest->points, rsize + 1);
            bzero(dest->disallowed, rsize + 1); 
            explodecpy(&dest->points[rsize], src->points, rsize, src->d.x + 2);
            explodecpy(&dest->disallowed[rsize], src->disallowed, rsize, src->d.x + 2);
        }
    }
    if (p.x == 0)
        p.x = 1;
    if (p.y == 0)
        p.y = 1;
    shape_mark_point(dest, p);
    shape_mark_disallowed(dest, p);
    return dest;
}

/* see also https://stackoverflow.com/questions/3534535/whats-a-time-efficient-algorithm-to-copy-unaligned-bit-arrays */
void bitmove(
    uchar *dptr, uint doff, uchar *sptr, uint soff, uint len
) {
    while (len) {
        uint need = 8 - (doff & 7);
        uint have = 8 - (soff & 7);
        if (have > len)
            have = len;
        if (need == 8)
            dptr[doff >> 3] = 0;
        uint try = (need < have) ? need : have;
        dptr[doff >> 3] |= (
            (sptr[soff >> 3] >> (soff & 7)) & ((1 << try) - 1)
        ) << (doff & 7);
        doff += try;
        soff += try;
        len -= try;
    }
}

void mem_shift_or(
    uchar *dest, uchar *src, size_t length, sint shift
) {
    if (shift == -1) {
        for (size_t i = 0; i < length - 1; ++i)
            dest[i] |= (src[i] >> 1) | ((src[i + 1] & 1) << 7);
        dest[length - 1] |= src[length - 1] >> 1;
    } else {
        dest[0] |= (src[0] & 0x7f) << 1;
        for (size_t i = 1; i < length; ++i)
            dest[i] |= ((src[i] & 0x7f) << 1) | (src[i - 1] >> 7);
    }
}

void mem_or(uchar *dest, uchar *src, size_t length) {
    while (length >= sizeof(uint)) {
        *(uint *)dest |= *(uint *)src;
        dest += sizeof(uint);
        src += sizeof(uint);
        length -= sizeof(uint);
    }
    while (length--)
        *dest++ |= *src++;
}

void mem_andn(uchar *dest, uchar *src, size_t length) {
    while (length >= sizeof(uint)) {
        *(uint *)dest &= ~*(uint *)src;
        dest += sizeof(uint);
        src += sizeof(uint);
        length -= sizeof(uint);
    }
    while (length--)
        *dest++ &= ~*src++;
}

/* initialize the provided iterator with the neighbouring points to the
 * provided shape. */
void shape_neighbours(t_iter *it, t_shape *s) {
    it->d = s->d;
    it->pos = (t_point){ 0, 1 };
    uchar *nb = &it->p[0];
    if (s->d.x == 1 && s->d.y == 1 && !shape_test_point(s, (t_point){ 1, 1 })) {
        bzero(nb, poly_size(it->d));
        poly_mark(it->d, nb, (t_point){ 1, 1 });
        return;
    }

    size_t rowsize = row_size(it->d);
    /* the memcpy() below means we only need clear the first two rows */
    bzero(nb, rowsize * 2);
    for (uint i = 1; i <= s->d.x; ++i) {
        uchar *src = &s->points[i * rowsize];
        uchar *dest = &nb[i * rowsize];
        memcpy(dest + rowsize, src, rowsize);   /* neighbour (+1, 0) */
        mem_or(dest - rowsize, src, rowsize);   /* neighbour (-1, 0) */
        mem_shift_or(dest, src, rowsize, +1);   /* neighbour (0, +1) */
        mem_shift_or(dest, src, rowsize, -1);   /* neighbour (0, -1) */
    }
    for (uint i = 0; i <= s->d.x + 1; ++i) {
        uchar *src = &s->disallowed[i * rowsize];
        uchar *dest = &nb[i * rowsize];
        mem_andn(dest, src, rowsize);
    }
}

void shape_write_pack2(t_shape *s, uchar *p, uchar *dptr) {
    uint size = row_size(s->d);
    uint doff = 0;
    for (uint i = 1; i <= s->d.x; ++i) {
        bitmove(dptr, doff, (uchar *)&p[i * size], 1, s->d.y);
        doff += s->d.y;
    }
}

void shape_poly_invert(t_shape *s, uchar *src, uchar *dest) {
    uint size = row_size(s->d);
    if (src == dest) {
        uchar temp[size];
        for (uint i = 0; i + i < s->d.x + 1; ++i) {
            uint j = s->d.x + 1 - i;
            memcpy(temp, src + i * size, size);
            memcpy(dest + i * size, src + j * size, size);
            memcpy(dest + j * size, temp, size);
        }
    } else {
        for (uint i = 0; i <= s->d.x + 1; ++i) {
            uint j = s->d.x + 1 - i;
            memcpy(dest + j * size, src + i * size, size);
        }
    }
}

void shape_poly_reverse(t_shape *s, uchar *src, uchar *dest) {
    uint size = row_size(s->d);
    uchar temp[size + 1];
    uint shift = 7 - ((s->d.y + 1) & 7);
    for (uint i = 0; i <= s->d.x + 1; ++i) {
        if (shift) {
            temp[0] = 0;
            bitmove(temp, shift, src + i * size, 0, s->d.y + 1);
        } else
            memmove(temp, src + i * size, size);
        for (uint j = 0; j < size; ++j) {
            uint k = size - 1 - j;
            dest[i * size + k] = rev_lookup[ temp[j] ];
        }
    }
}

void shape_poly_trans(t_shape *s, uchar *src, uchar *dest) {
    t_dim td = (t_dim){ s->d.y, s->d.x };
    bzero(dest, poly_size(td));
    for (uint i = 1; i <= s->d.x; ++i)
        for (uint j = 1; j <= s->d.y; ++j)
            if (shape_test_point(s, (t_point){ i, j }))
                poly_mark(td, dest, (t_point){ j, i });
}

size_t shape_write_canonical(t_shape *s, uchar *p2) {
    t_dim d = s->d;
    bool do_trans = (d.x > d.y) ? 1 : 0;
    size_t size = poly_size(d);
    if (do_trans) {
        /* poly size is asymmetric in (x, y), make sure we have room */
        size_t tsize = poly_size((t_dim){ d.y, d.x });
        if (tsize > size)
            size = tsize;
    }
    *p2++ = (uchar)(do_trans ? d.y : d.x);
    *p2++ = (uchar)(do_trans ? d.x : d.y);

    size_t p2size = pack2_size(d);
    uchar p2cur[p2size], cur[size];

    if (d.x <= d.y) {
        shape_write_pack2(s, s->points, p2);

        shape_poly_invert(s, s->points, cur);
        shape_write_pack2(s, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_reverse(s, s->points, cur);
        shape_write_pack2(s, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_invert(s, cur, cur);
        shape_write_pack2(s, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);
    }
    if (d.x >= d.y) {
        t_shape *t = transed;
        t->d = (t_dim){ d.y, d.x };
        shape_poly_trans(s, s->points, t->points); 

        if (d.x > d.y) {
            shape_write_pack2(t, t->points, p2);
        } else {
            shape_write_pack2(t, t->points, p2cur);
            if (memcmp(p2, p2cur, p2size) > 0)
                memcpy(p2, p2cur, p2size);
        }

        shape_poly_invert(t, t->points, cur);
        shape_write_pack2(t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_reverse(t, cur, cur);
        shape_write_pack2(t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_invert(t, cur, cur);
        shape_write_pack2(t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);
    }
    assert(bitcount(p2size, p2) == bitcount(poly_size(d), s->points));
    return 2 + p2size; 
}

t_point shape_iter(t_iter *i) {
    t_dim d = i->d;
    t_point pos = i->pos;
    uchar *p = &i->p[0];
    while (pos.x <= d.x + 1) {
        while (pos.y <= d.y + 1) {
            if (poly_test(d, p, pos)) {
                i->pos = (t_point){ pos.x, pos.y + 1 };
                return pos;
            }
            ++pos.y;
        }
        ++pos.x;
        pos.y = 0;
    }
    return (t_point){ 0, 0 };
}

/* Returns false if the suggested point would form any line of more than lim
 * points in the shape, else true.
 */
bool test_colinear(t_shape *s, t_point p, uint lim) {
    /* existing points are in the range {1, 1} .. {s->d.x, s->d.y} */
    uint exp = (p.x < s->d.x) ? s->d.x - p.x : 0;
    uint exm = (p.x > 1) ? p.x - 1 : 0;
    uint eyp = (p.y < s->d.y) ? s->d.y - p.y : 0;
    uint eym = (p.y > 1) ? p.y - 1 : 0;
    t_point q;
    uint seen;
    for (uint fx = 1; (exp / fx) + (exm / fx) >= lim; ++fx) {
        for (uint fy = 1; (eyp / fy) + (eym / fy) >= lim; ++fy) {
            seen = 0;
            for (uint m = 1; p.x + m * fx <= s->d.x && p.y + m * fy <= s->d.y; ++m) {
                q = (t_point){ p.x + m * fx, p.y + m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= lim)
                    return 0;
            }
            for (uint m = 1; p.x >= m * fx + 1 && p.y >= m * fy + 1; ++m) {
                q = (t_point){ p.x - m * fx, p.y - m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= lim)
                    return 0;
            }

            seen = 0;
            for (uint m = 1; p.x + m * fx <= s->d.x && p.y >= m * fy + 1; ++m) {
                q = (t_point){ p.x + m * fx, p.y - m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= lim)
                    return 0;
            }
            for (uint m = 1; p.x >= m * fx + 1 && p.y + m * fy <= s->d.y; ++m) {
                q = (t_point){ p.x - m * fx, p.y + m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= lim)
                    return 0;
            }
        }
    }
    if (p.x >= 1 && p.x <= s->d.x) {
        seen = 0;
        for (uint qy = 1; qy <= s->d.y; ++qy) {
            q = (t_point){ p.x, qy };
            if (qy == p.y || !shape_test_point(s, q))
                continue;
            ++seen;
            if (seen >= lim)
                return 0;
        }
    }
    if (p.y >= 1 && p.y <= s->d.y) {
        seen = 0;
        for (uint qx = 1; qx <= s->d.x; ++qx) {
            q = (t_point){ qx, p.y };
            if (qx == p.x || !shape_test_point(s, q))
                continue;
            ++seen;
            if (seen >= lim)
                return 0;
        }
    }

    return 1;
}

