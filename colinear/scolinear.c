#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long ullong;
typedef unsigned char bool;

uint opt_n;
bool opt_check = 0;
bool opt_recover = 0;
char path_buf[4096];
/* we probably want a size around 1Gb, but will test with 10kb */
/* #define MAX_SORT 10240 */
#define MAX_SORT 1073741824

uchar rev_lookup[256];
struct stat statbuf;

typedef enum {
    /* size, pack3 points+disallowed */
    REC_A = 0xf3,
    /* index + REC_A */
    REC_B,
    /* (partial) canonical + REC_B */
    REC_C,
    /* alignment marker */
    REC_ALIGN
} e_rec;

typedef struct s_reader {
    uint fake;
    e_rec type;
    FILE *f;
} t_reader;
typedef struct s_writer {
    e_rec type;
    FILE *f;
} t_writer;

typedef uchar t_poly;
typedef struct s_point {
    uint x;
    uint y;
} t_point;
typedef struct s_shape {
    uint x;
    uint y;
    size_t index;
    t_poly *points;
    t_poly *disallowed;
} t_shape;

typedef struct s_iter {
    uint x;
    uint y;
    uint i;
    uint j;
    t_poly p[0];
} t_iter;

typedef struct s_buf {
    size_t bsize;
    size_t csize;
    size_t bcount;
    size_t ccount;
    uchar *buf;
    ssize_t *ptr;
} t_buf;

double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000;
}

double seconds(double t1) {
    return (t1 - t0);
}

double elapsed(void) {
    return seconds(utime());
}

char *path_base(uint k) {
    snprintf(path_buf, sizeof(path_buf), "log/c%u.%u", opt_n, k);
    return path_buf;
}

char *path_unsorted(uint k) {
    snprintf(path_buf, sizeof(path_buf), "log/c%u.%uu", opt_n, k);
    return path_buf;
}

void reader_close(t_reader *r) {
    if (r->fake == 0)
        fclose(r->f);
    free(r);
}

void writer_close(t_writer *w) {
    fclose(w->f);
    free(w);
}

size_t reader_size(t_reader *r) {
    if (fstat(fileno(r->f), &statbuf)) {
        fprintf(stderr, "read size %02x: %s\n",
                r->type, strerror(errno));
        exit(1);
    }
    return (size_t)statbuf.st_size;
}

t_reader *reader_previous(uint k) {
    t_reader *r = malloc(sizeof(t_reader));
    if (k == 1) {
        r->fake = 1;
    } else {
        char *path = path_base(k - 1);
        r->f = fopen(path, "rb");
        if (!r->f) {
            fprintf(stderr, "open r %s: %s\n", path, strerror(errno));
            exit(1);
        }
        r->fake = 0;
    }
    r->type = REC_A;
    return r;
}

t_reader *reader_cur(uint k) {
    t_reader *r = malloc(sizeof(t_reader));
    char *path = path_unsorted(k);
    r->f = fopen(path, "rb");
    if (!r->f) {
        fprintf(stderr, "open r %s: %s\n", path, strerror(errno));
        exit(1);
    }
    r->fake = 0;
    r->type = REC_A;
    return r;
}

t_writer *writer_new(uint k, e_rec type) {
    t_writer *w = malloc(sizeof(t_writer));
    char *path = path_unsorted(k);
    w->f = fopen(path, "wb");
    if (!w->f) {
        fprintf(stderr, "open w %s: %s\n", path, strerror(errno));
        exit(1);
    }
    w->type = type;
    return w;
}

t_writer *writer_uniq(uint k) {
    t_writer *w = malloc(sizeof(t_writer));
    char *path = path_base(k);
    w->f = fopen(path, "wb");
    if (!w->f) {
        fprintf(stderr, "open w %s: %s\n", path, strerror(errno));
        exit(1);
    }
    w->type = REC_A;
    return w;
}

void shape_free(t_shape *s) {
    free(s->points);
    free(s->disallowed);
    free(s);
}

size_t raw_row_size(uint x, uint y) {
    uint bits = y + 2;
    return (bits + 7) / 8;
}

size_t shape_row_size(t_shape *s) {
    return raw_row_size(s->x, s->y);
}

size_t raw_poly_size(uint x, uint y) {
    return (x + 2) * raw_row_size(x, y);
}

size_t shape_poly_size(t_shape *s) {
    return raw_poly_size(s->x, s->y);
}

size_t raw_pack3_size(uint x, uint y) {
    uint bits = (x + 2) * (y + 2);
    return (bits + 4) / 5;
}

size_t shape_pack3_size(t_shape *s) {
    return raw_pack3_size(s->x, s->y);
}

size_t raw_pack2_size(uint x, uint y) {
    return (x * y + 7) / 8;
}

size_t shape_pack2_size(t_shape *s) {
    return raw_pack2_size(s->x, s->y);
}

size_t raw_rec_size(uint x, uint y, e_rec type) {
    uint size = 1;  /* type marker */
    switch (type) {
        default:
            fprintf(stderr, "panic, unknown type %02x\n", type);
            exit(1);
        case REC_C:
            size += 2;                      /* canonical size */
            /* pack2 size is symmetric in (x, y) */
            size += raw_pack2_size(x, y);   /* canonical data */
        case REC_B:
            size += 5;                      /* index */
        case REC_A:
            size += 2;                      /* oriented size */
            size += raw_pack3_size(x, y);   /* full data */
        case REC_ALIGN:
            ;                               /* align is type marker only */
    }
    return size;
}

size_t shape_rec_size(t_shape *s, e_rec type) {
    return raw_rec_size(s->x, s->y, type);
}

void shape_init(t_shape *s) {
    uint size = shape_poly_size(s);
    s->points = calloc(size, sizeof(t_poly));
    s->disallowed = calloc(size, sizeof(t_poly));
}

t_shape *shape_empty(void) {
    t_shape *s = malloc(sizeof(t_shape));
    s->x = 1;
    s->y = 1;
    shape_init(s);
    return s;
}

bool shape_test(t_shape *s, t_poly *py, t_point p) {
    uint bit = p.x * shape_row_size(s) * 8 + p.y;
    return (
        py[bit / (8 * sizeof(t_poly))] & (1 << (bit % (8 * sizeof(t_poly))))
    ) ? 1 : 0;
}

bool shape_test_point(t_shape *s, t_point p) {
    return shape_test(s, s->points, p);
}

bool shape_test_disallowed(t_shape *s, t_point p) {
    return shape_test(s, s->disallowed, p);
}

void shape_unmark(t_shape *s, t_poly *py, t_point p) {
    uint bit = p.x * shape_row_size(s) * 8 + p.y;
    py[bit / (8 * sizeof(t_poly))] &= ~(1 << (bit % (8 * sizeof(t_poly))));
}

void shape_mark(t_shape *s, t_poly *py, t_point p) {
    uint bit = p.x * shape_row_size(s) * 8 + p.y;
    py[bit / (8 * sizeof(t_poly))] |= 1 << (bit % (8 * sizeof(t_poly)));
}

void shape_mark_point(t_shape *s, t_point p) {
    shape_mark(s, s->points, p);
}

void shape_mark_disallowed(t_shape *s, t_point p) {
    shape_mark(s, s->disallowed, p);
}

t_iter *shape_neighbours(t_shape *s) {
    uint size = shape_poly_size(s);
    t_iter *it = calloc(1, sizeof(t_iter) + size);
    it->x = s->x;
    it->y = s->y;
    it->j = 1;
    t_poly *nb = &it->p[0];
    bool any = 0;
    t_point p, q;
    for (p.x = 1; p.x <= s->x; ++p.x) {
        for (p.y = 1; p.y <= s->y; ++p.y) {
            if (!shape_test_point(s, p))
                continue;
            any = 1;
            q = (t_point){ p.x + 1, p.y };
            if (!shape_test(s, nb, q) && !shape_test_disallowed(s, q))
                shape_mark(s, nb, q);
            q = (t_point){ p.x - 1, p.y };
            if (!shape_test(s, nb, q) && !shape_test_disallowed(s, q))
                shape_mark(s, nb, q);
            q = (t_point){ p.x, p.y + 1 };
            if (!shape_test(s, nb, q) && !shape_test_disallowed(s, q))
                shape_mark(s, nb, q);
            q = (t_point){ p.x, p.y - 1 };
            if (!shape_test(s, nb, q) && !shape_test_disallowed(s, q))
                shape_mark(s, nb, q);
        }
    }
    if (!any)
        shape_mark(s, nb, (t_point){ 1, 1 });
    return it;
}

void read_error(t_reader *r, char *legend) {
    fprintf(stderr, "read %02x %s: %s\n",
                r->type, legend, ferror(r->f) ? strerror(errno) : "got eof");
    exit(1);
}

bool read_head(t_reader *r) {
    uchar in[1];
    while (1) {
        size_t read = fread(in, 1, 1, r->f);
        if (read != 1) {
            if (!ferror(r->f))
                return 0;   /* eof */
            read_error(r, "head");
        }
        if (in[0] == (uchar)r->type)
            return 1;
        if (in[0] == (uchar)REC_ALIGN)
            continue;
        fprintf(stderr, "read head expected %02x got %02x\n",
                r->type, in[0]);
        exit(1);
    }
}

void read_xy(t_reader *r, t_shape *s) {
    uchar in[2];
    size_t read = fread(in, 1, 2, r->f);
    if (read != 2)
        read_error(r, "xy");
    s->x = in[0];
    s->y = in[1];
}

void read_pack3(t_reader *r, t_shape *s) {
    uint size = shape_pack3_size(s);
    uchar in[size];
    if (fread(in, 1, size, r->f) != size)
        read_error(r, "pack3");
    shape_init(s);
    t_point p = { 0, 0 };
    uint off = 0, rot = 0;
    while (1) {
        uint v = in[off] % 3;
        switch (v) {
          case 1:
            shape_mark_point(s, p);
            /* fall through */
          case 2:
            shape_mark_disallowed(s, p);
        }
        if (++rot >= 5) {
            ++off;
            rot = 0;
        } else
            in[off] /= 3;
        if (++p.y > s->y + 1) {
            p.y = 0;
            if (++p.x > s->x + 1)
                break;
        }
    }
}

t_shape *shape_fetch(t_reader *r) {
    if (r->fake) {
        if (r->fake == 1) {
            r->fake = 2;
            return shape_empty();
        }
        return NULL;
    }
    if (!read_head(r))
        return NULL;

    t_shape *s = malloc(sizeof(t_shape));
    uchar in[2];
    switch (r->type) {
      case REC_A:
        read_xy(r, s);
        read_pack3(r, s);
        break;
      default:
        fprintf(stderr, "shape_fetch type todo\n");
        exit(1);
    }
    return s;
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

void diag_shape(t_shape *s, char *legend) {
    fprintf(stderr, "%s [%u x %u]\n", legend, s->x, s->y);
    t_point p;
    for (p.x = 0; p.x <= s->x + 1; ++p.x) {
        fprintf(stderr, "  ");
        for (p.y = 0; p.y <= s->y + 1; ++p.y) {
            if (shape_test_point(s, p))
                fprintf(stderr, "*");
            else if (shape_test_disallowed(s, p))
                fprintf(stderr, "o");
            else
                fprintf(stderr, ".");
        }
        fprintf(stderr, "\n");
    }
}

t_shape *shape_append(t_shape *s0, t_point p) {
    bool shiftx = 0, shifty = 0;
    if (p.x == 0 || p.x > s0->x)
        shiftx = 1;
    if (p.y == 0 || p.y > s0->y)
        shifty = 1;
    t_shape *s = malloc(sizeof(t_shape));
    s->x = s0->x + shiftx;
    s->y = s0->y + shifty;
    shape_init(s);
    uint psize = shape_poly_size(s0);
    uint rsize = shape_row_size(s0);
    bool shiftr = (shape_row_size(s) != rsize);
    if (!shifty || (!shiftr && p.y != 0)) {
        if (!shiftx) {
            /* straight copy */
            memcpy(s->points, s0->points, psize);
            memcpy(s->disallowed, s0->disallowed, psize);
        } else if (p.x != 0) {
            /* straight copy, append extra row */
            memcpy(s->points, s0->points, psize);
            memcpy(s->disallowed, s0->disallowed, psize);
            bzero(&s->points[psize], rsize);
            bzero(&s->disallowed[psize], rsize);
        } else {
            /* straight copy, prepend extra row */
            bzero(s->points, rsize);
            bzero(s->disallowed, rsize);
            memcpy(&s->points[rsize], s0->points, psize);
            memcpy(&s->disallowed[rsize], s0->disallowed, psize);
        }
    } else if (!shiftr) {
        /* row shifts 1 bit, but same width */
        if (!shiftx) {
            /* just shift */
            shiftcpy(s->points, s0->points, psize);
            shiftcpy(s->disallowed, s0->disallowed, psize);
        } else if (p.x != 0) {
            /* shift, then append row */
            shiftcpy(s->points, s0->points, psize);
            shiftcpy(s->disallowed, s0->disallowed, psize);
            bzero(&s->points[psize], rsize);
            bzero(&s->disallowed[psize], rsize);
        } else {
            /* prepend row, then shift */
            bzero(s->points, rsize);
            bzero(s->disallowed, rsize);
            shiftcpy(&s->points[rsize], s0->points, psize);
            shiftcpy(&s->disallowed[rsize], s0->disallowed, psize);
        }
    } else if (p.y != 0) {
        /* rows are a byte longer, but no shift */
        if (!shiftx) {
            /* just expand row */
            expandcpy(s->points, s0->points, rsize, s0->x + 2);
            expandcpy(s->disallowed, s0->disallowed, rsize, s0->x + 2);
        } else if (p.x != 0) {
            /* expand, then append row */
            expandcpy(s->points, s0->points, rsize, s0->x + 2);
            expandcpy(s->disallowed, s0->disallowed, rsize, s0->x + 2);
            bzero(&s->points[psize], rsize + 1);
            bzero(&s->disallowed[psize], rsize + 1);
        } else {
            /* prepend row, then expand */
            bzero(s->points, rsize + 1);
            bzero(s->disallowed, rsize + 1); 
            expandcpy(&s->points[rsize], s0->points, rsize, s0->x + 2);
            expandcpy(&s->disallowed[rsize], s0->disallowed, rsize, s0->x + 2);
        }
    } else {
        /* rows are a byte longer, and shift */
        if (!shiftx) {
            /* just explode row */
            explodecpy(s->points, s0->points, rsize, s0->x + 2);
            explodecpy(s->disallowed, s0->disallowed, rsize, s0->x + 2);
        } else if (p.x != 0) {
            /* explode, then append row */
            explodecpy(s->points, s0->points, rsize, s0->x + 2);
            explodecpy(s->disallowed, s0->disallowed, rsize, s0->x + 2);
            bzero(&s->points[psize], rsize + 1);
            bzero(&s->disallowed[psize], rsize + 1);
        } else {
            /* prepend row, then explode */
            bzero(s->points, rsize + 1);
            bzero(s->disallowed, rsize + 1); 
            explodecpy(&s->points[rsize], s0->points, rsize, s0->x + 2);
            explodecpy(&s->disallowed[rsize], s0->disallowed, rsize, s0->x + 2);
        }
    }
    if (p.x == 0)
        p.x = 1;
    if (p.y == 0)
        p.y = 1;
    shape_mark_point(s, p);
    shape_mark_disallowed(s, p);
    return s;
}

void write_xy(t_writer *w, t_shape *s) {
    uchar out[2];
    out[0] = (uchar)s->x;
    out[1] = (uchar)s->y;
    if (fwrite(out, 1, 2, w->f) != 2) {
        fprintf(stderr, "write_xy: %s\n", strerror(errno));
        exit(1);
    }
}

void write_pack3(t_writer *w, t_shape *s) {
    uint size = shape_pack3_size(s);
    uchar out[size];
    t_point p = { 0, 0 };
    uint off = 0, rot = 0, mul = 1;
    while (1) {
        uint v = shape_test_point(s, p) ? 1
            : shape_test_disallowed(s, p) ? 2
            : 0;
        if (rot == 0)
            out[off] = v;
        else
            out[off] += v * mul;
        if (++rot >= 5) {
            ++off;
            rot = 0;
            mul = 1;
        } else
            mul *= 3;
        if (++p.y > s->y + 1) {
            p.y = 0;
            if (++p.x > s->x + 1)
                break;
        }
    }
    if (fwrite(out, 1, size, w->f) != size) {
        fprintf(stderr, "write pack3 %02x: %s\n", w->type, strerror(errno));
        exit(1);
    }
}

size_t shape_write_index(t_shape *s, uchar *p) {
    for (uint i = 0; i < 5; ++i)
        p[i] = (uchar)((s->index >> (8 * (5 - i - 1))) & 0xff);
    return 5;
}

size_t shape_write_data(t_shape *s, uchar *p) {
    size_t off = 0;
    p[off++] = s->x;
    p[off++] = s->y;
    t_point pt = { 0, 0 };
    uint rot = 0, mul = 1;
    while (1) {
        uint v = shape_test_point(s, pt) ? 1
            : shape_test_disallowed(s, pt) ? 2
            : 0;
        if (rot == 0)
            p[off] = v;
        else
            p[off] += v * mul;
        if (++rot >= 5) {
            ++off;
            rot = 0;
            mul = 1;
        } else
            mul *= 3;
        if (++pt.y > s->y + 1) {
            pt.y = 0;
            if (++pt.x > s->x + 1)
                break;
        }
    }
    return (rot) ? off + 1 : off;
}

void init_rev(void) {
    bzero(rev_lookup, 256);
    for (uint i = 0; i < 8; ++i) {
        uint from = 1 << i;
        uint to = 1 << (7 - i);
        for (uint j = 0; j < 256; ++j) {
            if (j & from)
                rev_lookup[j] |= to;
        }
    }
}

/* see also https://stackoverflow.com/questions/3534535/whats-a-time-efficient-algorithm-to-copy-unaligned-bit-arrays */
void bitmove(uchar *dptr, uint doff, uchar *sptr, uint soff, uint len) {
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

void shape_write_pack2(t_shape *s, t_poly *p, uchar *dptr) {
    uint size = shape_row_size(s);
    uint doff = 0;
    for (uint i = 1; i <= s->x; ++i) {
        bitmove(dptr, doff, (uchar *)&p[i * size], 1, s->y);
        doff += s->y;
    }
}

void shape_poly_invert(t_shape *s, t_poly *src, t_poly *dest) {
    uint size = shape_row_size(s);
    if (src == dest) {
        uchar temp[size];
        for (uint i = 0; i + i < s->x + 1; ++i) {
            uint j = s->x + 1 - i;
            memcpy(temp, src + i * size, size);
            memcpy(dest + i * size, src + j * size, size);
            memcpy(dest + j * size, temp, size);
        }
    } else {
        for (uint i = 0; i <= s->x + 1; ++i) {
            uint j = s->x + 1 - i;
            memcpy(dest + j * size, src + i * size, size);
        }
    }
}

void shape_poly_reverse(t_shape *s, t_poly *src, t_poly *dest) {
    uint size = shape_row_size(s);
    uchar temp[size + 1];
    uint shift = 7 - ((s->y + 1) & 7);
    for (uint i = 0; i <= s->x + 1; ++i) {
        if (shift) {
            temp[0] = 0;
            bitmove(temp, shift, src + i * size, 0, s->y + 1);
        } else
            memmove(temp, src + i * size, size);
        for (uint j = 0; j < size; ++j) {
            uint k = size - 1 - j;
            dest[i * size + k] = rev_lookup[ temp[j] ];
        }
    }
}

void shape_poly_trans(t_shape *s, t_poly *src, t_poly *dest) {
    t_shape t;
    t.x = s->y;
    t.y = s->x;
    bzero(dest, shape_poly_size(&t));
    for (uint i = 1; i <= s->x; ++i)
        for (uint j = 1; j <= s->y; ++j)
            if (shape_test_point(s, (t_point){ i, j }))
                shape_mark(&t, dest, (t_point){ j, i });
}

uint bitcount(size_t size, uchar *poly) {
    uint count = 0;
    for (uint i = 0; i < size; ++i)
        count += __builtin_popcount(poly[i]);
    return count;
}

size_t shape_write_canonical(t_shape *s, uchar *p2) {
    uint x = s->x, y = s->y;
    bool do_trans = (x > y) ? 1 : 0;
    size_t size = raw_poly_size(x, y);
    if (do_trans) {
        /* poly size is asymmetric in (x, y), make sure we have room */
        size_t tsize = raw_poly_size(y, x);
        if (tsize > size)
            size = tsize;
    }
    *p2++ = (uchar)(do_trans ? y : x);
    *p2++ = (uchar)(do_trans ? x : y);

    size_t p2size = shape_pack2_size(s);
    uchar p2cur[p2size], cur[size];

    if (x <= y) {
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
    if (x >= y) {
        t_shape t;
        t.x = y;
        t.y = x;
        shape_init(&t);

        shape_poly_trans(s, s->points, t.points);
        if (s->x > s->y) {
            shape_write_pack2(&t, t.points, p2);
        } else {
            shape_write_pack2(&t, t.points, p2cur);
            if (memcmp(p2, p2cur, p2size) > 0)
                memcpy(p2, p2cur, p2size);
        }

        shape_poly_invert(&t, t.points, cur);
        shape_write_pack2(&t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_reverse(&t, cur, cur);
        shape_write_pack2(&t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);

        shape_poly_invert(&t, cur, cur);
        shape_write_pack2(&t, cur, p2cur);
        if (memcmp(p2, p2cur, p2size) > 0)
            memcpy(p2, p2cur, p2size);
    }
    assert(bitcount(p2size, p2) == bitcount(shape_poly_size(s), s->points));
    return 2 + p2size;
}

void shape_write(t_shape *s, t_writer *w) {
    e_rec type = w->type;
    uchar out[1];
    out[0] = (uchar)type;
    if (fwrite(out, 1, 1, w->f) != 1) {
        fprintf(stderr, "write type %02x: %s\n", type, strerror(errno));
        exit(1);
    }
    switch (type) {
      case REC_A:
        write_xy(w, s);
        write_pack3(w, s);
        break;
      default:
        fprintf(stderr, "shape_write type todo\n");
        exit(1);
    }
}

void iter_free(t_iter *i) {
    free(i);
}

t_point shape_iter(t_iter *i) {
    t_shape t;
    t.x = i->x;
    t.y = i->y;
    t_poly *p = &i->p[0];
    while (i->i <= i->x + 1) {
        while (i->j <= i->y + 1) {
            if (shape_test(&t, p, (t_point){ i->i, i->j }))
                return (t_point){ i->i, i->j++ };
            ++i->j;
        }
        ++i->i;
        i->j = 0;
    }
    return (t_point){ 0, 0 };
}

void iter_remove(t_iter *i, t_point p) {
    t_shape t;
    t.x = i->x;
    t.y = i->y;
    shape_unmark(&t, &i->p[0], p);
}

void iter_reset(t_iter *i) {
    i->i = 0;
    i->j = 1;
}

/* Returns false if the suggested point would form any line of more than opt_n
 * points in the shape, else true.
 */
bool test_colinear(t_shape *s, t_point p) {
    /* existing points are in the range {1, 1} .. {s->x, s->y} */
    uint exp = (p.x < s->x) ? s->x - p.x : 0;
    uint exm = (p.x > 1) ? p.x - 1 : 0;
    uint eyp = (p.y < s->y) ? s->y - p.y : 0;
    uint eym = (p.y > 1) ? p.y - 1 : 0;
    t_point q;
    uint seen;
    for (uint fx = 1; (exp / fx) + (exm / fx) >= opt_n; ++fx) {
        for (uint fy = 1; (eyp / fy) + (eym / fy) >= opt_n; ++fy) {
            seen = 0;
            for (uint m = 1; p.x + m * fx <= s->x && p.y + m * fy <= s->y; ++m) {
                q = (t_point){ p.x + m * fx, p.y + m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= opt_n)
                    return 0;
            }
            for (uint m = 1; p.x >= m * fx + 1 && p.y >= m * fy + 1; ++m) {
                q = (t_point){ p.x - m * fx, p.y - m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= opt_n)
                    return 0;
            }

            seen = 0;
            for (uint m = 1; p.x + m * fx <= s->x && p.y >= m * fy + 1; ++m) {
                q = (t_point){ p.x + m * fx, p.y - m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= opt_n)
                    return 0;
            }
            for (uint m = 1; p.x >= m * fx + 1 && p.y + m * fy <= s->y; ++m) {
                q = (t_point){ p.x - m * fx, p.y + m * fy };
                if (!shape_test_point(s, q))
                    continue;
                ++seen;
                if (seen >= opt_n)
                    return 0;
            }
        }
    }
    if (p.x >= 1 && p.x <= s->x) {
        seen = 0;
        for (uint qy = 1; qy <= s->y; ++qy) {
            q = (t_point){ p.x, qy };
            if (qy == p.y || !shape_test_point(s, q))
                continue;
            ++seen;
            if (seen >= opt_n)
                return 0;
        }
    }
    if (p.y >= 1 && p.y <= s->y) {
        seen = 0;
        for (uint qx = 1; qx <= s->x; ++qx) {
            q = (t_point){ qx, p.y };
            if (qx == p.x || !shape_test_point(s, q))
                continue;
            ++seen;
            if (seen >= opt_n)
                return 0;
        }
    }

    return 1;
}

/* TODO: try to make a sane guess at number of records given data size,
 * to reduce the number of reallocs needed.
 */
size_t buf_heuristic(size_t size) {
    return 100;
}

t_buf *buf_new(size_t size) {
    t_buf *b = malloc(sizeof(t_buf));
    b->buf = malloc(size);
    b->bsize = size;
    b->csize = buf_heuristic(size);
    b->ptr = malloc(sizeof(ssize_t) * b->csize);
    b->bcount = 0;
    b->ccount = 0;
    return b;
}

void buf_free(t_buf *b) {
    free(b->ptr);
    free(b->buf);
    free(b);
}

void buf_grow(t_buf *b, size_t extra) {
    if (b->ccount >= b->csize) {
        size_t newsize = b->csize * 3 / 2;
        b->ptr = realloc(b->ptr, sizeof(ssize_t) * newsize);
        b->csize = newsize;
    }
    if (b->bcount + extra > b->bsize) {
        size_t newsize = b->bsize * 3 / 2;
        if (b->bcount + extra > b->bsize)
            newsize = b->bcount + extra;
        b->buf = realloc(b->buf, newsize);
        b->bsize = newsize;
    }
}

/* GNU stdlib offers qsort_r, but it is not portable */
uchar *sort_buf_nonreentrant;
#define bytecmp(x) if (a[x] != b[x]) return (a[x] < b[x]) ? -1 : 1;
int canon_comparator(const void *va, const void *vb) {
    uchar *context = sort_buf_nonreentrant;
    uchar *a = &context[*(ssize_t *)va];
    uchar *b = &context[*(ssize_t *)vb];
    bytecmp(1);
    bytecmp(2);
    uint p2size = raw_pack2_size(a[1], a[2]);
    for (uint i = 0; i < p2size; ++i)
        bytecmp(3 + i);
    /* shapes are canonically identical, order by index */
    for (uint i = 0; i < 5; ++i)
        bytecmp(3 + p2size + i);
    fprintf(stderr, "panic: comparing items with same index\n");
    exit(1);
}

/* restore order of pointers, moving NULLed pointers to the end */
int ptr_comparator(const void *va, const void *vb) {
    ssize_t a = *(ssize_t *)va;
    ssize_t b = *(ssize_t *)vb;
    if (a == (ssize_t)-1)
        return (b == (ssize_t)-1) ? 0 : 1;
    if (b == (ssize_t)-1)
        return -1;
    return (a < b) ? -1 : 1;
}

void diag_buf(t_buf *b, char *legend, bool full) {
    fprintf(stderr, "%s\n", legend);
    if (full) {
        for (size_t i = 0; i < b->bcount; ++i)
            fprintf(stderr, " %02x", b->buf[i]);
        fprintf(stderr, "\n");
    }
    for (size_t i = 0; i < b->ccount; ++i)
        fprintf(stderr, " %zd", b->ptr[i]);
    fprintf(stderr, "\n");
}

void buf_uniq(t_buf *b) {
    size_t count = b->ccount;
    uchar *buf = b->buf;
    sort_buf_nonreentrant = buf;
    qsort(b->ptr, count, sizeof(ssize_t), canon_comparator);
    for (size_t i = 1; i < b->ccount; ++i) {
        ssize_t s1 = b->ptr[i - 1];
        ssize_t s2 = b->ptr[i];
        if (buf[s1 + 1] != buf[s2 + 1] || buf[s1 + 2] != buf[s2 + 2])
            continue;
        uint p2size = raw_pack2_size(buf[s1 + 1], buf[s1 + 2]);
        if (memcmp(&buf[s1 + 3], &buf[s2 + 3], p2size))
            continue;
        /* keep first pointer, but move it down so loop can continue */
        b->ptr[i] = s1;
        b->ptr[i - 1] = (ssize_t)-1;
        --count;
    }
    qsort(b->ptr, b->ccount, sizeof(ssize_t), ptr_comparator);
    b->ccount = count;
}

void buf_write(t_buf *b, t_writer *w) {
    uchar *buf = b->buf;
    for (size_t i = 0; i < b->ccount; ++i) {
        uchar *s = &buf[b->ptr[i]];
        uint x = s[1], y = s[2];
        uint rsize = raw_rec_size(x, y, w->type);
        switch (w->type) {
          case REC_C:
            /* ok, write rsize bytes from s */
            break;
          case REC_B:
            /* skip (x, y, pack2) */
            s += 2 + raw_pack2_size(x, y);
            s[0] = (uchar)REC_B;
            break;
          case REC_A:
            /* skip (x, y, pack2, index) */
            s += 2 + raw_pack2_size(x, y) + 5;
            s[0] = (uchar)REC_A;
            break;
          default:
            fprintf(stderr, "panic: buf_write got type %02x\n", w->type);
            exit(1);
        }
        if (fwrite(s, 1, rsize, w->f) != rsize) {
            fprintf(stderr, "buf write record %zu for %u of %02x: %s\n",
                    i, rsize, w->type, strerror(errno));
            exit(1);
        }
    }
}

void buffer_shape_canonical(t_buf *b, t_shape *s) {
    buf_grow(b, shape_rec_size(s, REC_C));
    uchar *bp = b->buf;
    size_t off = b->bcount;
    b->ptr[b->ccount++] = (ssize_t)off;
    bp[off++] = (uchar)REC_C;
    off += shape_write_canonical(s, &bp[off]);
    off += shape_write_index(s, &bp[off]);
    off += shape_write_data(s, &bp[off]);
    b->bcount = off;
}

ulong gen_appended(uint k) {
    t_reader *r = reader_previous(k);
    t_writer *w = writer_new(k, REC_A);
    ulong count = 0;
    t_shape *s;
    while (s = shape_fetch(r)) {
        assert(bitcount(shape_poly_size(s), s->points) == k - 1);
        t_iter *i = shape_neighbours(s);
        while (1) {
            t_point p = shape_iter(i);
            if (p.x == 0 && p.y == 0)
                break;
            if (test_colinear(s, p)) {
                /* ok */
            } else {
                shape_mark_disallowed(s, p);
                iter_remove(i, p);
            }
        }
        iter_reset(i);
        while (1) {
            t_point p = shape_iter(i);
            if (p.x == 0 && p.y == 0)
                break;
            shape_mark_disallowed(s, p);
            t_shape *s2 = shape_append(s, p);
            assert(bitcount(shape_poly_size(s2), s2->points) == k);
            shape_write(s2, w);
            shape_free(s2);
            ++count;
        }
        iter_free(i);
    }
    writer_close(w);
    reader_close(r);
    return count;
}

void uniq_direct(uint k, t_reader *r, t_writer *w) {
    size_t count = 0;
    t_shape *s;
    t_buf *b = buf_new(reader_size(r) * 2);
    while (s = shape_fetch(r)) {
        s->index = count++;
        buffer_shape_canonical(b, s);
    }
    printf("k=%u: read %zu unsorted records (%.2fs)\n", k, count, elapsed());
    buf_uniq(b);
    printf("k=%u: writing %zu unique records (%.2fs)\n", k, b->ccount, elapsed());
    buf_write(b, w);
    buf_free(b);
}

void gen_uniq(uint k) {
    t_reader *r = reader_cur(k);
    size_t size = reader_size(r);
    if (size <= MAX_SORT) {
        t_writer *w = writer_uniq(k);
        uniq_direct(k, r, w);
        writer_close(w);
    } else {
        fprintf(stderr, "gen_uniq todo\n");
        exit(1);
    }
    reader_close(r);
}

void run(void) {
    for (uint k = 1; 1; ++k) {
        ulong count = gen_appended(k);
        if (count == 0)
            break;
        gen_uniq(k);
    }
}

void init(void) {
    init_rev();
}

void done(void) {
    ;
}

/* We run in multiple modes:
 *  - normal: ignore/replace any existing data
 *  - recover: continue from where we left off
 *  - check: report known results
 *
 * When running normally (or with recovery), there are multiple phases.
 * For polyominoes of size k:
 * - read previously generated polyominoes of size k-1, find each way to
 *   add another connected point, write results in condensed non-canonical
 *   form (type A);
 * - set that data as the input
 * - while the input is too large to handle in memory
 *   - if type A, split to multiple files on canonicalized size with index
 *     (type B)
 *   - if type B, canonicalize and split to multiple files on first byte
 *     (type C)
 *   - if type C, split to multiple files on next row of canonical data
 *     (type C)
 * - read the input (canonicalizing if type A or B), sort on canonical,
 *   uniq, re-sort in index order, writing type B if more merging is
 *   needed else type A;
 * - merge type B results maintaining the index sort, writing type B if
 *   more merging is needed else type A.
 * - profit!
 */
int main(int argc, char **argv) {
    uint argi = 1;
    while (argi < argc && *argv[argi] == '-') {
        char *arg = argv[argi++];
        if (strcmp(arg, "--") == 0)
            break;
        if (strcmp(arg, "-c") == 0) {
            opt_check = 1;
            continue;
        }
        if (strcmp(arg, "-r") == 0) {
            opt_recover = 1;
            continue;
        }
        fprintf(stderr, "Unknown option '%s'\n", arg);
        exit(1);
    }
    if (argi + 1 != argc) {
        fprintf(stderr, "Usage: %s <options> <n>\n", argv[0]);
        exit(1);
    }
    opt_n = strtoul(argv[argi], 0, 10);
    if (opt_n < 2 || opt_n > 5) {
        fprintf(stderr, "Expected 2 <= n <= 5, not %u\n", opt_n);
        exit(1);
    }
    init();
    if (opt_check) {
        t_reader r = { 0, REC_A, stdin };
        t_shape *s;
        size_t count = 0;
        while (s = shape_fetch(&r)) {
            fprintf(stderr, "index %zu", count++);
            diag_shape(s, ":");
            shape_free(s);
        }
    } else {
        run();
    }
    done();
    return 0;
}
