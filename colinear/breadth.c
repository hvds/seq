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

#include "poly.h"

uint opt_n;
bool opt_check = 0;
bool opt_recover = 0;
char path_buf[4096];
/* we probably want a size around 1Gb, but will test with 10kb */
/* #define MAX_SORT 10240 */
#define MAX_SORT 1073741824

struct stat statbuf;
t_shape *fetched, *appended, *emptied;
t_iter *ga_iter;

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

size_t pack3_size(t_dim d) {
    uint bits = (d.x + 2) * (d.y + 2);
    return (bits + 4) / 5;
}

size_t rec_size(t_dim d, e_rec type) {
    uint size = 1;  /* type marker */
    switch (type) {
        default:
            fprintf(stderr, "panic, unknown type %02x\n", type);
            exit(1);
        case REC_C:
            size += 2;                      /* canonical size */
            /* pack2 size is symmetric in (x, y) */
            size += pack2_size(d);          /* canonical data */
        case REC_B:
            size += 5;                      /* index */
        case REC_A:
            size += 2;                      /* oriented size */
            size += pack3_size(d);          /* full data */
        case REC_ALIGN:
            ;                               /* align is type marker only */
    }
    return size;
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
    s->d = (t_dim){ in[0], in[1] };
}

void read_pack3(t_reader *r, t_shape *s) {
    uint p3size = pack3_size(s->d);
    uchar in[p3size];
    if (fread(in, 1, p3size, r->f) != p3size)
        read_error(r, "pack3");
    uint size = poly_size(s->d);
    bzero(s->points, size);
    bzero(s->disallowed, size);
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
        if (++p.y > s->d.y + 1) {
            p.y = 0;
            if (++p.x > s->d.x + 1)
                break;
        }
    }
}

t_shape *shape_fetch(t_reader *r) {
    if (r->fake) {
        if (r->fake == 1) {
            r->fake = 2;
            shape_empty(fetched);
            return fetched;
        }
        return NULL;
    }
    if (!read_head(r))
        return NULL;

    t_shape *s = fetched;
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

void write_xy(t_writer *w, t_shape *s) {
    uchar out[2];
    out[0] = (uchar)s->d.x;
    out[1] = (uchar)s->d.y;
    if (fwrite(out, 1, 2, w->f) != 2) {
        fprintf(stderr, "write_xy: %s\n", strerror(errno));
        exit(1);
    }
}

void write_pack3(t_writer *w, t_shape *s) {
    uint size = pack3_size(s->d);
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
        if (++p.y > s->d.y + 1) {
            p.y = 0;
            if (++p.x > s->d.x + 1)
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
    p[off++] = s->d.x;
    p[off++] = s->d.y;
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
        if (++pt.y > s->d.y + 1) {
            pt.y = 0;
            if (++pt.x > s->d.x + 1)
                break;
        }
    }
    return (rot) ? off + 1 : off;
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
    uint p2size = pack2_size((t_dim){ a[1], a[2] });
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
        uint p2size = pack2_size((t_dim){ buf[s1 + 1], buf[s1 + 2] });
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
        t_dim d = (t_dim){ s[1], s[2] };
        uint rsize = rec_size(d, w->type);
        switch (w->type) {
          case REC_C:
            /* ok, write rsize bytes from s */
            break;
          case REC_B:
            /* skip (x, y, pack2) */
            s += 2 + pack2_size(d);
            s[0] = (uchar)REC_B;
            break;
          case REC_A:
            /* skip (x, y, pack2, index) */
            s += 2 + pack2_size(d) + 5;
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
    buf_grow(b, rec_size(s->d, REC_C));
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
        assert(bitcount(poly_size(s->d), s->points) == k - 1);
        shape_neighbours(ga_iter, s);
        t_iter *i = ga_iter;
        while (1) {
            t_point p = shape_iter(i);
            if (p.x == 0 && p.y == 0)
                break;
            if (test_colinear(s, p, opt_n)) {
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
            t_shape *s2 = shape_append(appended, s, p);
            assert(bitcount(poly_size(s2->d), s2->points) == k);
            shape_write(s2, w);
            ++count;
        }
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
    init_poly();
    ga_iter = malloc(sizeof(t_iter) + max_poly_size());
    fetched = shape_new((t_dim){ MAXDIM, MAXDIM });
    appended = shape_new((t_dim){ MAXDIM, MAXDIM });
}

void done(void) {
    shape_free(appended);
    shape_free(fetched);
    free(ga_iter);
    done_poly();
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
        }
    } else {
        run();
    }
    done();
    return 0;
}
