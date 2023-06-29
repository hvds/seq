#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>
#include <sys/resource.h>

typedef unsigned char bool;

typedef uint t_vec;     /* any unsigned integer type should work */
t_vec vmax = 0;         /* set during init() */

/* Count the bits set in v1 */
static inline uint vbits(t_vec v1) {
    return __builtin_popcount(v1);
}

/* Return index of first one-bit in v1.
 * Note: this is undefined for v1 == 0. */
static inline uint vfirst(t_vec v1) {
    return __builtin_ctz(v1);
}

/* Return t_vec with only bit vfirst(v1) set.
 * Note: this is undefined for v1 == 0. */
static inline t_vec vfirstv(t_vec v1) {
    /* CHECKME: which is faster? */
#if 0
    return v1 & (v1 ^ (v1 - 1));
#else
    return 1 << vfirst(v1);
#endif
}

static inline int vcmp(const void *vp1, const void *vp2) {
    const t_vec v1 = *(const t_vec *)vp1;
    const t_vec v2 = *(const t_vec *)vp2;
    return (v1 > v2) ? 1 : (v1 < v2) ? -1 : 0;
}

void set_vmax(uint d) {
    if (d > sizeof(t_vec) * 8) {
        fprintf(stderr, "Error: too many dimensions for this build\n");
        exit(1);
    }
    if (d == sizeof(t_vec) * 8) {
        vmax = ~(t_vec)0;
    } else {
        vmax = ((t_vec)1 << d) - 1;
    }
}

/* set to utime at start of run */
double t0 = 0;
struct rusage rusage_buf;
static inline double utime(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000 - t0;
}

uint n;             /* number of points */
uint d;             /* number of dimensions */
uint best = 0;      /* highest number of squares seen for any set */
ulong recurse_iter = 0;  /* number of iterations */
typedef struct s_v {
    t_vec v;        /* the value */
    t_vec group;    /* bits set for all grouped dimensions */
    t_vec gfirst;   /* bits set for start dimension of each group */
    t_vec best;     /* value in best set */
} t_v;
t_v *v = NULL;    /* The points themselves */

static inline int vcmpv(const void *vp1, const void *vp2) {
    return vcmp(&((t_v *)vp1)->v, &((t_v *)vp2)->v);
}

/* after n and d are known */
void init(void) {
    t0 = utime();
    set_vmax(d);
    v = calloc(n, sizeof(t_v));
}

/* Returns true if point vp is in ordered list v[] in the range start..n-1.
 */
bool findv(t_vec p, uint start, uint end) {
    uint low = start;
    uint high = end - 1;
    while (low <= high) {
        uint m = (low + high) >> 1;
        if (v[m].v > p)
            high = m - 1;
        else if (v[m].v < p)
            low = m + 1;
        else
            return 1;
    }
    return 0;
}

/* Return the number of squares generated by the n d-dimensional points
 * in v[].
 */
uint count_squares(void) {
    uint count = 0;
    for (uint i = 0; i < n - 3; ++i) {
        t_vec vi = v[i].v;
        for (uint j = i + 1; j < n - 2; ++j) {
            t_vec vj = v[j].v;
            t_vec xij = vi ^ vj;
            uint cij = vbits(xij);
            for (uint k = j + 1; k < n - 1; ++k) {
                t_vec vk = v[k].v;
                t_vec xik = vi ^ vk;
                /* if xij and xik have bits in common, then JIK is not a right
                 * angle */
                if (xij & xik)
                    continue;
                uint cik = vbits(xik);
                /* if cij and cik differ then IJ and IK have different length */
                if (cij != cik)
                    continue;
                /* now if the fourth point is present, it makes a square */
                if (findv(vk ^ xij, k + 1, n)) {
                    ++count;
#if 0
                    printf("%2$0*1$x %3$0*1$x %4$0*1$x %5$0*1$x\n",
                            (d + 3) / 4, vi, vj, vk, vk ^ xij);
#endif
                }
            }
        }
    }
    return count;
}

void recurse(void) {
    uint sp = 0, dmin = 0;
    t_v *prev, *cur;
    t_vec next;

    /* first element is always (t_vec)0, leaving a single group of
     * dimensions starting at d_0 = 1, extending through all dimensions */
    v[0].v = (t_vec)0;
    v[0].gfirst = (t_vec)1;
    v[0].group = vmax;

    while (1) {
        ++sp;
        if (sp == n) {
            uint count = count_squares();
            if (count > best) {
                best = count;
#if 0
                printf("%u:", best);
                for (uint i = 0; i < n; ++i) {
                    v[i].best = v[i].v;
                    printf(" %u", v[i].v);
                }
                printf(" (%.2fs)\n", utime());
#else
                for (uint i = 0; i < n; ++i)
                    v[i].best = v[i].v;
#endif
            }
            goto derecurse;
        }
        prev = &v[sp - 1];
        cur = &v[sp];
        cur->v = prev->v;
        goto next_value;
      derecurse:
        --sp;
        if (sp == 0)
            break;
        prev = &v[sp - 1];
        cur = &v[sp];
        goto next_value;
      retry:
        cur->v = next;
      next_value:
        ++recurse_iter;
        if (cur->v == vmax)
            goto derecurse;
        next = cur->v + 1;
      retry_value: ;
        uint bits = vbits(next);
        if (bits < dmin)
            goto retry;
        /* eg {0,1,3} can be reflected to {0,1,2} */
        if (dmin == 1 && bits == 2 && (next & 1) && sp > 1
                && !findv(next ^ 1, 2, sp))
            goto retry;
        for (uint i = 1; i < sp; ++i)
            if (vbits(next ^ v[i].v) < dmin)
                goto retry;
        t_vec group = prev->group;
        t_vec gfirst = prev->gfirst;
        t_vec grouped = next & group;
        while (grouped) {
            t_vec gbu = vfirstv(grouped);
            t_vec gbl = gbu;
            t_vec gb = gbu;
            grouped = grouped & ~gbu;
            while ((gfirst & gb) == 0) {
                gbl >>= 1;
                gb |= gbl;
            }
            if (gb == gbu) {
                group = group & ~gbu;
                gfirst = gfirst & ~gbu;
            } else {
                next = (next | gb) & ~(gbl - 1);
                goto retry_value;
            }
            /* FIXME: do this if (gb << 2) is set in group and clear in gfirst,
             * or if (gb << 2) is out of range; else clear (gb << 1) in group */
            gfirst = gfirst | (group & (gbu << 1));
        }
        cur->v = next;
        cur->group = group;
        cur->gfirst = gfirst;
#if 0
for (uint i = 0; i <= sp; ++i) { printf("[%u %u %u] ", v[i].v, v[i].gfirst, v[i].group); } printf("\n");
#endif
        if (sp == 1)
            dmin = vbits(next);
        /* loop to select next point */
    }
    return;
}

int main(int argc, char **argv) {
    int argi = 1;
    bool do_count_squares = 0;
    uint base = 10;
    while (argi < argc) {
        char *arg = argv[argi];
        if (arg[0] != '-')
            break;
        ++argi;
        if (strcmp(arg, "--") == 0)
            break;
        if (strcmp(arg, "-c") == 0) {
            do_count_squares = 1;
            continue;
        }
        if (strncmp(arg, "-b", 2) == 0) {
            base = atoi(&arg[2]);
            continue;
        }
        if (strncmp(arg, "-d", 2) == 0) {
            d = atoi(&arg[2]);
            continue;
        }
        fprintf(stderr, "Unknown option '%s'\n", arg);
        return 1;
    }
    if (do_count_squares) {
        n = argc - argi;
        if (d == 0)
            d = sizeof(t_vec) * 8;
        init();
        for (uint i = 0; i < n; ++i)
            v[i].v = strtol(argv[i + argi], NULL, base);
        qsort(&v[0], n, sizeof(t_v), &vcmpv);
/*
        for (uint i = 0; i < n; ++i) { printf("%u ", v[i].v); } printf("\n");
*/
        printf("%u\n", count_squares());
        return 0;
    }

    if (argc - argi != 2) {
        fprintf(stderr, "Wrong number of arguments\n");
        return 1;
    }
    n = atoi(argv[argi++]);
    d = atoi(argv[argi++]);
    init();
    recurse();
    printf("f(%u, %u) = %u:", n, d, best);
    for (uint i = 0; i < n; ++i)
        printf(" %u", v[i].best);
    printf(" (%lu iter %.2fs)", recurse_iter, utime());
    printf("\n");
    return 0;
}
