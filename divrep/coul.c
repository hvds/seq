#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#ifdef HAVE_SETPROCTITLE
#   include <sys/types.h>
#endif
#include <sys/times.h>

#include "coul.h"
#include "coulfact.h"
#include "diag.h"
#include "rootmod.h"
#include "coultau.h"
#include "pell.h"

/* from MPUG */
#include "factor.h"
#include "gmp_main.h"
#include "utility.h"
#include "primality.h"

/* primary parameters - we are searching for D(n, k), the least d such
 * that tau(d + i) = n for all 0 <= i < k.
 */
uint n, k;

/* mpz_t passed as function parameter decays to pointer in a way that
 * allows it to be used as mpz_t, but cannot be converted to a pointer
 * in a typesafe manner. Given a function called as foo(z), use this as
 *   mpz_t *zp = PARAM_TO_PTR(z);
 */
static inline mpz_t *PARAM_TO_PTR(__mpz_struct *z) {
    return (mpz_t *)z;
}

/* stash of mpz_t, initialized once at start */
typedef enum {
    zero, zone,                 /* constants */
    sqm_t, sqm_q, sqm_b, sqm_z, sqm_x,  /* sqrtmod_t */
    uc_minusvi, uc_px,          /* update_chinese */
    ur_a, ur_m, ur_ipg,         /* update_residues */
    asq_o, asq_qq, asq_m,       /* alloc_square */
    wv_ati, wv_end, wv_cand,    /* walk_v */
    wv_startr, wv_endr, wv_qqr, wv_r, wv_rx, wv_temp,
    wv_x, wv_y, wv_x2, wv_y2,
    w1_v, w1_j, w1_r,           /* walk_1 */
    lp_x, lp_mint,              /* limit_p */
    r_walk,                     /* recurse */

    sdm_p, sdm_r,               /* small_divmod (TODO) */
    dm_r,                       /* divmod */
    np_p,                       /* next_prime (TODO) */

    MAX_ZSTASH
} t_zstash;
mpz_t *zstash;
static inline mpz_t *ZP(t_zstash e) { return &zstash[e]; }
#define Z(e) *ZP(e)
/* additional arrays of mpz_t initialized once at start */
mpz_t *wv_o = NULL, *wv_qq = NULL;  /* wv_o[k], wv_qq[k] */

typedef unsigned char bool;

/* used to store disallowed inverses in walk_v() */
typedef struct s_mod {
    ulong v;
    ulong m;
} t_mod;

t_divisors *divisors = NULL;

/* For prime p < k, we "force" allocation of powers in a batch to ensure
 * that multiple allocations of the same prime are coherent. Whereas normal
 * allocation considers only p^{x-1} where x is divisible by the highest
 * prime dividing t_i, forced primes may allocate any p^{x-1} with x | t_i.
 *
 * For cases where two or more v_i are divisible by p, we always force
 * every possible case. For cases where only one v_i is divisible by p,
 * we force them only if n == 2 (mod 4) (heuristically, since allocations
 * are always either at least as powerful as normal allocations _or_ they
 * have x=2, leaving v_i.t odd) or if requested by the -f option (force_all).
 *
 * Each batch describes the location and magnitude of the highest power of p
 * that divides any v_i (using the lower i if there are more than one).
 */
typedef struct s_forcebatch {
    uint vi;    /* allocate p^{x-1} at v_{vi} */
    uint x;
} t_forcebatch;
typedef struct s_forcep {
    uint p;
    uint count;
    t_forcebatch *batch;
} t_forcep;
t_forcep *forcep = NULL;
uint forcedp;
uint force_all = 0;

/* When allocation forces the residue of some v_i to be square, we want
 * to calculate the roots mod every allocation (before and after this one),
 * first to check if a solution is possible, and second to avoid duplicate
 * work when we actually use the roots in walk_v().
 * The set of roots lives in resarray(level), but here we track what power
 * root they are: the allocation at value[sq0].alloc[i] leaves gcddm sqg[i].
 */
uint sq0 = 0;
uint *sqg = NULL;   /* size maxfact */

/* Each level of "recursion" allocates one prime power p^{x-1} with x | n
 * to one value v_i. It may be "forced", in which case it is part of a
 * batch of simultaneous allocations of the same p to different i (and
 * in which case derecursion should unwind the whole batch), or "unforced",
 * in which case no other v_j will be divisible by p.
 */
typedef struct s_level {
    uint level;     /* index of this entry */
    bool is_forced;
    uint vi;        /* allocation of p^x into v_i */
    ulong p;
    uint x;
    uint have_square;   /* number of v_i residues forced square so far */
    uint nextpi;    /* index of least prime not yet allocated */
    ulong maxp;     /* highest prime allocated so far */
    /* (optional) union */
        uint bi;    /* batch index, if forced */
    /* .. with */
        uint di;    /* divisor index, if unforced */
        uint ti;    /* target tau for v_i */
        ulong limp; /* limit for p */
        uint max_at;/* which max value used for limp calculation */
    /* end union */
    mpz_t aq;       /* running LCM of allocated p^x */
    mpz_t rq;       /* running CRT of (-i) % p^x */
} t_level;
t_level *levels = NULL;
uint level = 0;

/* Each value v_0 to v_{k-1} has had 'level' allocations of prime powers
 * p_i^{x_i-1}; q tracks the product of those prime powers, and t tracks
 * our target tau - starting at n, and divided by x_i on each allocation.
 */
typedef struct s_allocation {
    ulong p;
    uint x;
    mpz_t q;
    uint t;
} t_allocation;
typedef struct s_value {
    uint vlevel;
    t_allocation *alloc;    /* size maxfact */
} t_value;
t_value *value = NULL;

/* saved counts of allocations in each value before applying nth forced prime */
uint *vlevels = NULL;
uint vl_forced = 0;
static inline uint *VLP(uint vlevel) { return &vlevels[vlevel * k]; }
static inline void STOREVL(uint vli) {
    uint *vlp = VLP(vli);
    for (uint vi = 0; vi < k; ++vi)
        vlp[vi] = value[vi].vlevel;
}
static inline void FETCHVL(uint vli) {
    uint *vlp = VLP(vli);
    for (uint vi = 0; vi < k; ++vi)
        value[vi].vlevel = vlp[vi];
}

/* list of some small primes, at least enough for one per allocation  */
uint *sprimes = NULL;
uint nsprimes;

long ticks_per_second;
/* set to utime at start of run, minus last timestamp of recovery file */
clock_t ticks = 0;
struct tms time_buf;
static inline clock_t utime(void) {
    times(&time_buf);
    return time_buf.tms_utime;
}

mpz_t min, max;     /* limits to check for v_0 */
uint seen_best = 0; /* number of times we've improved max */
ulong gain = 0;     /* used to fine-tune balance of recursion vs. walk */
ulong antigain = 0;
/* maxp is the greatest prime we should attempt to allocate; minp is the
 * threshold that at least one allocated prime should exceed (else we can
 * skip the walk)
 */
uint minp = 0, maxp = 0;
uint rough = 0;     /* test roughness if tau >= rough */
bool opt_print = 0; /* print candidates instead of fully testing them */
/* If opt_alloc is true and opt_batch < 0, just show the forced-prime
 * allocation; if opt_alloc is true and opt_batch >= 0, just process
 * the specified batch_alloc.
 */
bool opt_alloc = 0;
int opt_batch = -1;
int batch_alloc = 0;   /* index of forced-prime allocations */

int debug = 0;     /* diag and keep every case seen */
ulong randseed = 1; /* for ECM, etc */

char *rpath = NULL; /* path to log file */
FILE *rfp = NULL;   /* file handle to log file */
bool start_seen = 0;    /* true if log file has been written to before */
t_fact *rstack = NULL;  /* point reached in recovery log file */
bool have_rwalk = 0;    /* true if recovery is mid-walk */
mpz_t rwalk_from;
mpz_t rwalk_to;

t_fact nf;      /* factors of n */
uint tn;        /* tau(n) */
uint maxfact;   /* count of prime factors dividing n, with multiplicity */
uint maxodd;    /* as above for odd prime factors */
uint *maxforce = NULL;  /* max prime to force at v_i */
mpz_t px;       /* p^x */

#define DIAG 1
#define LOG 600
clock_t diag_delay, log_delay, diagt, logt;
ulong countr, countw, countwi;
#define DIAG_BUFSIZE (3 + k * maxfact * (20 + 1 + 5 + 1) + 1)
char *diag_buf = NULL;

void prep_show_v(void) {
    uint offset = 0;
    for (uint vi = 0; vi < k; ++vi) {
        t_value *vp = &value[vi];
        if (vi)
            diag_buf[offset++] = ' ';
        if (vp->vlevel == 0)
            diag_buf[offset++] = '.';
        else {
            for (uint ai = 0; ai < vp->vlevel; ++ai) {
                t_allocation *ap = &vp->alloc[ai];
                if (ai)
                    diag_buf[offset++] = '.';
                offset += sprintf(&diag_buf[offset], "%lu", ap->p);
                if (ap->x > 2)
                    offset += sprintf(&diag_buf[offset], "^%u", ap->x - 1);
            }
        }
    }
    diag_buf[offset] = 0;
}

void report(char *format, ...) {
    keep_diag();
    va_list ap;
    va_start(ap, format);
    gmp_vfprintf(stdout, format, ap);
    va_end(ap);

    if (rfp) {
        va_start(ap, format);
        gmp_vfprintf(rfp, format, ap);
        va_end(ap);
        fflush(rfp);
    }
}

double seconds(clock_t t1) {
    return (double)(t1 - ticks) / ticks_per_second;
}

void diag_plain(void) {
    clock_t t1 = utime();

    prep_show_v();  /* into diag_buf */
    diag("%s", diag_buf);
    if (debug)
        keep_diag();
    diagt = t1 + diag_delay;

    if (rfp && t1 >= logt) {
        fprintf(rfp, "305 %s (%.2fs)\n", diag_buf, seconds(t1));
        logt = t1 + log_delay;
    }
}

void diag_walk_v(ulong ati, ulong end) {
    clock_t t1 = utime();

    prep_show_v();  /* into diag_buf */
    if (!(debug == 1 && ati))
        diag("%s: %lu / %lu", diag_buf, ati, end);
    if (debug)
        keep_diag();
    diagt = t1 + diag_delay;

    if (rfp && t1 >= logt) {
        fprintf(rfp, "305 %s: %lu / %lu (%.2fs)\n",
                diag_buf, ati, end, seconds(t1));
        logt = t1 + log_delay;
    }
}

void diag_walk_zv(mpz_t ati, mpz_t end) {
    clock_t t1 = utime();

    prep_show_v();  /* into diag_buf */
    if (!(debug && mpz_sgn(ati)))
        diag("%s: %Zu / %Zu", diag_buf, ati, end);
    if (debug)
        keep_diag();
    diagt = t1 + diag_delay;

    if (rfp && t1 >= logt) {
        gmp_fprintf(rfp, "305 %s: %Zu / %Zu (%.2fs)\n",
                diag_buf, ati, end, seconds(t1));
        logt = t1 + log_delay;
    }
}

void diag_walk_pell(uint pc) {
    clock_t t1 = utime();

    prep_show_v();  /* into diag_buf */
    if (!(debug && pc))
        diag("%s: P%u", diag_buf, pc);
    if (debug)
        keep_diag();
    diagt = t1 + diag_delay;

    if (rfp && t1 >= logt) {
        fprintf(rfp, "305 %s: P%u (%.2fs)\n",
                diag_buf, pc, seconds(t1));
        logt = t1 + log_delay;
    }
}

void disp_batch(t_level *lp) {
    prep_show_v();      /* into diag_buf */
    printf("203 b%u: %s", batch_alloc, diag_buf);
    if (lp->have_square)
        printf(" [sq=%u]", lp->have_square);
    printf("\n");
}

void candidate(mpz_t c) {
    keep_diag();
    clock_t t1 = utime();
    report("202 Candidate %Zu (%.2fs)\n", c, seconds(t1));
    if (mpz_cmp(c, max) <= 0) {
        mpz_set(max, c);
        ++seen_best;
    }
}

void free_levels(void) {
    for (uint i = 0; i < k * maxfact + 1; ++i) {
        t_level *l = &levels[i];
        mpz_clear(l->aq);
        mpz_clear(l->rq);
    }
    free(levels);
}

void init_levels(void) {
    /* CHECKME: can this be maxodd * k + forcedp? */
    levels = (t_level *)calloc(k * maxfact + 1, sizeof(t_level));
    for (uint i = 0; i < k * maxfact + 1; ++i) {
        t_level *l = &levels[i];
        l->level = i;
        mpz_init(l->aq);
        mpz_init(l->rq);
    }
    mpz_set_ui(levels[0].aq, 1);
    mpz_set_ui(levels[0].rq, 0);
    levels[0].have_square = 0;
    levels[0].nextpi = 0;
    levels[0].maxp = 0;
    level = 1;
}

void free_value(void) {
    for (int i = 0; i < k; ++i) {
        t_value *v = &value[i];
        for (int j = 0; j < maxfact; ++j)
            mpz_clear(v->alloc[j].q);
        free(v->alloc);
    }
    free(value);
}

void init_value(void) {
    value = (t_value *)malloc(k * sizeof(t_value));
    for (int i = 0; i < k; ++i) {
        t_value *v = &value[i];
        v->vlevel = 0;
        v->alloc = (t_allocation *)malloc(maxfact * sizeof(t_allocation));
        for (uint j = 0; j < maxfact; ++j)
            mpz_init(v->alloc[j].q);
    }
}

void done(void) {
    free(diag_buf);
    if (wv_qq)
        for (uint i = 0; i < k; ++i)
            mpz_clear(wv_qq[i]);
    free(wv_qq);
    if (wv_o)
        for (uint i = 0; i < k; ++i)
            mpz_clear(wv_o[i]);
    free(wv_o);
    free(sqg);
    free(vlevels);
    free_value();
    free_levels();
    free(sprimes);
    if (forcep)
        for (int i = 0; i < forcedp; ++i)
            free(forcep[i].batch);
    free(forcep);
    free(maxforce);
    if (divisors)
        for (int i = 0; i <= n; ++i)
            free(divisors[i].div);
    free(divisors);
    if (have_rwalk) {
        mpz_clear(rwalk_from);
        mpz_clear(rwalk_to);
    }
    if (rstack)
        for (int i = 0; i < k; ++i)
            free_fact(&rstack[i]);
    free(rstack);
    if (rfp)
        fclose(rfp);
    free(rpath);
    for (t_zstash i = 0; i < MAX_ZSTASH; ++i)
        mpz_clear(Z(i));
    free(zstash);
    mpz_clear(px);
    free_fact(&nf);
    mpz_clear(max);
    mpz_clear(min);
    done_pell();
    done_rootmod();
    done_tau();
    _GMP_destroy();
}

void fail(char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    /* we accept leaks on fatal error, but should close the log file */
    if (rfp)
        fclose(rfp);
    exit(1);
}

void init_pre(void) {
    _GMP_init();
    /* reseed immediately to retain reproducibility */
    clear_randstate();
    init_randstate(1);
    /* we may do this again after options handled, to select real seed */

    init_pell();
    ticks_per_second = sysconf(_SC_CLK_TCK);
    ticks = utime();
    mpz_init_set_ui(min, 0);
    mpz_init_set_ui(max, 0);
    init_fact(&nf);
    mpz_init(px);
    zstash = (mpz_t *)malloc(MAX_ZSTASH * sizeof(mpz_t));
    for (t_zstash i = 0; i < MAX_ZSTASH; ++i)
        mpz_init(Z(i));
    mpz_set_ui(Z(zero), 0);
    mpz_set_ui(Z(zone), 1);
}

/* Parse a "305" log line for initialization.
 * Input string should point after the initial "305 ".
 */
void parse_305(char *s) {
    double dtime;
    t_ppow pp;

    rstack = (t_fact *)malloc(k * sizeof(t_fact));
    for (int i = 0; i < k; ++i)
        init_fact(&rstack[i]);

    for (int i = 0; i < k; ++i) {
        if (i) {
            assert(s[0] == ' ');
            ++s;
        }
        if (s[0] == '.') {
            ++s;
            continue;
        }
        while (1) {
            pp.p = strtoul(s, &s, 10);
            pp.e = (s[0] == '^') ? strtoul(&s[1], &s, 10) : 1;
            add_fact(&rstack[i], pp);
            if (s[0] != '.')
                break;
            ++s;
        }
        /* reverse them, so we can pop as we allocate */
        reverse_fact(&rstack[i]);
    }
    if (s[0] == ':') {
        assert(s[1] == ' ');
        s += 2;
        if (strncmp("t=1", s, 3) == 0)
            s += 3; /* ignore */
        else {
            int from_start, from_end, to_start, to_end;
            have_rwalk = 1;
            if (EOF == sscanf(s, "%n%*[0-9]%n / %n%*[0-9]%n ",
                    &from_start, &from_end, &to_start, &to_end))
                fail("could not parse 305 from/to: '%s'", s);
            s[from_end] = 0;
            mpz_init_set_str(rwalk_from, &s[from_start], 10);
            s[to_end] = 0;
            mpz_init_set_str(rwalk_to, &s[to_start], 10);
            have_rwalk = 1;
            s[to_end] = ' ';
            s = &s[to_end];
        }
    }
    if (EOF == sscanf(s, " (%lfs)\n", &dtime))
        fail("could not parse 305 time: '%s'", s);
    ticks -= (clock_t)dtime * ticks_per_second;
}

void recover(void) {
    char *last305 = NULL;
    char *last202 = NULL;
    char *curbuf = NULL;
    size_t len = 120, len305 = 0, len202 = 0;

    while (1) {
        ssize_t nread = getline(&curbuf, &len, rfp);
        if (nread < 0) {
            if (errno == 0)
                break;
            fail("error reading %s: %s", rpath, strerror(errno));
        }
        if (curbuf[nread - 1] != '\n'
                || memchr(curbuf, 0, nread) != NULL) {
            /* corrupt line, file should be truncated */
            off_t offset = ftello(rfp);
            if (offset == -1)
                fail("could not ask offset: %s", strerror(errno));
            if (ftruncate(fileno(rfp), offset - nread) != 0)
                fail("could not truncate %s to %lu: %s", rpath, offset - nread,
                        strerror(errno));
            if (freopen(NULL, "a+", rfp) == NULL)
                fail("could not reopen %s: %s", rpath, strerror(errno));
            break;
        }
        if (strncmp("305 ", curbuf, 4) == 0) {
            char *t = last305;
            last305 = curbuf;
            curbuf = t;
            size_t lt = len305;
            len305 = len;
            len = lt;
        } else if (strncmp("202 ", curbuf, 4) == 0) {
            char *t = last202;
            last202 = curbuf;
            curbuf = t;
            size_t lt = len202;
            len202 = len;
            len = lt;
        } else if (strncmp("001 ", curbuf, 4) == 0) {
            /* TODO: parse and check for consistent options */
            start_seen = 1;
        } else if (strncmp("000 ", curbuf, 4) == 0)
            ;   /* comment */
        else
            fail("unexpected log line %.3s in %s", curbuf, rpath);
    }
    if (last202) {
        int start, end;
        mpz_t cand;
        if (EOF == sscanf(last202, "202 Candidate %n%*[0-9]%n (%*[0-9.]s)\n",
                &start, &end))
            fail("error parsing 202 line '%s'", last202);
        last202[end] = 0;
        mpz_init_set_str(cand, &last202[start], 10);
        if (mpz_sgn(max) == 0 || mpz_cmp(max, cand) >= 0) {
            mpz_set(max, cand);
            ++seen_best;
        }
        mpz_clear(cand);
    }
    if (last305)
        parse_305(last305 + 4);
    free(curbuf);
    free(last305);
    free(last202);
}

int cmp_high(const void *va, const void *vb) {
    uint a = *(uint *)va, b = *(uint *)vb;
    return (int)divisors[b].high - (int)divisors[a].high;
}

/* TODO: use MPU code for ulong next_prime */
ulong next_prime(ulong cur) {
    mpz_set_ui(Z(np_p), cur);
    _GMP_next_prime(Z(np_p));
    if (mpz_fits_ulong_p(Z(np_p)))
        return mpz_get_ui(Z(np_p));
    diag_plain();
    keep_diag();
    report("002 next_prime overflow\n");
    exit(1);
}

/* recurse() wants the list of powers to try: each divisor of t_i (which
 * itself divides n) that is divisible by the highest odd prime factor
 * dividing t_i, in increasing order.
 * prep_forcep() wants the full list of divisors, but in similar order.
 * For each power, recurse() also wants to know which powers to skip
 * if the previous power was a given value, but that's simply:
 * skip x' if x' < x and high(x') == high(x).
 * mintau() wants sumpm, sum{p_j - 1} of the primes dividing t_i with
 * multiplicity.
 * When a square is fixed, walk_v() wants gcddm, the gcd{d_j-1} of all
 * divisors d_j of t_i.
 */
void prep_fact(void) {
    t_fact f;

    divisors = (t_divisors *)calloc(n + 1, sizeof(t_divisors));
    init_fact(&f);
    for (uint i = 1; i <= n; ++i) {
        if (n % i)
            continue;
        t_divisors *dp = &divisors[i];
        f.count = 0;
        simple_fact(i, &f);
        uint td = simple_tau(&f);
        dp->high = (f.count) ? f.ppow[f.count - 1].p : 1;
        dp->sumpm = dp->high - 1;
        dp->sumpm += divisors[i / dp->high].sumpm;
        uint nd = 0;
        dp->div = (uint *)malloc(td * sizeof(uint));
        for (uint j = 1; j <= i; ++j) {
            if (i % j)
                continue;
            dp->div[dp->alldiv++] = j;
            if ((j % dp->high) == 0)
                ++dp->highdiv;
        }
        qsort(dp->div, dp->alldiv, sizeof(uint), &cmp_high);

        uint g = dp->div[0] - 1;
        for (uint di = 1; di < dp->alldiv; ++di)
            g = tiny_gcd(g, dp->div[di] - 1);
        dp->gcddm = g;
    }
    free_fact(&f);
}

void prep_mintau(void) {
    /* max tau we care about */
    uint maxt = n / (nf.count ? nf.ppow[nf.count - 1].p : 1);

    /* max factors in a factorization of maxt */
    int cp = -1;
    for (int i = 0; i < nf.count; ++i)
        cp += nf.ppow[i].e;

    /* max allocation per element is number of odd prime factors */
    int cop = cp + 1 - (nf.count ? 0 : nf.ppow[0].e);

    /* max total allocations exclduding current one, ignoring forced primes */
    int max_alloc = cop * k - 1;

    /* we need at most 1 per factor of a factorization, plus 1 for each
     * allocation, plus 1 per forced prime (which can be used without
     * filling a normal allocation).
     * Note: prime_count(k) will be forcedp.
     * Note: this could be reduced by the min number of forced allocations
     * that will use allocatable factors (ie divisible by odd prime).
     */
    int need = cp + max_alloc + simple_prime_count(k);

    /* TODO: work out first what we want mintau() to do; we probably want
     * to fill this cache lazily, or not cache at all.
     */
}

void prep_maxforce(void) {
    maxforce = (uint *)malloc(k * sizeof(uint));
    if ((n & 3) != 0) {
        for (int i = 0; i < k; ++i)
            maxforce[i] = k;
    } else
        for (int i = 0; i < k; ++i) {
            int mf = k - i - 1;
            if (i > mf)
                mf = i;
            if (force_all > mf)
                mf = force_all;
            maxforce[i] = mf;
        }
}

void prep_primes(void) {
    /* We can certainly not allocate more than (each of the forced primes)
     * plus (one per odd prime factor for each v_i); in practice it will
     * usually be less. */
    nsprimes = maxodd * k + forcedp;
    sprimes = (uint *)malloc(nsprimes * sizeof(uint));
    uint p = 1;
    for (uint i = 0; i < nsprimes; ++i) {
        p = next_prime(p);
        sprimes[i] = p;
    }
}

void prep_forcep(void) {
    mpz_t pz;
    uint p;
    uint pi[k];

    forcedp = 0;
    mpz_init_set_ui(pz, 1);
    while (1) {
        _GMP_next_prime(pz);
        p = mpz_get_ui(pz);
        if (p > k)
            break;
        pi[forcedp++] = p;
    }
    mpz_clear(pz);

    uint first_bad = n + 1;
    for (uint div = 2; div < n; ++div)
        if (n % div) {
            first_bad = div;
            break;
        }

    forcep = (t_forcep *)malloc(forcedp * sizeof(t_forcep));
    t_divisors *d = &divisors[n];
    for (uint fpi = 0; fpi < forcedp; ++fpi) {
        t_forcep *fp = &forcep[fpi];
        fp->p = p = pi[fpi];
        fp->count = 0;
        fp->batch = (t_forcebatch *)malloc(tn * k * sizeof(t_forcebatch));

        uint x = 1, px = 1;
        while (px * p < k) {
            ++x;
            px *= p;
        }
        uint bad_pow = k + 1;
        if (first_bad <= x) {
            bad_pow = 1;
            for (uint i = 1; i < first_bad; ++i)
                bad_pow *= p;
        }

        bool keep_single = ((n & 3) || p <= force_all) ? 1 : 0;
        bool skipped = 0;
        for (uint vi = 0; vi < k; ++vi) {
            /* skip if forcing p^{x-1} at v_i would require p^{y-1} at v_j
             * such that not y | n.
             */
            if (vi + bad_pow < k || vi >= bad_pow)
                continue;
            if (p > maxforce[vi] && !keep_single) {
                skipped = 1;
                continue;
            }
            for (uint di = 0; di < d->alldiv; ++di) {
                uint fx = d->div[di];
                if (fx < x)
                    continue;   /* would not be highest power */
                if (fx == x) {
                    if (vi >= px)
                        continue;   /* would not be first highest power */
                    if (vi + px * (p - 1) < k)
                        continue;   /* would not be highest power */
                }
                fp->batch[fp->count++] = (t_forcebatch){ .vi = vi, .x = fx };
            }
        }
        if (fp->count == 0) {
            forcedp = fpi;
            free(fp->batch);
            break;
        }
        if (skipped)
            fp->batch[fp->count++] = (t_forcebatch){ .x = 0 };
        fp->batch = (t_forcebatch *)realloc(fp->batch,
                fp->count * sizeof(t_forcebatch));
    }
}

void init_post(void) {
    init_tau(rough);
    alloc_taum(k);
    if (randseed != 1) {
        /* hard to guarantee we haven't used any randomness before this.
         * note also that this will give different results for a run that
         * is stopped and recovered.
         */
        clear_randstate();
        init_randstate(randseed);
    }
    if (rpath) {
        printf("path %s\n", rpath);
        rfp = fopen(rpath, "a+");
        if (rfp == NULL)
            fail("%s: %s", rpath, strerror(errno));
        setlinebuf(rfp);
        recover();
    }
#ifdef HAVE_SETPROCTITLE
    setproctitle("oul(%lu %lu)", n, k);
#endif
    simple_fact(n, &nf);
    tn = simple_tau(&nf);
    maxfact = 0;
    for (int i = 0; i < nf.count; ++i)
        maxfact += nf.ppow[i].e;
    maxodd = maxfact - nf.ppow[0].e;    /* n is always even */
    init_rootmod(k * maxfact + 1);

    prep_fact();
    prep_maxforce();
    prep_forcep();
    prep_primes();  /* needs forcedp */
    prep_mintau();
    sqg = (uint *)malloc(maxfact * sizeof(uint));

    diag_delay = (debug) ? 0 : DIAG * ticks_per_second;
    log_delay = (debug) ? 0 : LOG * ticks_per_second;
    diagt = diag_delay;
    logt = log_delay;

    init_levels();
    init_value();
    vlevels = (uint *)malloc(forcedp * k * sizeof(uint));
    countr = 0;
    countw = 0;
    countwi = 0;

    wv_o = (mpz_t *)malloc(k * sizeof(mpz_t));
    wv_qq = (mpz_t *)malloc(k * sizeof(mpz_t));
    for (uint i = 0; i < k; ++i) {
        mpz_init(wv_o[i]);
        mpz_init(wv_qq[i]);
    }
    diag_buf = (char *)malloc(DIAG_BUFSIZE);
    init_diag();    /* ignore result: worst case we lose ^Z handling */
}

void report_init(FILE *fp, char *prog) {
    fprintf(fp, "001 %scoul(%u %u)", (start_seen ? "recover " : ""), n, k);
    if (opt_print)
        fprintf(fp, " -o");
    if (minp || maxp) {
        fprintf(fp, " -p");
        if (minp)
            fprintf(fp, "%u:", minp);
        if (maxp)
            fprintf(fp, "%u", maxp);
    }
    if (force_all)
        fprintf(fp, " -f%u", force_all);
    if (gain > 1 || antigain > 1) {
        fprintf(fp, " -g");
        if (antigain > 1)
            fprintf(fp, "%lu:", antigain);
        if (gain > 1)
            fprintf(fp, "%lu", gain);
    }
    if (mpz_sgn(min) || mpz_sgn(max)) {
        fprintf(fp, " -x");
        if (mpz_sgn(min))
            gmp_fprintf(fp, "%Zu:", min);
        if (mpz_sgn(max))
            gmp_fprintf(fp, "%Zu", max);
    }
    if (randseed != 1)
        fprintf(fp, " -s%lu", randseed);
    if (rough)
        fprintf(fp, " -h%u", rough);
    fprintf(fp, "\n");
}

void set_minmax(char *s) {
    char *t = strchr(s, ':');
    if (t) {
        *t = 0;
        if (*s)
            mpz_set_str(min, s, 10);
        else
            mpz_set_ui(min, 0);
        if (t[1])
            mpz_set_str(max, &t[1], 10);
        else
            mpz_set_ui(max, 0);
    } else {
        mpz_set_ui(min, 0);
        mpz_set_str(max, s, 10);
    }
}

void set_gain(char *s) {
    char *t = strchr(s, ':');
    if (t) {
        *t = 0;
        antigain = *s ? strtoul(s, NULL, 10) : 0;
        gain = t[1] ? strtoul(&t[1], NULL, 10) : 0;
    } else {
        antigain = 0;
        gain = strtoul(s, NULL, 10);
    }
}

void set_cap(char *s) {
    char *t = strchr(s, ':');
    if (t) {
        *t = 0;
        minp = *s ? strtoul(s, NULL, 10) : 0;
        maxp = t[1] ? strtoul(&t[1], NULL, 10) : 0;
    } else {
        minp = 0;
        maxp = strtoul(s, NULL, 10);
    }
}

/* Return p if no inverse exists.
 * TODO: mod to ulong, use the fast 64-bit mulmod from MPU-mulmod.h.
 */
ulong small_divmod(mpz_t za, mpz_t zb, ulong p) {
    mpz_set_ui(Z(sdm_p), p);
    mpz_mod_ui(Z(sdm_r), zb, p);
    if (!mpz_invert(Z(sdm_r), Z(sdm_r), Z(sdm_p)))
        return p;
    mpz_mul(Z(sdm_r), Z(sdm_r), za);
    mpz_mod_ui(Z(sdm_r), Z(sdm_r), p);
    return mpz_get_ui(Z(sdm_r));
}

/* Return FALSE if no inverse exists, else sets result = (a / b) mod m.
 */
bool divmod(mpz_t result, mpz_t a, mpz_t b, mpz_t m) {
    mpz_mod(Z(dm_r), b, m);
    if (!mpz_invert(Z(dm_r), Z(dm_r), m))
        return 0;
    mpz_mul(result, Z(dm_r), a);
    mpz_mod(result, result, m);
    return 1;
}

/* This allocation uses what was the next unused prime, so find the
 * index of the new next unused prime.
 */
uint find_nextpi(t_level *cur) {
    uint pi = cur->nextpi;
    while (1) {
        ++pi;
        uint p = sprimes[pi];
        for (uint i = 1; i < level; ++i)
            if (levels[i].p == p)
                goto NEXT_PI;
        return pi;
      NEXT_PI:
        ;
    }
}

bool test_prime(mpz_t qq, mpz_t o, ulong ati) {
    mpz_mul_ui(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return _GMP_is_prob_prime(Z(wv_cand));
}

bool test_zprime(mpz_t qq, mpz_t o, mpz_t ati) {
    mpz_mul(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return _GMP_is_prob_prime(Z(wv_cand));
}

bool test_multi(uint *need, uint nc, ulong ati, uint *t) {
    for (uint i = 0; i < nc; ++i) {
        uint vi = need[i];
        t_tm *tm = &taum[i];
        mpz_mul_ui(tm->n, wv_qq[vi], ati);
        mpz_add(tm->n, tm->n, wv_o[vi]);
        tm->t = t[vi];
        if (!tau_multi_prep(i))
            return 0;
    }
    return tau_multi_run(nc);
}

bool test_zmulti(uint *need, uint nc, mpz_t ati, uint *t) {
    for (uint i = 0; i < nc; ++i) {
        uint vi = need[i];
        t_tm *tm = &taum[i];
        mpz_mul(tm->n, wv_qq[vi], ati);
        mpz_add(tm->n, tm->n, wv_o[vi]);
        tm->t = t[vi];
        if (!tau_multi_prep(i))
            return 0;
    }
    return tau_multi_run(nc);
}

bool test_other(mpz_t qq, mpz_t o, ulong ati, uint t) {
    mpz_mul_ui(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return is_taux(Z(wv_cand), t, 1);
}
bool test_zother(mpz_t qq, mpz_t o, mpz_t ati, uint t) {
    mpz_mul(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return is_taux(Z(wv_cand), t, 1);
}

/* Sort need_prime by qq multiplier ascending */
int prime_comparator(const void *va, const void *vb) {
    uint a = *(uint *)va;
    uint b = *(uint *)vb;
    return mpz_cmp(wv_qq[a], wv_qq[b]);
}

/* Sort need_other by tau (power of 2 ascending, then size ascending),
 * then by qq multiplier ascending
 * TODO: if we implement trial division to higher levels to show
 * a given tau is not reachable, we may want tau by size descending.
 */
uint *oc_t;
int other_comparator(const void *va, const void *vb) {
    uint a = *(uint *)va;
    uint b = *(uint *)vb;
    uint at2 = oc_t[a] ^ (oc_t[a] - 1);
    uint bt2 = oc_t[b] ^ (oc_t[b] - 1);
    if (at2 < bt2)
        return -1;
    if (at2 > bt2)
        return 1;
    if (oc_t[a] < oc_t[b])
        return -1;
    if (oc_t[a] > oc_t[b])
        return 1;
    return mpz_cmp(wv_qq[a], wv_qq[b]);
}

/* Sort inverses by modulus ascending */
int inv_comparator(const void *va, const void *vb) {
    t_mod *a = (t_mod *)va;
    t_mod *b = (t_mod *)vb;
    return (b->m < a->m) - (a->m < b->m);
}

void walk_v(t_level *cur_level, mpz_t start) {
#ifdef SQONLY
    if (!cur_level->have_square)
        return;
#endif
    if (minp && cur_level->maxp <= minp)
        return;

    mpz_t *q[k];
    mpz_t *m = &cur_level->rq;
    uint t[k];
    mpz_t *aq = &cur_level->aq;
    t_mod inv[maxfact * k];
    uint inv_count = 0;
    uint need_prime[k];
    uint need_square[k];
    uint need_other[k];
    uint npc = 0, nqc = 0, noc = 0;

    ++countw;

    mpz_sub(Z(wv_end), max, *m);
    mpz_fdiv_q(Z(wv_end), Z(wv_end), *aq);
    if (mpz_sgn(Z(wv_end)) < 0)
        return;

    if (mpz_sgn(start)) {
        mpz_set(Z(wv_ati), start);
    } else {
        mpz_sub(Z(wv_ati), min, *m);
        mpz_cdiv_q(Z(wv_ati), Z(wv_ati), *aq);
    }

    for (uint vi = 0; vi < k; ++vi) {
        t_value *vp = &value[vi];
        q[vi] = vp->vlevel ? &vp->alloc[vp->vlevel - 1].q : ZP(zone);
        t[vi] = vp->vlevel ? vp->alloc[vp->vlevel - 1].t : n;
        mpz_divexact(wv_qq[vi], *aq, *q[vi]);
        mpz_add_ui(wv_o[vi], *m, vi);
        mpz_divexact(wv_o[vi], wv_o[vi], *q[vi]);
        for (uint ai = 0; ai < vp->vlevel; ++ai) {
            t_allocation *ap = &vp->alloc[ai];
            /* the case for p=2 is handled in advance by update_chinese */
            if (ap->p == 2)
                continue;
            ulong inverse = small_divmod(wv_o[vi], wv_qq[vi], ap->p);
            if (inverse < ap->p) {
                t_mod *ip = &inv[inv_count++];
                ip->v = (inverse) ? ap->p - inverse : 0;
                ip->m = ap->p;
            }
        }
        if (t[vi] == 2)
            need_prime[npc++] = vi;
        else if (t[vi] & 1)
            need_square[nqc++] = vi;
        else
            need_other[noc++] = vi;
    }

#if 0
    qsort(inv, inv_count, sizeof(t_mod), &inv_comparator);
    qsort(need_prime, npc, sizeof(uint), &prime_comparator);
#endif
    oc_t = t;
    qsort(need_other, noc, sizeof(uint), &other_comparator);

    if (nqc) {
        uint sqi = need_square[0];
        mpz_t *oi = &wv_o[sqi];
        mpz_t *qi = q[sqi];
        mpz_t *qqi = &wv_qq[sqi];
        uint ti = t[sqi];
        if (nqc > 1) {
            uint sqj = need_square[1];
            mpz_t *oj = &wv_o[sqj];
            mpz_t *qj = q[sqj];
            mpz_t *qqj = &wv_qq[sqj];
            uint tj = t[sqj];

            mpz_fdiv_q(Z(wv_endr), max, *qi);
            mpz_root(Z(wv_endr), Z(wv_endr), 2);
            /* solve Ax^2 - By^2 = C with x <= D */
            new_pell(*qi, *qj, (int)sqi - (int)sqj, Z(wv_endr));
            uint pc = 0;
            while (next_pell(Z(wv_x), Z(wv_y))) {
                /* v_{sqi} = x^2 . q_i; v_{sqj} = y^2 . q_j */
                mpz_mul(Z(wv_x2), Z(wv_x), Z(wv_x));

                /* verify limit */
                mpz_mul(Z(wv_temp), Z(wv_x2), *qi);
                mpz_sub_ui(Z(wv_temp), Z(wv_temp), sqi);
                if (mpz_cmp(Z(wv_temp), max) > 0)
                    break;  /* CHECKME: should this be an error? */

                ++countwi;
                ++pc;
                if (utime() >= diagt)
                    diag_walk_pell(pc);

                /* verify mod, coprime and tau */
                mpz_fdiv_r(Z(wv_temp), Z(wv_x2), *qqi);
                if (mpz_cmp(Z(wv_temp), *oi) != 0)
                    continue;   /* CHECKME: should this be an error? */
                mpz_mul(Z(wv_y2), Z(wv_y), Z(wv_y));
                mpz_fdiv_r(Z(wv_temp), Z(wv_y2), *qqj);
                if (mpz_cmp(Z(wv_temp), *oj) != 0)
                    continue;   /* CHECKME: should this be an error? */
                mpz_gcd(Z(wv_temp), Z(wv_x), *qqi);
                if (mpz_cmp_ui(Z(wv_temp), 1) != 0)
                    continue;
                mpz_gcd(Z(wv_temp), Z(wv_y), *qqj);
                if (mpz_cmp_ui(Z(wv_temp), 1) != 0)
                    continue;
                /* Note: assume the square roots are small enough to
                 * factorize without fuss */
                if (!is_taux(Z(wv_x), ti, 2))
                    continue;
                if (!is_taux(Z(wv_y), tj, 2))
                    continue;

                mpz_sub(Z(wv_ati), Z(wv_x2), *oi);
                mpz_divexact(Z(wv_ati), Z(wv_ati), *qqi);
                mpz_sub(Z(wv_temp), Z(wv_y2), *oj);
                mpz_divexact(Z(wv_temp), Z(wv_temp), *qqj);
                if (mpz_cmp(Z(wv_ati), Z(wv_temp)) != 0)
                    continue;   /* CHECKME: should this be an error? */

                /* now test the rest */
                for (uint ii = 0; ii < inv_count; ++ii) {
                    t_mod *ip = &inv[ii];
                    if (mpz_fdiv_ui(Z(wv_ati), ip->m) == ip->v)
                        goto next_pell;
                }
                /* note: we may have had more than 2 squares */
                for (uint i = 2; i < nqc; ++i) {
                    uint vi = need_square[i];
                    mpz_mul(Z(wv_cand), wv_qq[vi], Z(wv_ati));
                    mpz_add(Z(wv_cand), Z(wv_cand), wv_o[vi]);
                    if (!is_taux(Z(wv_cand), t[vi], 1))
                        goto next_pell;
                }
                for (uint i = 0; i < npc; ++i) {
                    uint vi = need_prime[i];
                    if (!test_zprime(wv_qq[vi], wv_o[vi], Z(wv_ati)))
                        goto next_pell;
                }
                /* TODO: bail and print somewhere here if 'opt_print' */
#ifdef PARALLEL
                if (!test_zmulti(need_other, noc, Z(wv_ati), t))
                    goto next_pell;
#else
                for (uint i = 0; i < noc; ++i) {
                    uint vi = need_other[i];
                    if (!test_zother(wv_qq[vi], wv_o[vi], Z(wv_ati), t[vi]))
                        goto next_pell;
                }
#endif
                /* have candidate: calculate and apply it */
                mpz_mul(Z(wv_cand), wv_qq[0], Z(wv_ati));
                mpz_add(Z(wv_cand), Z(wv_cand), wv_o[0]);
                mpz_mul(Z(wv_cand), Z(wv_cand), *q[0]);
                candidate(Z(wv_cand));
                break;

              next_pell:
                ;
            }
            /* clear_pell(); */
            return;
        }
        /* gcd(d - 1) for all divisors d of ti */
        uint xi = divisors[ti].gcddm;
        t_results *xr = res_array(cur_level->level);
        qsort(xr->r, xr->count, sizeof(mpz_t), &_mpz_comparator);
        uint rindex = 0;
        if (mpz_sgn(Z(wv_ati)) > 0) {
            mpz_mul(Z(wv_startr), Z(wv_ati), *qqi);
            mpz_add(Z(wv_startr), Z(wv_startr), *oi);
            mpz_root(Z(wv_startr), Z(wv_startr), xi);
            mpz_fdiv_qr(Z(wv_startr), Z(wv_temp), Z(wv_startr), *qqi);
            /* Note: on recover, we expect an exact match here, but on
             * normal entry we don't. */
            for (uint xmi = 0; xmi < xr->count; ++xmi) {
                int cmp = mpz_cmp(xr->r[xmi], Z(wv_temp));
                if (cmp == 0) {
                    rindex = xmi;
                    break;
                } else if (cmp > 0) {
                    if (mpz_sgn(start) == 0) {
                        rindex = xmi;
                        break;
                    }
                    gmp_fprintf(stderr,
                        "from restart %Zu no match found for mod %Zu < %Zu\n",
                        Z(wv_ati), Z(wv_temp), xr->r[xmi]
                    );
                    exit(1);
                }
                if (xmi + 1 == xr->count) {
                    if (mpz_sgn(start) == 0) {
                        rindex = 0;
                        mpz_add_ui(Z(wv_startr), Z(wv_startr), 1);
                        break;
                    }
                    gmp_fprintf(stderr,
                        "from start %Zu no match found for mod %Zu > %Zu\n",
                        Z(wv_ati), Z(wv_temp), xr->r[xmi]
                    );
                    exit(1);
                }
            }
            mpz_mul(Z(wv_qqr), *qqi, Z(wv_startr));
        } else {
            mpz_set_ui(Z(wv_qqr), 0);
        }
        mpz_sub(Z(wv_end), max, *m);
        mpz_fdiv_q(Z(wv_end), Z(wv_end), *aq);
        mpz_add_ui(Z(wv_endr), max, sqi);
        mpz_fdiv_q(Z(wv_endr), Z(wv_endr), *qi);
        mpz_root(Z(wv_endr), Z(wv_endr), xi);

        while (1) {
            mpz_add(Z(wv_r), Z(wv_qqr), xr->r[rindex]);
            if (mpz_cmp(Z(wv_r), Z(wv_endr)) > 0)
                return;
            ++countwi;
            mpz_pow_ui(Z(wv_rx), Z(wv_r), xi);
            mpz_sub(Z(wv_ati), Z(wv_rx), *oi);
            mpz_fdiv_q(Z(wv_ati), Z(wv_ati), *qqi);
            if (utime() >= diagt)
                diag_walk_zv(Z(wv_ati), Z(wv_end));
            for (uint ii = 0; ii < inv_count; ++ii) {
                t_mod *ip = &inv[ii];
                if (mpz_fdiv_ui(Z(wv_ati), ip->m) == ip->v)
                    goto next_sqati;
            }
            if (!is_taux(Z(wv_r), ti, xi))
                goto next_sqati;
            /* note: we have no more squares */
            for (uint i = 0; i < npc; ++i) {
                uint vi = need_prime[i];
                if (!test_zprime(wv_qq[vi], wv_o[vi], Z(wv_ati)))
                    goto next_sqati;
            }
            /* TODO: bail and print somewhere here if 'opt_print' */
#ifdef PARALLEL
            if (!test_zmulti(need_other, noc, Z(wv_ati), t))
                goto next_sqati;
#else
            for (uint i = 0; i < noc; ++i) {
                uint vi = need_other[i];
                if (!test_zother(wv_qq[vi], wv_o[vi], Z(wv_ati), t[vi]))
                    goto next_sqati;
            }
#endif
            /* have candidate: calculate and apply it */
            mpz_mul(Z(wv_cand), wv_qq[0], Z(wv_ati));
            mpz_add(Z(wv_cand), Z(wv_cand), wv_o[0]);
            mpz_mul(Z(wv_cand), Z(wv_cand), *q[0]);
            candidate(Z(wv_cand));
            return;
          next_sqati:
            ++rindex;
            if (rindex >= xr->count) {
                mpz_add(Z(wv_qqr), Z(wv_qqr), *qqi);
                rindex = 0;
            }
        }
    }

    if (!mpz_fits_ulong_p(Z(wv_end)))
        fail("TODO: walk_v.end > 2^64");
    ulong end = mpz_get_ui(Z(wv_end));
    for (ulong ati = mpz_get_ui(Z(wv_ati)); ati <= end; ++ati) {
        ++countwi;
        if (utime() >= diagt)
            diag_walk_v(ati, end);
        for (uint ii = 0; ii < inv_count; ++ii) {
            t_mod *ip = &inv[ii];
            if (ati % ip->m == ip->v)
                goto next_ati;
        }
        /* note: we have no squares */
        for (uint i = 0; i < npc; ++i) {
            uint vi = need_prime[i];
            if (!test_prime(wv_qq[vi], wv_o[vi], ati))
                goto next_ati;
        }
        /* TODO: bail and print somewhere here if 'opt_print' */
#ifdef PARALLEL
        if (!test_multi(need_other, noc, ati, t))
            goto next_ati;
#else
        for (uint i = 0; i < noc; ++i) {
            uint vi = need_other[i];
            if (!test_other(wv_qq[vi], wv_o[vi], ati, t[vi]))
                goto next_ati;
        }
#endif
        /* have candidate: calculate and apply it */
        mpz_mul_ui(Z(wv_cand), wv_qq[0], ati);
        mpz_add(Z(wv_cand), Z(wv_cand), wv_o[0]);
        mpz_mul(Z(wv_cand), Z(wv_cand), *q[0]);
        candidate(Z(wv_cand));
        return;
      next_ati:
        ;
    }
}

void walk_1(t_level *cur_level, uint vi) {
#ifdef SQONLY
    if (!cur_level->have_square)
        return;
#endif
    if (minp && cur_level->maxp <= minp)
        return;

    t_value *vip = &value[vi];
    t_allocation *aip = &vip->alloc[vip->vlevel - 1];
    mpz_sub_ui(Z(w1_v), aip->q, vi);

    if (mpz_cmp(Z(w1_v), min) < 0)
        return;
    ++countw;

    uint t[k];
    uint need_prime[k];
    uint need_square[k];
    uint need_other[k];
    uint npc = 0, nqc = 0, noc = 0;
    for (uint vj = 0; vj < k; ++vj) {
        if (vi == vj)
            continue;
        t_value *vjp = &value[vj];
        t_allocation *ajp = (vjp->vlevel) ? &vjp->alloc[vjp->vlevel - 1] : NULL;
        mpz_add_ui(Z(w1_j), Z(w1_v), vj);
        if (ajp) {
            /* FIXME: replace this with a single initial check of
             * v_0 == rq mod aq, then use divexact */
            mpz_fdiv_qr(Z(w1_j), Z(w1_r), Z(w1_j), ajp->q);
            if (mpz_sgn(Z(w1_r)) != 0)
                return;
            mpz_gcd(Z(w1_r), Z(w1_j), ajp->q);
            if (mpz_cmp(Z(w1_r), Z(zone)) != 0)
                return;
        }
        t[vj] = ajp ? ajp->t : n;
        if (t[vj] == 1) {
            if (mpz_cmp_ui(Z(w1_j), 1) != 0)
                return;
        } else if (t[vj] == 2)
            need_prime[npc++] = vj;
        else if (t[vj] & 1)
            need_square[nqc++] = vj;
        else
            need_other[noc++] = vj;
        mpz_set(wv_o[vj], Z(w1_j));
    }
    ++countwi;
    for (uint i = 0; i < npc; ++i)
        if (!_GMP_is_prob_prime(wv_o[need_prime[i]]))
            return;
    for (uint i = 0; i < nqc; ++i)
        if (!is_taux(wv_o[need_square[i]], t[need_square[i]], 1))
            return;
    oc_t = t;
    qsort(need_other, noc, sizeof(uint), &other_comparator);
    for (uint i = 0; i < noc; ++i)
        if (!is_taux(wv_o[need_other[i]], t[need_other[i]], 1))
            return;
    candidate(Z(w1_v));
    return;
}

/* When some v_j is known to be of the form m.z^g, we keep a running set
 * of possible values of z: * ( (rq + i)/q_j )^{ 1/g } mod aq/q_j.
 * If no such residues exist, this returns FALSE to signal that this
 * case cannot lead to a solution.
 */
bool update_residues(t_level *old, t_level *new,
        uint vi, ulong p, uint x, mpz_t px, bool secondary) {
    if (!old->have_square)
        return 1;

    uint vj = sq0;
    t_value *vjp = &value[vj];
    uint jlevel = vjp->vlevel - 1;
    if (vi == vj) {
        /* Another allocation on our square, so q_j changes but aq/q_j does
         * not. We must divide the known residues by p^((x-1)/g) mod aq/q_j,
         * and if the required power g has changed take roots again. */
        uint oldg = sqg[jlevel - 1];
        uint newg = divisors[vjp->alloc[jlevel].t].gcddm;
        uint divpow = (x - 1) / oldg;

        /* note: old and new give same answer, but old has smaller values */
        mpz_divexact(Z(ur_m), old->aq, vjp->alloc[jlevel - 1].q);
        mpz_ui_pow_ui(Z(ur_ipg), p, divpow);
        if (!mpz_invert(Z(ur_ipg), Z(ur_ipg), Z(ur_m)))
            fail("Cannot find mandatory inverse for %lu^%u", p, divpow);

        t_results *rsrc = res_array(old->level);
        t_results *rdest = res_array(new->level);
        resize_results(rdest, rsrc->count);
        for (uint i = 0; i < rsrc->count; ++i) {
            mpz_mul(rdest->r[i], rsrc->r[i], Z(ur_ipg));
            mpz_mod(rdest->r[i], rdest->r[i], Z(ur_m));
        }
        rdest->count = rsrc->count;

        sqg[jlevel] = newg;
        if (oldg == newg)
            return 1;

        /* we want to upgrade the roots from oldg to newg
         * It is guaranteed that oldg | newg. */
        root_extract(new->level, new->level, newg / oldg, Z(ur_m));
        if (res_array(new->level)->count == 0)
            return 0;
        return 1;
    }
    /* secondary allocation cannot affect roots unless vi == vj */
    if (secondary)
        return 1;

    /* propagate the existing roots */
    uint g = sqg[jlevel];
    mpz_set_si(Z(ur_a), (int)vj - (int)vi);
    if (!divmod(Z(ur_a), Z(ur_a), vjp->alloc[jlevel].q, px))
        return 0;
    mpz_divexact(Z(ur_m), old->aq, vjp->alloc[jlevel].q);
    root_extend(new->level, old->level, Z(ur_m), Z(ur_a), g, p, x - 1, px);
    if (res_array(new->level)->count == 0)
        return 0;
    return 1;
}

void update_chinese(t_level *old, t_level *new, uint vi, mpz_t px) {
    mpz_t zarray[4];
    mpz_t *pxp = PARAM_TO_PTR(px);
    mpz_set_si(Z(uc_minusvi), -(long)vi);

    /* v_0 == -i (mod 2^e) can be upgraded to v_0 = 2^e - i (mod 2^{e + 1}) */
    if (mpz_even_p(px)) {
        mpz_add(Z(uc_minusvi), Z(uc_minusvi), px);
        mpz_mul_2exp(Z(uc_px), px, 1);
        pxp = ZP(uc_px);
    }

    /* TODO: write a custom chinese() */
    memcpy(&zarray[0], old->rq, sizeof(mpz_t));
    memcpy(&zarray[1], Z(uc_minusvi), sizeof(mpz_t));
    memcpy(&zarray[2], old->aq, sizeof(mpz_t));
    memcpy(&zarray[3], *pxp, sizeof(mpz_t));
    if (chinese(new->rq, new->aq, &zarray[0], &zarray[2], 2))
        return;
    fail("chinese failed");
}

/* Record a new square at v_i; return FALSE if invalid.
 * Finds the quadratic (or higher-order) residues; stash them for later
 * propagation if this is the first square; fail if there are none.
 * Note: we just allocated to vi, so at least one allocation exists.
 */
bool alloc_square(t_level *cur, uint vi) {
    t_value *v = &value[vi];
    t_allocation *ap = &v->alloc[v->vlevel - 1];
    uint sqi = cur->have_square++;
    uint g = divisors[ap->t].gcddm;
    if (sqi == 0) {
        sq0 = vi;
        sqg[v->vlevel - 1] = g;
    }

    /* if this is first square, store in the level for further propagation;
     * else use level 0 as temporary */
    uint stash_level = (sqi == 0) ? cur->level : 0;
    /* o = (rq + i) / q */
    mpz_add_ui(Z(asq_o), cur->rq, vi);
    mpz_divexact(Z(asq_o), Z(asq_o), ap->q);
    /* qq = aq / q */
    mpz_divexact(Z(asq_qq), cur->aq, ap->q);

    allrootmod(stash_level, Z(asq_o), g, Z(asq_qq));
    t_results *rp = res_array(stash_level);
    if (rp->count == 0)
        return 0;
    return 1;
}

/* Allocate p^{x-1} to v_{vi}. Returns FALSE if it is invalid, or if no
 * work to do.
 * Updates value[vi]; checks for a new square; calls walk_1() directly
 * if remaining tau == 1.
 */
bool apply_allocv(t_level *prev_level, t_level *cur_level,
        uint vi, ulong p, uint x, mpz_t px, bool secondary) {
    t_value *v = &value[vi];
    t_allocation *prev = (v->vlevel) ? &v->alloc[v->vlevel - 1] : NULL;
    t_allocation *cur = &v->alloc[v->vlevel];
    ++v->vlevel;
    uint prevt = prev ? prev->t : n;
    /* invalid if x does not divide remaining tau */
    if (prevt % x)
        return 0;

    cur->p = p;
    cur->x = x;
    cur->t = prevt / x;
    if (prev)
        mpz_mul(cur->q, prev->q, px);
    else
        mpz_set(cur->q, px);
    /* nothing more to do if we walk_1() */
    if (cur->t == 1) {
        walk_1(cur_level, vi);
        return 0;
    }

    if ((cur->t & 1) && !(prevt & 1))
        if (!alloc_square(cur_level, vi))
            return 0;
    if (prev_level->have_square)
        if (!update_residues(prev_level, cur_level, vi, p, x, px, secondary))
            return 0;
    return 1;
}

/* Allocate p^{x-1} to v_{vi}. Returns FALSE if it is invalid, or if
 * no work to do.
 * Updates level and value, updates or propagates square residues.
 */
bool apply_alloc(t_level *prev, t_level *cur, uint vi, ulong p, uint x) {
    cur->vi = vi;
    cur->p = p;
    cur->x = x;
    cur->have_square = prev->have_square;
    cur->nextpi = prev->nextpi;
    if (p == sprimes[cur->nextpi])
        cur->nextpi = find_nextpi(cur);
    cur->maxp = (p > prev->maxp) ? p : prev->maxp;
    mpz_ui_pow_ui(px, p, x - 1);

    /* Note: the work for this call is wasted if x does not divide t; but
     * that can only happen for forced primes, which is a tiny proportion
     * of the calls (and involves the smallest numbers). Avoiding it means
     * either duplication of other efforts or a more complicated workflow. */
    update_chinese(prev, cur, vi, px);

/* this appears to cost more than it saves in almost all cases */
#ifdef CHECK_OVERFLOW
    /* if rq > max, no solution <= max is possible */
    if (mpz_cmp(cur->rq, max) > 0) {
        /* caller expects vlevel to have been incremented on failure */
        ++value[vi].vlevel;
        return 0;
    }
#endif
    if (!apply_allocv(prev, cur, vi, p, x, px, 0))
        return 0;
    return 1;
}

bool apply_secondary(t_level *prev, t_level *cur, uint vi, ulong p, uint x) {
    mpz_ui_pow_ui(px, p, x - 1);
    return apply_allocv(prev, cur, vi, p, x, px, 1);
}

bool apply_batch(t_level *prev, t_level *cur, t_forcep *fp, uint bi) {
    assert(fp->count > bi);
    t_value *vp;
    cur->is_forced = 1;
    cur->bi = bi;

    t_forcebatch *bp = &fp->batch[bi];
    if (!apply_alloc(prev, cur, bp->vi, fp->p, bp->x))
        return 0;
    /* check if we overshot */
    vp = &value[bp->vi];
    if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
        return 0;

    /* TODO: prep this */
    for (uint i = fp->p; i <= bp->vi; i += fp->p) {
        uint x = simple_valuation(i, fp->p) + 1;
        if (!apply_secondary(prev, cur, bp->vi - i, fp->p, x))
            return 0;
        vp = &value[bp->vi - i];
        if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
            return 0;
    }
    for (uint i = fp->p; bp->vi + i < k; i += fp->p) {
        uint x = simple_valuation(i, fp->p) + 1;
        if (!apply_secondary(prev, cur, bp->vi + i, fp->p, x))
            return 0;
        vp = &value[bp->vi + i];
        if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
            return 0;
    }

    if (opt_alloc && next_prime(fp->p) > maxforce[bp->vi]) {
        if (batch_alloc == opt_batch) {
            /* this is the one batch we want to process */
            ++batch_alloc;
            return 1;
        }
        if (opt_batch < 0)
            disp_batch(cur);
        ++batch_alloc;
        return 0;
    }
    return 1;
}

/* find the best entry to progress: the one with the highest tau()
 * still to fulfil, or (on equality) with the highest q, but having
 * at least one factor to allocate.
 * If there is no best entry, returns k.
 */
uint best_v(void) {
    uint vi, ti = 0;
    mpz_t *qi;
    for (uint vj = 0; vj < k; ++vj) {
        t_value *vpj = &value[vj];
        t_allocation *apj = (vpj->vlevel) ? &vpj->alloc[vpj->vlevel - 1] : NULL;
        uint tj = apj ? apj->t : n;
        mpz_t *qj = apj ? &apj->q : ZP(zone);

        /* skip if no odd prime factor */
        if (divisors[tj].high <= 2)
            continue;
        /* skip prime powers when capped */
        if (maxp && (tj & 1) && divisors[tj].alldiv == 2)
            continue;
        if (ti) {
            /* skip if not higher tau, or same tau with higher q */
            if (tj < ti)
                continue;
            if (tj == ti && mpz_cmp(*qj, *qi) <= 0)
                continue;
        }
        vi = vj;
        ti = tj;
        qi = qj;
    }
    return ti ? vi : k;
}

/* Calculate the minimum contribution from primes satisfying the given tau.
 */
void mintau(mpz_t mint, uint vi, uint t) {
    /* quick version: given the minimum prime p that can be used, we
     * calculate p^k where k = sum{p_i - 1} over the primes dividing
     * t _with multiplicity_.
     */
    uint minnext = sprimes[levels[level - 1].nextpi];
    mpz_ui_pow_ui(mint, minnext, divisors[t].sumpm);
}

/* return the maximum prime to iterate to */
ulong limit_p(uint vi, uint x, uint nextt) {
    t_value *vp = &value[vi];
    t_allocation *ap = (vp->vlevel) ? &vp->alloc[vp->vlevel - 1] : NULL;
    mpz_add_ui(Z(lp_x), max, vi);
    if (ap)
        mpz_div(Z(lp_x), Z(lp_x), ap->q);
    /* else notional ap->q is 1 */

    if (x == nextt && divisors[x].high == x) {
        /* We are allocating p^{x-1} with x prime, leaving x as the
         * remaining tau; so this and the remaining allocation will
         * be of the form p^{x-1}.q^{x-1}, and we can set the limit
         * to max^{1/2(x-1)}.
         */
        mpz_root(Z(lp_x), Z(lp_x), 2 * (x - 1));
    } else {
        /* divide through by the minimum contribution that could supply the
         * remaining tau */
        if (nextt > 1) {
            mintau(Z(lp_mint), vi, nextt);
            mpz_div(Z(lp_x), Z(lp_x), Z(lp_mint));
        }
        mpz_root(Z(lp_x), Z(lp_x), x - 1);
    }

    if (maxp && mpz_cmp_ui(Z(lp_x), maxp) > 0)
        return maxp;
    if (mpz_fits_ulong_p(Z(lp_x)))
        return mpz_get_ui(Z(lp_x));
    return 0;
}

/*
 * return 0 if nothing more to do at this level for any x;
 * return 1 if nothing more to do for this x;
 * return 2 if prepped for this x with work to do.
 */
uint prep_unforced_x(t_level *prev, t_level *cur, ulong p) {
    uint ti = cur->ti;
    uint x = divisors[ti].div[cur->di];
    uint vi = cur->vi;
    t_value *vp = &value[vi];
    t_allocation *ap = (vp->vlevel) ? &vp->alloc[vp->vlevel - 1] : NULL;

    /* pick up any previous unforced x */
    uint nextt = ti / x;
    if (p == 0) {
        uint prevx = (ap && ap->p > maxforce[vi]) ? ap->x : 0;
        if (x == prevx)
            p = ap->p;      /* skip smaller p, we already did the reverse */
        else if (x <= prevx && divisors[x].high == divisors[prevx].high)
            return 1;   /* skip this x, we already did the reverse */
        else if (x > nextt && divisors[x].high == divisors[nextt].high)
            /* skip this x, we already did any possible continuation in
             * reverse. */
            return 1;
        else
            p = maxforce[vi];
    } /* else we're continuing from known p */

    /* try p^{x-1} for all p until q_i . p^{x-1} . minrest > max + i */
    ulong limp = limit_p(vi, x, nextt);
    if (limp == 0) {
        if (!prev->have_square) {
            diag_plain();
            keep_diag();
            report("002 %s: outsize limit %Zu without square\n",
                    diag_buf, Z(lp_x));
            return 1;   /* skip this cycle */
        }
        /* force walk */
#ifdef SQONLY
        if (prev->have_square)
            walk_v(prev, Z(zero));
#else
        walk_v(prev, Z(zero));
#endif
        return 0;
    } else if (limp < p + 1)
        return 1;   /* nothing to do here */
    mpz_add_ui(Z(r_walk), max, vi);
    mpz_fdiv_q(Z(r_walk), Z(r_walk), prev->aq);
    if (prev->have_square) {
        if (prev->have_square == 1) {
            /* if we fix a square, expect to actually walk only the g'th
             * roots of rq mod aq */
            uint g = sqg[value[sq0].vlevel - 1];
            mpz_root(Z(r_walk), Z(r_walk), g);
            mpz_mul_ui(Z(r_walk), Z(r_walk), res_array(prev->level)->count);
        } else {
            /* if we fix multiple squares, we'll solve a Pell equation;
             * treat that as effectively free */
            mpz_set_ui(Z(r_walk), 0);
        }
    }
    if (gain > 1)
        mpz_mul_ui(Z(r_walk), Z(r_walk), gain);
    if (antigain > 1)
        mpz_fdiv_q_ui(Z(r_walk), Z(r_walk), antigain);
    if (mpz_fits_ulong_p(Z(r_walk))
        && mpz_get_ui(Z(r_walk)) < limp - p
    ) {
#ifdef SQONLY
        if (prev->have_square)
            walk_v(prev, Z(zero));
        else if (level > 1 && !prev->is_forced)
            prev->p = prev->limp;
#else
        walk_v(prev, Z(zero));
#endif
        return 0;
    }
    cur->p = p;
    cur->x = x;
    cur->limp = limp;
    cur->max_at = seen_best;
    /* TODO: do some constant alloc stuff in advance */
    /* TODO: special case for nextt == 1 */
    return 2;
}

/* On recovery, set up the recursion stack to the point we had reached.
 */
void insert_stack(void) {
    /* first insert forced primes */
    for (uint fpi = 0; fpi < forcedp; ++fpi) {
        t_forcep *fp = &forcep[fpi];
        uint p = fp->p;
        uint maxx = 0, mini;
        for (uint vi = 0; vi < k; ++vi) {
            t_fact *rs = &rstack[vi];
            if (rs->count && rs->ppow[rs->count - 1].p == p) {
                uint x = rs->ppow[rs->count - 1].e + 1;
                if (maxx < x) {
                    maxx = x;
                    mini = vi;
                }
            }
        }
        /* find the batch */
        if (maxx == 0) {
            if (fp->batch[fp->count -1].x != 0)
                goto insert_check;
            /* this prime unforced, so any remaining ones are too */
            break;
        }

        uint bi;
        for (bi = 0; bi < fp->count; ++bi) {
            t_forcebatch *b = &fp->batch[bi];
            if (b->x == maxx && b->vi == mini)
                break;
        }
        if (bi >= fp->count) {
            if (fp->batch[fp->count -1].x != 0)
                fail("no batch found for %u^{%u-1} at v_%u", p, maxx, mini);
            /* this prime unforced, so any remaining ones are too */
            break;
        }

        STOREVL(vl_forced++);
        --rstack[mini].count;
        t_level *prev = &levels[level - 1];
        t_level *cur = &levels[level];
        cur->is_forced = 1;
        cur->bi = bi;

        /* CHECKME: this returns 0 if t=1 */
        if (!apply_alloc(prev, cur, mini, p, maxx))
            fail("could not apply_alloc(%u, %lu, %u)", mini, p, maxx);
        /* TODO: prep this, per apply_batch */
        for (uint j = p; j <= mini; j += p) {
            uint vj = mini - j;
            uint x = simple_valuation(j, p) + 1;
            t_fact *rs = &rstack[vj];
            t_ppow *rsp = rs->count ? &rs->ppow[rs->count - 1] : NULL;
            if (!rsp || rsp->p != p || rsp->e + 1 != x)
                fail("missing secondary %u^%u at %u", p, x, vj);
            --rs->count;

            if (!apply_secondary(prev, cur, vj, p, x))
                fail("could not apply_secondary(%u, %lu, %u)", vj, p, x);
        }
        for (uint j = p; mini + j < k; j += p) {
            uint vj = mini + j;
            uint x = simple_valuation(j, p) + 1;
            t_fact *rs = &rstack[vj];
            t_ppow *rsp = rs->count ? &rs->ppow[rs->count - 1] : NULL;
            if (!rsp || rsp->p != p || rsp->e + 1 != x)
                fail("missing secondary %u^%u at %u", p, x, vj);
            --rs->count;

            if (!apply_secondary(prev, cur, vj, p, x))
                fail("could not apply_secondary(%u, %lu, %u)", vj, p, x);
        }
        ++level;
    }
    /* now insert the rest */
    while (1) {
        uint vi = best_v();
        if (vi == k)
            break;
        t_fact *rs = &rstack[vi];
        if (rs->count == 0)
            break;
        --rs->count;
        ulong p = rs->ppow[rs->count].p;
        uint x = rs->ppow[rs->count].e + 1;

        t_level *prev = &levels[level - 1];
        t_level *cur = &levels[level];
        t_value *vp = &value[vi];

        uint ti = (vp->vlevel) ? vp->alloc[vp->vlevel - 1].t : n;
        t_divisors *dp = &divisors[ti];
        if (dp->highdiv == 0)
            fail("best_v() returned %u, but nothing to do there", vi);

        uint di;
        for (di = 0; di <= dp->highdiv; ++di) {
            if (di == dp->highdiv)
                fail("x=%u is not a highdiv of t=%u\n", x, ti);
            if (dp->div[di] == x)
                break;
        }
        cur->vi = vi;
        cur->ti = ti;
        cur->di = di;

        /* note: must pass in p=0 for it to calculate limp */
        uint pux = prep_unforced_x(prev, cur, 0);
        if (pux != 2)
            fail("prep_nextt %u for %lu^%u at %u\n",
                    pux, p, x, vi);

        /* CHECKME: this returns 0 if t=1 */
        if (!apply_alloc(prev, cur, vi, p, x))
            fail("could not apply_alloc(%u, %lu, %u)", vi, p, x);
        ++level;
    }
  insert_check:
    /* check we found them all */
    for (uint vi = 0; vi < k; ++vi) {
        t_fact *rs = &rstack[vi];
        if (rs->count) {
            t_ppow pp = rs->ppow[rs->count - 1];
            fail("failed to inject %lu^%u at v_%u", pp.p, pp.e, vi);
        }
    }
}

/* we emulate recursive calls via the levels[] array */
void recurse(void) {
    ulong p;
    uint x;
    t_level *prev_level, *cur_level;

    if (have_rwalk) {
        prev_level = &levels[level - 1];
        walk_v(prev_level, rwalk_from);
        goto derecurse;
    }

    while (1) {
        ++countr;
        prev_level = &levels[level - 1];
        cur_level = &levels[level];

        /* recurse deeper */
        {
            uint fi = level - 1;
            if (fi < forcedp && (fi == 0 || prev_level->is_forced)) {
                t_forcep *fp = &forcep[fi];
                if (fp->count == 0)
                    goto unforced;
                STOREVL(vl_forced++);
                if (!apply_batch(prev_level, cur_level, fp, 0)) {
                    FETCHVL(vl_forced - 1);
                    goto continue_recurse;
                }
                ++level;
                if (utime() >= diagt)
                    diag_plain();
                continue;   /* deeper */
            }
        }
      unforced:
        {
            /* note: cur_level->is_forced is either always 0 by calloc(),
             * or gets set to 0 on the tail of a batch */
            /* cur_level->is_forced = 0; */
            uint vi = best_v();
            /* TODO: walk_v() directly at previous level, if best_v() would
             * give same result each time.
             */
            if (vi == k) {
#ifdef SQONLY
                if (prev_level->have_square)
                    walk_v(prev_level, Z(zero));
                else if (level > 1 && !prev_level->is_forced)
                    prev_level->p = prev_level->limp;
#else
                walk_v(prev_level, Z(zero));
#endif
                goto derecurse;
            }
            t_value *vp = &value[vi];
            uint ti = (vp->vlevel) ? vp->alloc[vp->vlevel - 1].t : n;
            t_divisors *dp = &divisors[ti];
            if (dp->highdiv == 0)
                fail("best_v() returned %u, but nothing to do there", vi);
            cur_level->vi = vi;
            cur_level->ti = ti;
            cur_level->di = 0;
            goto have_unforced_x;
        }
      continue_unforced_x:
        ++cur_level->di;
      have_unforced_x:
        {
            if (cur_level->di >= divisors[cur_level->ti].highdiv)
                goto derecurse;
            switch (prep_unforced_x(prev_level, cur_level, 0)) {
                case 0:
                    /* nothing to do for any x */
                    goto derecurse;
                case 1:
                    /* nothing to do for this x */
                    goto continue_unforced_x;
                case 2:
                    /* ok, continue for this x */
                    ;
            }
            goto continue_unforced;
        }
        break;
      /* entry point, must set prev_level/cur_level before using */
      derecurse:
        --level;
        if (level == 0)
            break;
        prev_level = &levels[level - 1];
        cur_level = &levels[level];
        if (cur_level->is_forced) {
            /* unapply the batch */
            FETCHVL(vl_forced - 1);
        } else
            --value[cur_level->vi].vlevel;
        /* goto continue_recurse; */
      continue_recurse:
        if (cur_level->is_forced) {
            uint fi = level - 1;
            t_forcep *fp = &forcep[fi];

            uint bi = cur_level->bi + 1;
            if (bi >= fp->count) {
                --vl_forced;
                goto derecurse;
            }
            if (fp->batch[bi].x == 0) {
                cur_level->is_forced = 0;
                FETCHVL(--vl_forced);
                if (opt_alloc) {
                    if (batch_alloc == opt_batch) {
                        /* this is the one batch we want to process */
                        ++batch_alloc;
                        goto unforced;
                    }
                    if (opt_batch < 0)
                        disp_batch(prev_level);
                    ++batch_alloc;
                    goto derecurse;
                }
                goto unforced;
            }
            if (!apply_batch(prev_level, cur_level, fp, bi)) {
                /* unapply a possible partial batch */
                FETCHVL(vl_forced - 1);
                goto continue_recurse;
            }
            ++level;
            if (utime() >= diagt)
                diag_plain();
            continue;
        }
      continue_unforced:
        {
            ulong p = cur_level->p;
            /* recalculate limit if we have an improved maximum */
            if (seen_best > cur_level->max_at)
                if (!prep_unforced_x(prev_level, cur_level, p))
                    goto continue_unforced_x;
            /* note: only valid to use from just below here */
          redo_unforced:
            p = next_prime(p);
            if (p > cur_level->limp)
                goto continue_unforced_x;
            /* TODO: save max_used_p, use it to short-circuit */
            for (uint li = 1; li < level; ++li)
                if (p == levels[li].p)
                    goto redo_unforced;
            /* note: this returns 0 if t=1 */
            if (!apply_alloc(prev_level, cur_level, cur_level->vi, p, cur_level->x)) {
                if (utime() >= diagt)
                    diag_plain();
                --value[cur_level->vi].vlevel;
                /* not redo_unforced, we may have improved max */
                goto continue_unforced;
            }
            ++level;
            if (utime() >= diagt)
                diag_plain();
            continue;   /* deeper */
        }
    }
}

int main(int argc, char **argv, char **envp) {
    int i = 1;
#ifdef HAVE_SETPROCTITLE
    setproctitle_init(argc, argv, envp);
#endif
    init_pre();
    while (i < argc && argv[i][0] == '-') {
        char *arg = argv[i++];
        if (strncmp("--", arg, 2) == 0)
            break;
        if (arg[1] == 'x')
            set_minmax(&arg[2]);
        else if (arg[1] == 'g')
            set_gain(&arg[2]);
        else if (arg[1] == 'p')
            set_cap(&arg[2]);
        else if (arg[1] == 'r') {
            rpath = (char *)malloc(strlen(&arg[2]) + 1);
            strcpy(rpath, &arg[2]);
        } else if (arg[1] == 'f')
            force_all = strtoul(&arg[2], NULL, 10);
        else if (arg[1] == 's')
            randseed = strtoul(&arg[2], NULL, 10);
        else if (arg[1] == 'h')
            rough = strtoul(&arg[2], NULL, 10);
        else if (strncmp("-a", arg, 2) == 0)
            opt_alloc = 1;
        else if (arg[1] == 'b') {
            opt_batch = strtoul(&arg[2], NULL, 10);
            opt_alloc = 1;
        } else if (strncmp("-o", arg, 2) == 0)
            opt_print = 1;
        else if (strncmp("-d", arg, 2) == 0)
            ++debug;
        else
            fail("unknown option '%s'", arg);
    }
    if (i + 2 == argc) {
        n = strtoul(argv[i++], NULL, 10);
        k = strtoul(argv[i++], NULL, 10);
    } else
        fail("wrong number of arguments");
    if (force_all > k)
        fail("require force_all <= k");

    init_post();
    report_init(stdout, argv[0]);
    if (rfp) report_init(rfp, argv[0]);
#if 0
    char s[] = "7^2 2.71^2 3^8 2^2.5^2 11^2.29^2 (0.00s)\n";
    parse_305(s);
#endif
    if (rstack)
        insert_stack();
    recurse();
    keep_diag();

    clock_t tz = utime();
    report("367 coul(%u, %u): recurse %lu, walk %lu, walkc %lu (%.2fs)\n",
            n, k, countr, countw, countwi, seconds(tz));
    if (seen_best)
        report("200 f(%u, %u) = %Zu (%.2fs)\n", n, k, max, seconds(tz));
    done();
    return 0;
}
