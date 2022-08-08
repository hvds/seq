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
#include "coulfact.h"
#include "diag.h"

/* from MPUG */
#include "factor.h"
#include "gmp_main.h"
#include "utility.h"
#include "primality.h"

/* primary parameters - we are searching for D(n, k), the least d such
 * that tau(d + i) = n for all 0 <= i < k.
 */
uint n, k;

typedef enum {
    zero, zone,                 /* constants */
    res_m, res_e, res_px,       /* check_residue */
    sqm_t, sqm_q, sqm_b, sqm_z, sqm_x,  /* sqrtmod_t */
    uc_minusvi,                 /* update_chinese */
    wv_ati, wv_end, wv_cand,    /* walk_v */
    lp_x, lp_mint,              /* limit_p */
    r_walk,                     /* recurse */

    sdm_p, sdm_r,               /* small_divmod (TODO) */
    np_p,                       /* next_prime (TODO) */

    MAX_ZSTASH
} t_zstash;
mpz_t *zstash;
static inline mpz_t *ZP(t_zstash e) { return &zstash[e]; }
#define Z(e) *ZP(e)

typedef unsigned char bool;

typedef struct s_mod {
    ulong v;
    ulong m;
} t_mod;

/* 'divisors[i].div' is a list of divisors of 'i' in descending order of
 * highest prime factor, then ascending. 'high' is the highest prime
 * factor of 'i'; 'alldiv' is the number of factors; 'highdiv' is the
 * number of factors that are a multiple of 'high'.
 */
typedef struct s_divisors {
    uint alldiv;
    uint highdiv;
    uint high;
    uint *div;
} t_divisors;
t_divisors *divisors = NULL;

/* For prime p < k, we "force" allocation of powers in a batch to ensure
 * that multiple allocations of the same prime are coherent. Whereas normal
 * allocation considers only p^{x-1} where x is divisible by the highest
 * prime dividing v_i.t, forced primes may allocate any p^{x-1} with
 * x | v_i.t.
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
uint forcedp;
t_forcep *forcep = NULL;

/* When allocation forces the residue of some v_i to be square, we want
 * to capture some information, and check if that is consistent with
 * existing allocations.
 */
#define MAX_SQUARE 2
typedef struct s_square {
    uint vi;
    mpz_t m;    /* v_{vi} = mz^2 for some z */
} t_square;
t_square *squares = NULL;
t_square *sqspare = NULL;

/* Each level of "recursion" allocates one prime power p^{x-1} with x | n
 * to one value v_i. It may be "forced", in which case it is part of a
 * batch of simultaneous allocations of the same p to different i (and
 * in which case derecursion should unwind the whole batch), or not, in
 * which case no other v_j will be divisible by p.
 */
typedef struct s_level {
    bool is_forced;
    uint vi;    /* allocation of p^x into v_i */
    ulong p;
    /* union */
        uint bi;    /* batch index, if forced */
    /* .. with */
        uint di;    /* divisor index, if unforced */
        uint ti;    /* target tau for v_i */
        ulong limp; /* limit for p */
        uint max_at;/* which max value used for limp calculation */
    /* end union */
    uint x;
    mpz_t aq;   /* running LCM of allocated p^x */
    mpz_t rq;   /* running CRT of (-i) % p^x */
    uint have_square;
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

/* allocations in each value before applying nth forced prime */
uint *vlevels = NULL;
uint vl_forced = 0;
static inline uint *VLP(uint level) { return &vlevels[level * k]; }
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

long ticks_per_second;
clock_t ticks = 0;

mpz_t min, max;
uint seen_best = 0;
ulong gain = 1;
/* TODO: implement minp */
uint minp = 0, maxp = 0;
uint runid = 0;
bool opt_print = 0;
uint force_all = 0;
bool debug = 0;

char *rpath = NULL;
FILE *rfp = NULL;
bool start_seen = 0;
t_fact *rstack = NULL;
bool have_rwalk = 0;
mpz_t rwalk_from;
mpz_t rwalk_to;

t_fact nf;
uint tn;
uint maxfact;
uint *maxforce = NULL;
mpz_t px;   /* p^x */
mpz_t *wv_o = NULL, *wv_qq = NULL;  /* wv_o[k], wv_qq[k] */

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
    clock_t t1 = times(NULL);

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
    clock_t t1 = times(NULL);

    prep_show_v();  /* into diag_buf */
    if (!(debug && ati))
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

void free_levels(void) {
    for (uint i = 0; i < k * maxfact + 1; ++i) {
        t_level *l = &levels[i];
        mpz_clear(l->aq);
        mpz_clear(l->rq);
    }
    free(levels);
}

void init_levels(void) {
    levels = (t_level *)calloc(k * maxfact + 1, sizeof(t_level));
    for (uint i = 0; i < k * maxfact + 1; ++i) {
        t_level *l = &levels[i];
        mpz_init(l->aq);
        mpz_init(l->rq);
    }
    mpz_set_ui(levels[0].aq, 1);
    mpz_set_ui(levels[0].rq, 0);
    levels[0].have_square = 0;
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
    free(vlevels);
    free_value();
    free_levels();
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
    for (uint i = 0; i < MAX_SQUARE + 1; ++i)
        mpz_clear(squares[i].m);
    free(squares);
    mpz_clear(px);
    free_fact(&nf);
    mpz_clear(max);
    mpz_clear(min);
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
    ticks_per_second = sysconf(_SC_CLK_TCK);
    ticks = times(NULL);
    mpz_init_set_ui(min, 0);
    mpz_init_set_ui(max, 0);
    init_fact(&nf);
    mpz_init(px);
    squares = (t_square *)malloc((MAX_SQUARE + 1) * sizeof(t_square));
    sqspare = &squares[MAX_SQUARE];
    for (uint i = 0; i < MAX_SQUARE + 1; ++i)
        mpz_init(squares[i].m);
    zstash = (mpz_t *)malloc(MAX_ZSTASH * sizeof(mpz_t));
    for (t_zstash i = 0; i < MAX_ZSTASH; ++i)
        mpz_init(Z(i));
    mpz_set_ui(Z(zero), 0);
    mpz_set_ui(Z(zone), 1);
}

void recover(void) {
    char *last305 = NULL;
    char *last202 = NULL;
    char *curbuf = NULL;
    size_t len = 120;

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
            ftruncate(fileno(rfp), offset - nread);
            if (freopen(NULL, "a+", rfp) == NULL)
                fail("could not reopen %s: %s", rpath, strerror(errno));
            break;
        }
        if (strncmp("305 ", curbuf, 4) == 0) {
            char *t = last305;
            last305 = curbuf;
            curbuf = t;
        } else if (strncmp("202 ", curbuf, 4) == 0) {
            char *t = last202;
            last202 = curbuf;
            curbuf = t;
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
        mpz_t candidate;
        if (EOF == sscanf(last202, "202 Candidate %n%*[0-9]%n (%*[0-9.]s)\n",
                &start, &end))
            fail("error parsing 202 line '%s'", last202);
        last202[end] = 0;
        mpz_init_set_str(candidate, &last202[start], 10);
        if (mpz_sgn(max) == 0 || mpz_cmp(max, candidate) > 0)
            mpz_set(max, candidate);
        mpz_clear(candidate);
        free(last202);
    }
    if (last305) {
        char *s = last305 + 3;
        double dtime;
        t_ppow pp;

        rstack = (t_fact *)malloc(k * sizeof(t_fact));
        for (int i = 0; i < k; ++i)
            init_fact(&rstack[i]);

        for (int i = 0; i < k; ++i) {
            assert(s[0] == ' ');
            ++s;
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
    free(curbuf);
    free(last305);
    free(last202);
}

int cmp_high(const void *va, const void *vb) {
    uint a = *(uint *)va, b = *(uint *)vb;
    return divisors[b].high - divisors[a].high;
}

/* recurse() wants the list of powers to try: each divisor of t_i (which
 * itself divides n) that is divisible by the highest odd prime factor
 * dividing t_i, in increasing order.
 * prep_forcep() wants the full list of divisors, but in similar order.
 * For each power, recurse() also wants to know which powers to skip
 * if the previous power was a given value, but that's simply:
 * skip x' if x' < x and high(x') == high(x)
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
        return;
    }
    for (int i = 0; i < k; ++i) {
        int mf = k - i - 1;
        if (i > mf)
            mf = i;
        if (force_all > mf)
            mf = force_all;
        maxforce[i] = mf;
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
    if (runid) {
        char buf[100];
        snprintf(buf, sizeof(buf), "%s/%u.%u-%u",
                "logs/o", n, k, runid);
        rpath = (char *)malloc(strlen(buf) + 1);
        strcpy(rpath, buf);
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

    prep_fact();
    prep_maxforce();
    prep_forcep();
    prep_mintau();

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
    fprintf(fp, "001 %s%s", (start_seen ? "recover " : ""), prog);
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
    if (gain > 1)
        fprintf(fp, " -g%lu", gain);
    if (mpz_sgn(min) || mpz_sgn(max)) {
        fprintf(fp, " -x");
        if (mpz_sgn(min))
            gmp_fprintf(fp, "%Zu:", min);
        if (mpz_sgn(max))
            gmp_fprintf(fp, "%Zu", max);
    }
    fprintf(fp, " %u %u\n", n, k);
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

/* return TRUE if a is a quadratic residue mod m
 * TODO: since we don't need the root, can we speed this up?
 */
bool has_sqrtmod(mpz_t a, mpz_t m) {
    return sqrtmod_t(Z(sqm_x), a, m, Z(sqm_t), Z(sqm_q), Z(sqm_b), Z(sqm_z));
}

/* Return TRUE if allocating p^{x-1} at v_i is consistent with the known
 * square v_j = mz^2.
 * TODO: we probably always have mpz_fits_ulong_p(m)
 * TODO: for p > 2 see modinvert() in MPU-util.c
 * TODO: p == 2 seems like it should be easy, but existing algorithm
 * doesn't special-case it.
 */
bool check_residue(uint vi, ulong p, uint x, uint vj, mpz_t m) {
    int d = vj - vi;
    if (d == 0)
        return 1;
    uint e1 = x - 1;
    uint e2 = 0;
    while ((d % p) == 0) {
        ++e2;
        d /= p;
    }

    mpz_set(Z(res_m), m);
    t_value *vpj = &value[vj];
    for (uint i = 0; i < vpj->vlevel; ++i)
        if (vpj->alloc[i].p == p) {
            if (e2 == 0)
                return 0;
            if (e1 == e2)
                return 1;   /* see below */
            if (mpz_divisible_ui_p(Z(res_m), p)) {
                mpz_divexact_ui(Z(res_m), Z(res_m), p);
                --e1;
                --e2;
            }
            break;
        }

    if (e1 == e2) {
        /* we know only that valuation(v_j, p) > e1, punt on this */
        return 1;
    }

    /* we have v_j = p^{min(e1, e2)} . m . z^2 */
    if (p == 2)
        mpz_mul_2exp(Z(res_px), Z(zone), (e1 < e2) ? e2 - e1 : e1 - e2);
    else
        mpz_set_ui(Z(res_px), p);
    /* return true iff d / m (mod px) is a quadratic residue (mod px) */
    if (!mpz_invert(Z(res_m), Z(res_m), Z(res_px)))
        fail("logic error, p ~| m, so m should be invertible");
    if (d != 1) {
        mpz_mul_si(Z(res_m), Z(res_m), d);
        mpz_mod(Z(res_m), Z(res_m), Z(res_px));
    }
    return has_sqrtmod(Z(res_m), Z(res_px));
}

void update_chinese(t_level *old, t_level *new, uint vi, mpz_t px) {
    mpz_t zarray[4];
    mpz_set_si(Z(uc_minusvi), -(long)vi);

    memcpy(&zarray[0], old->rq, sizeof(mpz_t));
    memcpy(&zarray[1], Z(uc_minusvi), sizeof(mpz_t));
    memcpy(&zarray[2], old->aq, sizeof(mpz_t));
    memcpy(&zarray[3], px, sizeof(mpz_t));
    if (chinese(new->rq, new->aq, &zarray[0], &zarray[2], 2))
        return;
    fail("chinese failed");
}

/* Record a new square at v_i; return FALSE if any v_j factor is not a
 * residue.
 */
bool alloc_square(uint vi) {
    t_value *v = &value[vi];
    t_level *lp = &levels[level - 1];
    uint sqi = lp->have_square++;
    t_square *sqp = (sqi < MAX_SQUARE) ? &squares[sqi] : sqspare;
    sqp->vi = vi;
    mpz_set_ui(sqp->m, 1);
    for (uint ai = 0; ai < v->vlevel; ++ai) {
        t_allocation *ap = &v->alloc[ai];
        if (ap->x & 1)
            continue;
        mpz_mul_ui(sqp->m, sqp->m, ap->p);
    }
    for (uint vj = 0; vj < k; ++vj) {
        if (vi == vj)
            continue;
        v = &value[vj];
        for (uint ai = 0; ai < v->vlevel; ++ai) {
            t_allocation *ap = &v->alloc[ai];
            if (!check_residue(vj, ap->p, ap->x, vi, sqp->m))
                return 0;
        }
    }
    return 1;
}

/* returns TRUE if this newly creates a square */
bool apply_allocv(uint vi, ulong p, uint x, mpz_t px) {
    t_value *v = &value[vi];
    t_allocation *prev = (v->vlevel) ? &v->alloc[v->vlevel - 1] : NULL;
    t_allocation *cur = &v->alloc[v->vlevel];
    ++v->vlevel;
    uint prevt = prev ? prev->t : n;
    if (prevt % x)
        return 0;

    cur->p = p;
    cur->x = x;
    cur->t = prevt / x;
    if (prev)
        mpz_mul(cur->q, prev->q, px);
    else
        mpz_set(cur->q, px);

    if ((cur->t & 1) && !(prevt & 1))
        if (!alloc_square(vi))
            return 0;
    return 1;
}

/* Allocate p^{x-1} to v_{vi}. Returns FALSE if it is invalid. */
bool apply_alloc(uint vi, ulong p, uint x) {
    assert(level > 0);
    t_level *prev = &levels[level - 1];
    t_level *cur = &levels[level];
    ++level;
    cur->vi = vi;
    cur->p = p;
    cur->x = x;
    cur->have_square = prev->have_square;
    mpz_set_ui(px, p);
    mpz_pow_ui(px, px, x - 1);
    update_chinese(prev, cur, vi, px);
    return apply_allocv(vi, p, x, px);
}

bool apply_secondary(t_level *cur, uint vi, ulong p, uint x) {
    mpz_set_ui(px, p);
    mpz_pow_ui(px, px, x - 1);
    return apply_allocv(vi, p, x, px);
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
        if (maxx == 0)
            fail("no forced prime %u found", p);
        if (!apply_alloc(mini, p, maxx))
            fail("could not apply_alloc(%u, %lu, %u)", mini, p, maxx);
        t_level *lp = &levels[level];
        for (uint vi = 0; vi < k; ++vi) {
            t_fact *rs = &rstack[vi];
            if (rs->count && rs->ppow[rs->count - 1].p == p) {
                --rs->count;
                if (vi == mini)
                    continue;
                uint x = rs->ppow[rs->count].e + 1;
                if (!apply_secondary(lp, vi, p, x))
                    fail("could not apply_secondary(%u, %lu, %u)", vi, p, x);
            }
        }

        uint bi;
        for (bi = 0; bi < fp->count; ++bi) {
            t_forcebatch *b = &fp->batch[bi];
            if (b->x == maxx && b->vi == mini)
                break;
        }
        if (bi >= fp->count)
            fail("no batch found for %u^{%u-1} at v_%u", p, maxx, mini);
        levels[level - 1].is_forced = 1;
        levels[level - 1].bi = bi;
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
        uint p = rs->ppow[rs->count].p;
        uint x = rs->ppow[rs->count].e + 1;
        if (!apply_alloc(vi, p, x))
            fail("could not apply_alloc(%u, %lu, %u)", vi, p, x);
    }
    /* check we found them all */
    for (uint vi = 0; vi < k; ++vi) {
        t_fact *rs = &rstack[vi];
        if (rs->count) {
            t_ppow pp = rs->ppow[rs->count - 1];
            fail("failed to inject %lu^%u at v_%u", pp.p, pp.e, vi);
        }
    }
}

/* return the maximum prime to iterate to */
ulong limit_p(uint vi, uint x, uint nextt) {
    t_value *vp = &value[vi];
    t_allocation *ap = (vp->vlevel) ? &vp->alloc[vp->vlevel - 1] : NULL;
    mpz_add_ui(Z(lp_x), max, vi);
    if (ap)
        mpz_div(Z(lp_x), Z(lp_x), ap->q);
    /* else notional ap->q is 1 */
    /* TODO: calculate mintau(nextt) */
/*  mintau(ZP(lp_mint), vi, nextt);
    mpz_div(Z(lp_x), Z(lp_x), Z(lp_mint));
*/
    mpz_root(Z(lp_x), Z(lp_x), x - 1);
    if (mpz_fits_ulong_p(Z(lp_x))) {
        ulong lim = mpz_get_ui(Z(lp_x));
        if (maxp && maxp < lim)
            return maxp;
        return lim;
    }
    diag_plain();
    keep_diag();
    report("002 %s: outsize limit %Zu\n", diag_buf, Z(lp_x));
    return 0;
}

/* TODO: mod to ulong, use the fast 64-bit mulmod from MPU-mulmod.h.
 * Return p if no inverse exists.
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

void candidate(mpz_t c) {
    keep_diag();
    clock_t t1 = times(NULL);
    report("202 Candidate %Zu (%.2fs)\n", c, seconds(t1));
    if (mpz_cmp(c, max) <= 0) {
        mpz_set(max, c);
        ++seen_best;
    }
}

bool test_prime(mpz_t qq, mpz_t o, ulong ati) {
    mpz_mul_ui(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return _GMP_is_prob_prime(Z(wv_cand));
}

bool test_other(mpz_t qq, mpz_t o, ulong ati, uint t) {
    mpz_mul_ui(Z(wv_cand), qq, ati);
    mpz_add(Z(wv_cand), Z(wv_cand), o);
    return is_tau(Z(wv_cand), t);
}

void walk_v(mpz_t start) {
    mpz_t *q[k];
    mpz_t *m = &levels[level - 1].rq;
    uint t[k];
    mpz_t *aq = &levels[level - 1].aq;
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
/*      else if (t[vi] == 4)
            need_semiprime[nsc++] = vi;
*/
        else
            need_other[noc++] = vi;
    }

#if 0
    if (nqc) {
        /* TODO: special case square walk */
    }
#endif

    if (!mpz_fits_ulong_p(Z(wv_end)))
        fail("TODO: walk_v.end > 2^64");
    ulong end = mpz_get_ui(Z(wv_end));
    for (ulong ati = mpz_get_ui(Z(wv_ati)); ati <= end; ++ati) {
        ++countwi;
        if (times(NULL) >= diagt)
            diag_walk_v(ati, end);
        for (uint ii = 0; ii < inv_count; ++ii) {
            t_mod *ip = &inv[ii];
            if (ati % ip->m == ip->v)
                goto next_ati;
        }
        /* TODO: bail and print somewhere here if 'opt_print' */
        /* note: we have no squares */
        /* TODO: remove me once we handle squares */
        for (uint i = 0; i < nqc; ++i) {
            uint vi = need_square[i];
            if (!test_other(wv_qq[vi], wv_o[vi], ati, t[vi]))
                goto next_ati;
        }
        for (uint i = 0; i < npc; ++i) {
            uint vi = need_prime[i];
            if (!test_prime(wv_qq[vi], wv_o[vi], ati))
                goto next_ati;
        }
        /* TODO: test these in parallel, with optional printme cutoff */
        for (uint i = 0; i < noc; ++i) {
            uint vi = need_other[i];
            if (!test_other(wv_qq[vi], wv_o[vi], ati, t[vi]))
                goto next_ati;
        }
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

bool apply_batch(t_forcep *fp, uint bi) {
    assert(fp->count > bi);
    t_value *vp;
    t_level *lp = &levels[level];
    lp->is_forced = 1;
    lp->bi = bi;

    t_forcebatch *bp = &fp->batch[bi];
    if (!apply_alloc(bp->vi, fp->p, bp->x))
        return 0;
    /* check if we overshot */
    vp = &value[bp->vi];
    if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
        return 0;

    /* TODO: prep this */
    for (uint i = fp->p; i <= bp->vi; i += fp->p) {
        uint x = simple_valuation(i, fp->p) + 1;
        if (!apply_secondary(lp, bp->vi - i, fp->p, x))
            return 0;
        vp = &value[bp->vi - i];
        if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
            return 0;
    }
    for (uint i = fp->p; bp->vi + i < k; i += fp->p) {
        uint x = simple_valuation(i, fp->p) + 1;
        if (!apply_secondary(lp, bp->vi + i, fp->p, x))
            return 0;
        vp = &value[bp->vi + i];
        if (mpz_cmp(vp->alloc[vp->vlevel - 1].q, max) > 0)
            return 0;
    }
    return 1;
}

uint prep_unforced_x(ulong p) {
    t_level *cur_level = &levels[level];
    uint ti = cur_level->ti;

/* TODO: for n=54, we should disallow eg x=9 when t=54, since we
 * will already have tried all x=3 and x=6 before that.
 */

    uint x = divisors[ti].div[cur_level->di];
    uint vi = cur_level->vi;
    t_value *vp = &value[vi];
    t_allocation *ap = (vp->vlevel) ? &vp->alloc[vp->vlevel - 1] : NULL;

    /* pick up any previous unforced x */
    uint unforced_base = (vl_forced) ? VLP(vl_forced - 1)[vi] : 0;
    uint prevx = (vp->vlevel > unforced_base) ? ap->x : 0;
    if (p == 0 && x <= prevx && (ti % prevx) == 0) {
        if (x < prevx)
            return 0;   /* skip this x, we already did the reverse */
        p = ap->p;      /* skip smaller p, we already did the reverse */
    } else if (p == 0)
        p = maxforce[vi];
    /* else we're continuing from known p */

    uint nextt = ti / x;
    /* try p^{x-1} for all p until q_i . p^{x-1} . minrest > max + i */
    ulong limp = limit_p(vi, x, nextt);
    if (limp < p + 1)
        return 0;   /* nothing to do here */
    t_level *prev_level = &levels[level - 1];
    mpz_add_ui(Z(r_walk), max, vi);
    mpz_fdiv_q(Z(r_walk), Z(r_walk), prev_level->aq);
#if 0
/* TODO: support square walk */
    if (prev_level->have_square) {
        /* If we fix a square, expect to actually walk sqrt(r_walk)
         * times number of roots mod cur_level->aq, typically 2^k
         * if there are k primes dividing aq.
         */
        mpz_root(Z(r_walk), Z(r_walk), 2);
        mpz_mul_2exp(Z(r_walk), Z(r_walk), primes_used());
#   if 0
/* TODO: support Pell */
        if (prev_level->have_square > 1)
            mpz_set_ui(Z(r_walk), 0);
#   endif
    }
#endif
    if (gain > 1)
        mpz_mul_ui(Z(r_walk), Z(r_walk), gain);
    if (mpz_fits_ulong_p(Z(r_walk))
        && mpz_get_ui(Z(r_walk)) < limp - p
    ) {
        walk_v(Z(zero));
        return 1;
    }
    cur_level->p = p;
    cur_level->x = x;
    cur_level->limp = limp;
    cur_level->max_at = seen_best;
    /* TODO: do some constant alloc stuff in advance */
    /* TODO: special case for nextt == 1 */
    return 2;
}

/* we emulate recursive calls via the levels[] array */
void recurse(void) {
    if (have_rwalk) {
        walk_v(rwalk_from);
        goto continue_recurse;
    }

    ulong p;
    uint x;
    while (1) {
        ++countr;
        /* recurse deeper */
        {
            uint fi = level - 1;
            t_level *prev_level = &levels[level - 1];
            if (fi < forcedp && (fi == 0 || prev_level->is_forced)) {
                t_forcep *fp = &forcep[fi];
                if (fp->count == 0)
                    goto unforced;
                STOREVL(vl_forced++);
                if (!apply_batch(fp, 0))
                    goto continue_recurse;
                if (times(NULL) >= diagt)
                    diag_plain();
                continue;   /* deeper */
            }
        }
      unforced:
        {
            uint vi = best_v();
            /* TODO: walk_v() directly at previous level, if best_v() would
             * give same result each time.
             */
            if (vi == k) {
                walk_v(Z(zero));
                goto continue_recurse;
            }
            t_value *vp = &value[vi];
            uint ti = (vp->vlevel) ? vp->alloc[vp->vlevel - 1].t : n;
            t_divisors *dp = &divisors[ti];
            if (dp->highdiv == 0) {
                /* is this an error? */
                walk_v(Z(zero));
                goto continue_recurse;
            }
            t_level *cur_level = &levels[level];
            cur_level->vi = vi;
            cur_level->ti = ti;
            cur_level->di = 0;
            goto have_unforced_x;
        }
      continue_unforced_x:
        ++levels[level].di;
      have_unforced_x:
        {
            t_level *cur_level = &levels[level];
            if (cur_level->di >= divisors[cur_level->ti].highdiv)
                goto derecurse;
            switch (prep_unforced_x(0)) {
                /* nothing to do for this x */
                case 0: goto continue_unforced_x;
                /* nothing to do for any x */
                case 1: goto derecurse;
                /* ok, continue for this x */
                case 2: ;
            }
            ++level;
            ++value[cur_level->vi].vlevel;
            goto continue_recurse;
        }
        break;
      derecurse:
        if (levels[level].is_forced)
            --vl_forced;
      continue_recurse:
        --level;
        if (level == 0)
            break;
        {
            t_level *prev_level = &levels[level - 1];
            t_level *cur_level = &levels[level];
            if (cur_level->is_forced) {
                uint fi = level - 1;
                t_forcep *fp = &forcep[fi];

                /* unapply the batch */
                FETCHVL(vl_forced - 1);
                uint bi = cur_level->bi + 1;
                if (bi >= fp->count)
                    goto derecurse;
                if (fp->batch[bi].x == 0) {
                    cur_level->is_forced = 0;
                    goto unforced;
                }
                if (!apply_batch(fp, bi))   /* ++level */
                    goto continue_recurse;
                if (times(NULL) >= diagt)
                    diag_plain();
            } else {
                --value[cur_level->vi].vlevel;
                ulong p = cur_level->p;
                /* recalculate limit if we have an improved maximum */
                if (seen_best > cur_level->max_at)
                    if (!prep_unforced_x(p))
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
                if (!apply_alloc(cur_level->vi, p, cur_level->x)) {
                    /* is this an error? */
                    if (times(NULL) >= diagt)
                        diag_plain();
                    goto continue_recurse;
                }
                if (times(NULL) >= diagt)
                    diag_plain();
            }
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
            gain = strtoul(&arg[2], NULL, 0);
        else if (arg[1] == 'p')
            set_cap(&arg[2]);
        else if (arg[1] == 'r')
            runid = strtoul(&arg[2], NULL, 10);
        else if (arg[1] == 'f')
            force_all = strtoul(&arg[2], NULL, 10);
        else if (strncmp("-o", arg, 2) == 0)
            opt_print = 1;
        else if (strncmp("-d", arg, 2) == 0)
            debug = 1;
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
    if (rstack)
        insert_stack();
    recurse();
    keep_diag();

    clock_t tz = times(NULL);
    report("367 coul(%u, %u): recurse %lu, walk %lu, walkc %lu (%.2fs)\n",
            n, k, countr, countw, countwi, seconds(tz));
    if (seen_best)
        report("200 f(%u, %u) = %Zu (%.2fs)\n", n, k, max, seconds(tz));
    done();
    return 0;
}
