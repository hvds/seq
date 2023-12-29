#include "unit.h"
#include "diag.h"
#include "sieve.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <signal.h>
#include <time.h>
#include <errno.h>
#include <assert.h>

#include <sys/time.h>
#include <sys/resource.h>

/* Search recursively for representations of p/q as a sum of <depth> distinct
   unit fractions. */

uint maxdepth;

/* a prime power p^k, used for representing a factorization of some q
 * and iterating divisors of q^2 */
typedef struct fac_s {
    uint p;     /* prime */
    uint k;     /* power */
    uint i;     /* position of divisors iterator */
    mpz_t pk;   /* power of divisors iterator */
} fac_t;

fac_t* f;       /* factors of q */
uint fcount;    /* number of factors in f */
uint fsize;     /* malloced size of f */

/* Structure recording information about a rational p/q at a given depth in
 * the recursion.
 * Note that we always set rat[ri].r = rat[ri-1].r - rat[ri].qmin: both
 * qmin and max in rat[ri] refer to the range of unit fractions to be
 * trial-subtracted from rat[ri-1].r. We leave rat[0].qmin = 1/0 to get the
 * correct limits.
 */
typedef struct rat_s {
    mpq_t r;        /* the rational itself, p/q */
    mpq_t qmin;     /* 1/min, the last unit fraction to be trial-subtracted */
    mpz_t max;      /* final den. for trial subtraction */
} rat_t;
rat_t *rat;

#define RI(ri) rat[ri].r
#define QI(ri) mpq_denref(RI(ri))
#define PI(ri) mpq_numref(RI(ri))
#define MINI(ri) mpq_denref(rat[ri].qmin)
#define MAXI(ri) rat[ri].max
uint recover_ri;

mpz_t f2_mod, f2_min, nextdiv, ztemp;
mpq_t qtemp;

#define SIEVESIZE 1000000

char *rpath = NULL; /* path to log file */
FILE *rfp = NULL;   /* file handle to log file */

#define DIAG 1
#define LOG 600
double diag_delay = DIAG, log_delay = LOG;
ulong count_s2 = 0, count_p3 = 0, count_q3 = 0;
char *diag_buf = NULL;
uint diag_buf_size = 0;
uint diag_bufp = 0;

double t0 = 0;
struct rusage rusage_buf;
static inline double seconds(void) {
    getrusage(RUSAGE_SELF, &rusage_buf);
    return (double)rusage_buf.ru_utime.tv_sec
            + (double)rusage_buf.ru_utime.tv_usec / 1000000
            - t0;
}
timer_t diag_timerid, log_timerid;
volatile bool need_work, need_diag, need_log;

void reset_buf(void) {
    diag_bufp = 0;
}

void append_buf(char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    uint size = diag_bufp + gmp_vsnprintf(
            diag_buf + diag_bufp, diag_buf_size - diag_bufp, fmt, ap);
    va_end(ap);
    if (size > diag_buf_size) {
        diag_buf = realloc(diag_buf, size + 32);
        diag_buf_size = size + 32;
        va_start(ap, fmt);
        gmp_vsprintf(diag_buf + diag_bufp, fmt, ap);
        va_end(ap);
    }
    diag_bufp = size;
    return;
}

void diag_plain(uint ri) {
    double t1 = seconds();
    reset_buf();
    append_buf("%Qu %lu/%lu/%lu [", RI(0), count_s2, count_p3, count_q3);
    for (uint i = 1; i <= ri; ++i) {
        if (i > 1) append_buf(" ");
        append_buf("%Zu", MINI(i));
    }
    append_buf("] (%.2fs)", t1);
    if (need_diag) {
        diag("%s", diag_buf);
        need_diag = 0;
    }
    if (rfp && need_log) {
        fprintf(rfp, "305 %s\n", diag_buf);
        need_log = 0;
    }
    need_work = 0;
    return;
}

void resize_fac(uint size) {
    f = (fac_t *)realloc(f, size * sizeof(fac_t));
    for (uint i = fsize; i < size; ++i)
        ZINIT(f[i].pk);
    fsize = size;
}

/* print details of the structure rat[ri] to stderr */
void diagnose(uint ri) {
    rat_t* r = &rat[ri];
    gmp_fprintf(stderr, "rat[%d]: r = %Qd; qmin = %Qd; max = %Zd\n",
            ri, r->r, r->qmin, r->max);
    gmp_fprintf(stderr, "\n");
}

void fail(char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    if (rfp)
        fclose(rfp);
    exit(1);
}

void handle_sig(int sig) {
    need_work = 1;
    if (sig == SIGUSR1)
        need_diag = 1;
    else
        need_log = 1;
}

void init_time(void) {
    struct sigaction sa;
    struct sigevent sev;
    struct itimerspec diag_timer, log_timer;

    sa.sa_handler = &handle_sig;
    sa.sa_flags = SA_RESTART;
    sigemptyset(&sa.sa_mask);

    if (diag_delay) {
        if (sigaction(SIGUSR1, &sa, NULL))
            fail("Could not set USR1 handler: %s\n", strerror(errno));
        sev.sigev_notify = SIGEV_SIGNAL;
        sev.sigev_signo = SIGUSR1;
        sev.sigev_value.sival_ptr = &diag_timerid;
        if (timer_create(CLOCK_PROCESS_CPUTIME_ID, &sev, &diag_timerid)) {
            /* guess that the CPUTIME clock is not supported */
            if (timer_create(CLOCK_REALTIME, &sev, &diag_timerid))
                fail("Could not create diag timer: %s\n", strerror(errno));
        }
        diag_timer.it_value.tv_sec = diag_delay;
        diag_timer.it_value.tv_nsec = 0;
        diag_timer.it_interval.tv_sec = diag_delay;
        diag_timer.it_interval.tv_nsec = 0;
        if (timer_settime(diag_timerid, 0, &diag_timer, NULL))
            fail("Could not set diag timer: %s\n", strerror(errno));
    }

    if (log_delay) {
        if (sigaction(SIGUSR2, &sa, NULL))
            fail("Could not set USR2 handler: %s\n", strerror(errno));
        sev.sigev_notify = SIGEV_SIGNAL;
        sev.sigev_signo = SIGUSR2;
        sev.sigev_value.sival_ptr = &log_timerid;
        if (timer_create(CLOCK_PROCESS_CPUTIME_ID, &sev, &log_timerid)) {
            /* guess that the CPUTIME clock is not supported */
            if (timer_create(CLOCK_REALTIME, &sev, &log_timerid))
                fail("Could not create log timer: %s\n", strerror(errno));
        }
        log_timer.it_value.tv_sec = log_delay;
        log_timer.it_value.tv_nsec = 0;
        log_timer.it_interval.tv_sec = log_delay;
        log_timer.it_interval.tv_nsec = 0;
        if (timer_settime(log_timerid, 0, &log_timer, NULL))
            fail("Could not set log timer: %s\n", strerror(errno));
    }
}

/* Parse a "305" log line for initialization.
 * Input string should point after the initial "305 ".
 */
void parse_305(char *s) {
    double dtime;

    uint p, q;
    p = strtoul(s, &s, 10);
    assert(s[0] == '/');
    q = strtoul(s + 1, &s, 10);
    assert(s[0] == ' ');
    count_s2 = strtoul(s + 1, &s, 10);
    assert(s[0] == '/');
    count_p3 = strtoul(s + 1, &s, 10);
    assert(s[0] == '/');
    count_q3 = strtoul(s + 1, &s, 10);
    assert(s[0] == ' ');
    assert(s[1] == '[');
    s += 2;

    recover_ri = 0;
    while (1) {
        if (s[0] == ']')
            break;
        ulong mini = strtoul(s, &s, 10);
        ++recover_ri;
        mpz_set_ui(MINI(recover_ri), mini);
        if (s[0] != ' ')
            break;
        ++s;
    }

    assert(s[0] == ']');
    assert(s[1] == ' ');
    assert(s[2] == '(');
    double t = strtod(s + 3, &s);
    t0 = -t;
    assert(s[0] == 's');
    assert(s[1] == ')');
}

char *recover(FILE *fp) {
    char *last305 = NULL;
    char *curbuf = NULL;
    size_t len = 120, len305 = 0;

    while (1) {
        ssize_t nread = getline(&curbuf, &len, fp);
        if (nread <= 0) {
            if (errno == 0)
                break;
            fail("error reading %s: %s", rpath, strerror(errno));
        }
        if (curbuf[nread - 1] != '\n'
                || memchr(curbuf, 0, nread) != NULL) {
            /* corrupt line, file should be truncated */
            off_t offset = ftello(fp);
            if (offset == -1)
                fail("could not ask offset: %s", strerror(errno));
            /* not ftruncate(), we are open only for reading */
            if (truncate(rpath, offset - nread) != 0)
                fail("could not truncate %s to %lu: %s", rpath, offset - nread,
                        strerror(errno));
            break;
        }
        if (strncmp("305 ", curbuf, 4) == 0) {
            char *t = last305;
            last305 = curbuf;
            curbuf = t;
            size_t lt = len305;
            len305 = len;
            len = lt;
        } else if (strncmp("000 ", curbuf, 4) == 0) {
            /* comment */
        } else
            fail("unexpected log line %.3s in %s", curbuf, rpath);
    }

    free(curbuf);
    return last305;
}

void init_unit(uint max, char* report_path) {
    char *recover_line = NULL;

    rpath = report_path;
    if (rpath) {
        FILE *fp = fopen(rpath, "r");
        if (fp) {
            recover_line = recover(fp);
            fclose(fp);
        }

        rfp = fopen(rpath, "a");
        if (rfp == NULL)
            fail("%s: %s", rpath, strerror(errno));
        setlinebuf(rfp);
    }
    t0 = seconds();
    init_diag();    /* ignore result: worst case we lose ^Z handling */
    init_time();
    extend_sieve(SIEVESIZE);
    maxdepth = max;
    rat = malloc(maxdepth * sizeof(rat_t));
    for (uint i = 0; i < maxdepth; ++i) {
        QINIT(rat[i].r);
        QINIT(rat[i].qmin);
        mpz_set_ui(mpq_numref(rat[i].qmin), 1);    /* never changes */
        ZINIT(rat[i].max);
    }
    mpz_set_ui(mpq_denref(rat[0].qmin), 0);     /* base, never changes */
    fsize = 0;
    f = NULL;
    resize_fac(10);
    ZINIT(f2_mod);
    ZINIT(f2_min);
    ZINIT(nextdiv);
    ZINIT(ztemp);
    QINIT(qtemp);

    if (recover_line) {
        parse_305(recover_line + 4);
        free(recover_line);
    }
}

void done_unit(void) {
    keep_diag();
    free(diag_buf);
    for (uint i = 0; i < maxdepth; ++i) {
        QCLEAR(rat[i].r);
        QCLEAR(rat[i].qmin);
        ZCLEAR(rat[i].max);
    }
    free(rat);
    for (uint i = 0; i < fsize; ++i)
        ZCLEAR(f[i].pk);
    free(f);
    ZCLEAR(f2_mod);
    ZCLEAR(f2_min);
    ZCLEAR(nextdiv);
    ZCLEAR(ztemp);
    QCLEAR(qtemp);
}

/* set rat[ri].f to the factors of QI(ri), and initialize to walk divisors */
void init_divs(uint ri) {
    uint ni = 0;
    mpz_set(ztemp, QI(ri));
    if (mpz_even_p(ztemp)) {
        f[ni].p = 2;
        f[ni].k = 0;
        f[ni].i = 0;
        mpz_set_ui(f[ni].pk, 1);
        while (mpz_even_p(ztemp)) {
            f[ni].k += 2;
            mpz_divexact_ui(ztemp, ztemp, 2);
        }
        ++ni;
    }
    uint p = 3;
    while (mpz_cmp_ui(ztemp, 1) > 0) {
        if (mpz_divisible_ui_p(ztemp, p)) {
            if (ni + 1 >= fsize)
                resize_fac(fsize + 8);
            f[ni].p = p;
            f[ni].k = 0;
            f[ni].i = 0;
            mpz_set_ui(f[ni].pk, 1);
            while (mpz_divisible_ui_p(ztemp, p)) {
                f[ni].k += 2;
                mpz_divexact_ui(ztemp, ztemp, p);
            }
            ++ni;
        }
        p += 2;
    }
    fcount = ni;
    if (fcount == 0) {
        mpz_set_ui(f[0].pk, 1);
        f[0].i = 0;
    }
    return;
}

/* For squares, we want divisors d of q^2 such that q/d is a square.
 * For each odd prime power in q, reduce the exponent by 1 and multiply
 * the prime into an accumulator; finally, set the accumulator to be the
 * base of all divisors.
 * Then replace (p, 2k) with (p^2, k) for the iteration.
 * If any prime p=4n+3 divides q to an odd power, abort and return false.
 */
bool init_sq_divs(uint ri) {
    uint ni = 0;
    mpz_set(ztemp, QI(ri));
    mpz_set_ui(nextdiv, 1); /* accumulate odd powers here */
    if (mpz_even_p(ztemp)) {
        uint k = 0;
        while (mpz_even_p(ztemp)) {
            ++k;
            mpz_divexact_ui(ztemp, ztemp, 2);
        }
        if (k & 1) {
            --k;
            mpz_mul_ui(nextdiv, nextdiv, 2);
        }
        if (k) {
            f[ni].p = 4;
            f[ni].k = k;
            f[ni].i = 0;
            mpz_set_ui(f[ni].pk, 1);
            ++ni;
        }
    }
    uint pi = 1;
    while (mpz_cmp_ui(ztemp, 1) > 0) {
        uint p = prime[pi];
        if (mpz_divisible_ui_p(ztemp, p)) {
          final_prime:
            if (ni + 1 >= fsize)
                resize_fac(fsize + 8);
            uint k = 0;
            while (mpz_divisible_ui_p(ztemp, p)) {
                ++k;
                mpz_divexact_ui(ztemp, ztemp, p);
            }
            if (k & 1) {
                if ((p & 3) == 3)
                    return 0;
                --k;
                mpz_mul_ui(nextdiv, nextdiv, p);
            }
            if (k) {
                f[ni].p = p * p;
                f[ni].k = k;
                f[ni].i = 0;
                mpz_set_ui(f[ni].pk, 1);
                ++ni;
            }
        } else if (mpz_cmp_ui(ztemp, p * p) < 0) {
            p = mpz_get_ui(ztemp);
            goto final_prime;
        }
        ++pi;
        if (pi >= nprime)
            extend_sieve(maxsieve + SIEVESIZE);
    }
    fcount = ni;
    if (fcount == 0) {
        mpz_set_ui(f[0].pk, 1);
        f[0].i = 0;
        return 1;
    }

    if (mpz_cmp_ui(nextdiv, 1) > 0) {
        ++fcount;
        if (fcount > fsize)
            resize_fac(fsize + 8);
        f[ni].p = 1;
        f[ni].k = 0;
        f[ni].i = 0;
        mpz_set(f[ni].pk, nextdiv);
        for (uint fi = 0; fi < ni; ++fi)
            mpz_set(f[fi].pk, nextdiv);
    }
    return 1;
}

/* Set 'nextdiv' to the next divisor of QI(ri)^2.
 * init_divs() must have been called before this; divisors are not returned
 * in order.
 * Returns false when all divisors have been walked.
 */
bool next_div(mpz_t max) {
    uint fi = 0;
    mpz_set(nextdiv, f[0].pk);
    while (1) {
        if (fi == fcount) {
            if (fi == 0)
                return f[0].i++ == 0 ? 1 : 0;
            if (f[0].i < f[0].k + 2) {
                f[0].i = f[0].k + 2;
                return 1;
            }
            return 0;
        }
        ++f[fi].i;
        if (f[fi].i > f[fi].k) {
            ++fi;
        } else {
            mpz_mul_ui(f[fi].pk, f[fi].pk, f[fi].p);
            if (mpz_cmp(f[fi].pk, max) <= 0)
                break;
            ++fi;
        }
    }
    while (fi > 0) {
        --fi;
        mpz_set(f[fi].pk, f[fi + 1].pk);
        f[fi].i = 0;
    }
    return 1;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators > m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q.
 */
bool find_s2(uint ri) {
    init_divs(ri);
    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    while (next_div(QI(ri))) {
        if (mpz_cmp(f2_min, nextdiv) >= 0)
            continue;
        if (mpz_cmp(QI(ri), nextdiv) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, nextdiv, PI(ri)))
            return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * unit fractions with denominators >= m = MINI(ri).
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q <= d < q.
 */
bool find_m2(uint ri) {
    init_divs(ri);
    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    while (next_div(QI(ri))) {
        if (mpz_cmp(f2_min, nextdiv) > 0)
            continue;
        if (mpz_cmp(QI(ri), nextdiv) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (mpz_congruent_p(f2_mod, nextdiv, PI(ri)))
            return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth> distinct
 * unit fractions with denominators > m = MINI(ri).
 * Requires depth >= 3.
 */
bool find_sn(uint ri, uint depth) {
    uint rj = ri + 1;
    if (mpz_cmp_ui(PI(ri), 1) == 0)
        return 1;
    mpz_mul_ui(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    if (mpz_cmp(MINI(ri), MINI(rj)) >= 0)
        mpz_add_ui(MINI(rj), MINI(ri), 1);
    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_sub(RI(rj), RI(ri), rat[rj].qmin);
        if ((depth == 3) ? find_s2(rj) : find_sn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth>
 * unit fractions, not necessarily distinct, with denominators >= m = MINI(ri).
 * Requires depth >= 3.
 */
bool find_mn(uint ri, uint depth) {
    uint rj = ri + 1;
    if (mpz_cmp_ui(PI(ri), 1) == 0)
        return 1;
    mpz_mul_ui(MAXI(rj), QI(ri), depth);
    mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
    mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
    if (mpz_cmp(MINI(ri), MINI(rj)) > 0)
        mpz_set(MINI(rj), MINI(ri));
    if (mpz_sgn(MINI(rj)) == 0)
        mpz_set_ui(MINI(rj), 1);
    while (mpz_cmp(MINI(rj), MAXI(rj)) <= 0) {
        mpq_sub(RI(rj), RI(ri), rat[rj].qmin);
        if ((depth == 3) ? find_m2(rj) : find_mn(rj, depth - 1))
            return 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
    }
    return 0;
}

/* Return TRUE if p/q can be represented as fewer than c distinct unit
 * fractions.
 */
bool better_set(mpq_t r, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_s2(0))
        return 1;
    if (c == 3)
        return 0;
    if (find_sn(0, c - 1))
        return 1;
    return 0;
}

/* Return TRUE if p/q can be represented as fewer than c unit fractions
 * (not necessarily distinct).
 */
bool better_multi(mpq_t r, uint c) {
    if (c <= 1)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    if (c == 2)
        return 0;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_m2(0))
        return 1;
    if (c == 3)
        return 0;
    if (find_mn(0, c - 1))
        return 1;
    return 0;
}

uint find_set(mpq_t r) {
    if (mpq_sgn(r) == 0)
        return 0;
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return 1;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_s2(0))
        return 2;
    for (uint c = 3; 1; ++c) {
        if (find_sn(0, c))
            return c;
    }
    /* not reached */
}

/* return floor(r), assumed to fit in uint */
uint mpq_int(mpq_t r) {
    if (mpz_cmp(mpq_numref(r), mpq_denref(r)) < 0)
        return 0;
    mpz_fdiv_q(ztemp, mpq_numref(r), mpq_denref(r));
    return mpz_get_ui(ztemp);
}

uint find_multi(mpq_t r) {
    if (mpq_sgn(r) == 0)
        return 0;
    uint unit = mpq_int(r);
    if (unit) {
        mpz_fdiv_r(mpq_numref(r), mpq_numref(r), mpq_denref(r));
        if (mpq_sgn(r) == 0)
            return unit;
    }
    if (mpz_cmp_ui(mpq_numref(r), 1) == 0)
        return unit + 1;
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (find_m2(0))
        return unit + 2;
    for (uint c = 3; 1; ++c) {
        if (find_mn(0, c))
            return unit + c;
    }
    /* not reached */
}

/* quick div-3 test: return true if the bit-pattern of z matches ...110*,
 * a quick-and-dirty way to spot most cases in which the square-free
 * residue of z is divisible by a prime of the form 4n+3.
 * In tests, this missed about 4.5% of cases.
 */
static inline bool test_and3(mpz_t z) {
    mp_limb_t z0 = mpz_getlimbn(z, 0);
    return (z0 & ((z0 ^ (z0 & (z0 - 1))) << 1)) ? 1 : 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of two distinct
 * square unit fractions with denominators > (m = MINI(ri))^2.
 * Given r = p/q, this is true precisely if there exists a divisor d of q^2
 * with d == -q (mod p) and mp-q < d < q such that a = (q + d)/p and
 * b = (q + q^2/d)/p are both perfect squares. Since this implies b = q/d a,
 * we require q/d itself to be a square.
 *
 * Note that no solution is possible if the square-free residue of p is
 * divisible by any prime of the form 4k+3; however, unlike q, we have no
 * constraint on how hard p may be to factorize, so we do only a quick
 * partial check for this.
 */
bool find_square_s2(uint ri) {
    if (need_work)
        diag_plain(ri);
    if (test_and3(PI(ri))) {
        ++count_p3;
        return 0;
    }
    if (test_and3(QI(ri)) || !init_sq_divs(ri)) {
        ++count_q3;
        return 0;
    }
    ++count_s2;

    mpz_ui_sub(f2_mod, 0, QI(ri));
    mpz_mul(f2_min, PI(ri), MINI(ri));
    mpz_mul(f2_min, f2_min, MINI(ri));
    mpz_sub(f2_min, f2_min, QI(ri));

    /* if the divisors were sorted, we could use a binary chop to find the
     * start point, and know that floor(dcount / 2) is the end point */
    while (next_div(QI(ri))) {
        if (mpz_cmp(f2_min, nextdiv) >= 0)
            continue;
        if (mpz_cmp(QI(ri), nextdiv) <= 0)
            continue;
        /* it might be faster to normalize both mod p, and check equality */
        if (!mpz_congruent_p(f2_mod, nextdiv, PI(ri)))
            continue;
        mpz_add(ztemp, QI(ri), nextdiv);
        mpz_divexact(ztemp, ztemp, PI(ri));
        if (!mpz_perfect_square_p(ztemp))
            continue;
        /* we have a solution, so store the relevant values for reporting */
        mpz_sqrt(MINI(ri + 1), ztemp);
        mpz_mul(ztemp, ztemp, QI(ri));
        mpz_divexact(ztemp, ztemp, nextdiv);
        mpz_sqrt(MINI(ri + 2), ztemp);
        keep_diag();
        gmp_printf("%Qu: %u [", RI(0), ri + 2);
        for (uint i = 1; i <= ri + 2; ++i) {
            if (i > 1) gmp_printf(" ");
            gmp_printf("%Zu", MINI(i));
        }
        printf("] (%.2fs)\n", seconds());
        return 1;
    }
    return 0;
}

/* Return TRUE if r = RI(ri) can be expressed as the sum of <depth> distinct
 * square unit fractions with denominators > (m = MINI(ri))^2.
 * Requires depth >= 3.
 */
bool find_square_sn(uint depth, uint recover) {
    uint ri = 0, rj;

    if (recover) {
        for (; ri + 1 <= recover; ++ri) {
            rj = ri + 1;
            /* max := floor(sqrt(q * depth / p)) */
            mpz_mul_ui(MAXI(rj), QI(ri), depth - ri);
            mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
            mpz_sqrt(MAXI(rj), MAXI(rj));

            mpq_mul(qtemp, rat[rj].qmin, rat[rj].qmin);
            mpq_sub(RI(rj), RI(ri), qtemp);
        }
        --ri;
        goto loop_test;
    }

    while (1) {
        rj = ri + 1;
        /* max := floor(sqrt(q * depth / p)) */
        mpz_mul_ui(MAXI(rj), QI(ri), depth - ri);
        mpz_cdiv_q(MAXI(rj), MAXI(rj), PI(ri));
        mpz_sqrt(MAXI(rj), MAXI(rj));

        /* min := max(ceil(sqrt(q / p)), prevmin + 1) */
        mpz_cdiv_q(MINI(rj), QI(ri), PI(ri));
        mpz_sqrtrem(MINI(rj), ztemp, MINI(rj));
        if (mpz_sgn(ztemp) > 0)
            mpz_add_ui(MINI(rj), MINI(rj), 1);
        if (mpz_cmp(MINI(rj), MINI(ri)) <= 0)
            mpz_add_ui(MINI(rj), MINI(ri), 1);
        goto loop_level;

      next_level:
        rj = ri + 1;
        mpz_add_ui(MINI(rj), MINI(rj), 1);
      loop_level:
        if (mpz_cmp(MINI(rj), MAXI(rj)) > 0) {
            if (ri--)
                goto next_level;
            break;
        }
        mpq_mul(qtemp, rat[rj].qmin, rat[rj].qmin);
        mpq_sub(RI(rj), RI(ri), qtemp);
      loop_test:
        if (depth - ri == 3) {
            if (find_square_s2(rj))
                return 1;
            goto next_level;
        }
        ++ri;
    }
    return 0;
}

uint find_square_set(mpq_t r, uint min_depth, uint max_depth) {
    if (min_depth <= 0) {
        if (mpq_sgn(r) == 0)
            return 0;
        keep_diag();
        gmp_printf("%Qu: not 0 (%.2fs)\n", r, seconds());
        t0 += seconds();
    }
    if (min_depth <= 1) {
        if (mpz_cmp_ui(mpq_numref(r), 1) == 0
            && mpz_perfect_square_p(mpq_denref(r))
        )
            return 1;
        keep_diag();
        gmp_printf("%Qu: not 1 (%.2fs)\n", r, seconds());
        t0 += seconds();
    }
    mpq_get_num(PI(0), r);
    mpq_get_den(QI(0), r);
    if (min_depth <= 2) {
        if (find_square_s2(0))
            return 2;
        keep_diag();
        gmp_printf("%Qu: not 2 (%.2fs) c=%lu/%lu/%lu\n", r, seconds(), count_s2, count_p3, count_q3);
        count_s2 = 0;
        count_p3 = 0;
        count_q3 = 0;
        t0 += seconds();
    }
    for (uint c = (min_depth > 3) ? min_depth : 3; c <= max_depth; ++c) {
        if (find_square_sn(c, recover_ri))
            return c;
        recover_ri = 0;
        keep_diag();
        gmp_printf("%Qu: not %u (%.2fs) c=%lu/%lu/%lu\n", r, c, seconds(), count_s2, count_p3, count_q3);
        count_s2 = 0;
        count_p3 = 0;
        count_q3 = 0;
        t0 += seconds();
    }
    return 0;
    /* not reached */
}

