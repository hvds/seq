#include "breadth.h"

extern mpq_t r;        /* the rational we're testing */
extern mpq_t rone;     /* handy 1 */
extern void report_breadth(ulong gen, ulong count);

mpq_t limit;    /* r + 1/r: no point queueing anything smaller than this */
mpq_t curq;     /* rational from queue under consideration */
mpq_t limit_calc;   /* scratch space for limit calculation */
mpz_t limit_div;    /* scratch space for limit calculation */

/* bitsizes for simplest in the generation */
ulong best_bits_num, best_bits_den;

rat_array_t ra1, ra2;

inline rat_array_t *choose_cur(ulong gen) {
    return (gen & 1) ? &ra1 : &ra2;
}
inline rat_array_t *choose_next(ulong gen) {
    return (gen & 1) ? &ra2 : &ra1;
}

void array_init(rat_array_t *a) {
    a->count = (size_t)0;
    a->actual = (size_t)0;
    a->size = (size_t)1024; /* arbitrary */
    a->space = (pack_t*)malloc(a->size * sizeof(pack_t));
}

void array_reset(rat_array_t *a) {
    a->count = (size_t)0;
    a->actual = (size_t)0;
}

void array_resize(rat_array_t *a, size_t newsize) {
    if (newsize > a->size) {
        size_t bigger = a->size * 3 / 2;
        size_t other_total = ra1.size + ra2.size - a->size;
        if (bigger < newsize)
            bigger = newsize;
        if ((other_total + bigger) > (size_t)HARD_LIMIT) {
            bigger = HARD_LIMIT - other_total;
            if (bigger < newsize) {
                fprintf(stderr, "request for %lu + %lu exceeds hard limit\n",
                        newsize, other_total);
                exit(1);
            }
        }
        a->space = (pack_t*)realloc(a->space, bigger * sizeof(pack_t));
        a->size = bigger;
        if (!a->space) {
            fprintf(stderr, "out of memory allocating %lu\n", (ulong)bigger);
            exit(1);
        }
    }
}

void array_free(rat_array_t *a) {
    free(a->space);
    a->space = (pack_t *)NULL;
    a->size = 0;
}

void array_push(rat_array_t *a, mpq_t new) {
    size_t s = mpq_packsize(new);
    array_resize(a, a->count + s);
    store_packed((mpq_pack_t*)&(a->space[a->count]), new);
    a->count += s;
    ++a->actual;
}

/* Initialization, must call this before any call to breadth_one */
void init_breadth(mpq_t r) {
    /* r + 1/r == (r^2 + 1)/r */
    QINIT(&limit, "limit");
    mpq_mul(limit, r, r);
    mpq_add_ui(limit, (ulong)1);
    mpq_div(limit, limit, r);
    QINIT(&curq, "curq");
    QINIT(&limit_calc, "limit_calc");
    ZINIT(&limit_div, "limit_div");
    array_init(&ra1);
    array_init(&ra2);
    array_push(choose_next(0), rone);  /* gen 0 next will be gen 1 current */
}

/* Cleanup, we should be valgrind leak-check clean. */
void finish_breadth(void) {
    array_free(&ra2);
    array_free(&ra1);
    ZCLEAR(&limit_div, "limit_div");
    QCLEAR(&limit_calc, "limit_calc");
    QCLEAR(&curq, "curq");
    QCLEAR(&limit, "limit");
}

/* Given the pending queue for generation I<gen> in the C<choose_cur(gen)>
 * array, check that list for solutions, pushing any new values to check
 * onto the C<choose_next(gen)> array.
 * Returns boolean TRUE if a solution was found (in which case the contents
 * of the C<choose_next(gen)> array may be incomplete), else FALSE.
 */
bool breadth_one(ulong gen) {
    bool solved = 0;
    rat_array_t *cur = choose_cur(gen);
    rat_array_t *next = choose_next(gen);
    ulong curp = 0;
    ulong bits_num, bits_den;

    array_reset(next);
    best_bits_num = (ulong)0;
    best_bits_den = (ulong)0;
    while (curp < cur->count) {
        ulong count = 0;
        mpq_pack_t *pack = (mpq_pack_t*)&(cur->space[curp]);
        fetch_packed(pack, curq);
        curp += packed_size(pack);
        mpq_inv(curq, curq);

        if (mpq_cmp(curq, limit) > 0) {
            /* count = floor((curq - limit) / r) */
            mpq_sub(limit_calc, curq, limit);
            mpq_div(limit_calc, limit_calc, r);
            mpz_fdiv_q(limit_div,
                    mpq_numref(limit_calc), mpq_denref(limit_calc));
            count = mpz_strict_get_ui(limit_div);

            /* curq -= count * r */
            mpq_set(limit_calc, r);
            mpz_mul_ui(mpq_numref(limit_calc),
                    mpq_numref(limit_calc), count);
            mpq_canonicalize(limit_calc);
            mpq_sub(curq, curq, limit_calc);
        }

        /* now we've guaranteed curq - r <= limit, so all values we find
         * in the loop are useful.
         */
        /* while ((curq -= r) > 0) { ... } */
        while (1) {
            mpq_sub(curq, curq, r);
            if (mpq_sgn(curq) <= 0)
                break;
            ++count;
            if (mpq_equal(curq, rone)) {
                report_breadth(gen, count);
                solved = 1;
            }
            if (!solved) {
                array_push(next, curq);
                bits_num = mpz_bitsize(mpq_numref(curq));
                if (!best_bits_num || best_bits_num >= bits_num) {
                    best_bits_num = bits_num;
                    bits_den = mpz_bitsize(mpq_denref(curq));
                    if (!best_bits_den || best_bits_den > bits_den)
                        best_bits_den = bits_den;
                }
            }
        }
    }
    return solved;
}
