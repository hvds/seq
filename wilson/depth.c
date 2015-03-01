#include "depth.h"

mpq_t limit_calc;
mpz_t limit_div;

extern mpq_t r;
extern mpq_t rone;
extern mpq_t limit;

typedef struct {
    mpq_t q;
    ulong count;
    ulong actual;
    ulong best_bits_num;
    ulong best_bits_den;
} frame_t;

frame_t *stack;
ulong max_depth;
bool solved;
ulong count;

void init_depth(ulong depth) {
    ulong i;
    max_depth = depth;
    /* we count from 0 to n */
    stack = (frame_t *)malloc((depth + 1) * sizeof(frame_t));
    for (i = 0; i <= max_depth; ++i) {
        QINIT(&stack[i].q, "stack[%lu].q", i);
        stack[i].actual = (ulong)0;
        stack[i].best_bits_num = (ulong)0;
        stack[i].best_bits_den = (ulong)0;
    }
    solved = 0;
    QINIT(&limit_calc, "depth limit_calc");
    ZINIT(&limit_div, "depth limit_div");
}

void finish_depth(void) {
    ulong i;
    ZCLEAR(&limit_div, "depth limit_div");
    QCLEAR(&limit_calc, "depth limit_calc");
    for (i = 0; i <= max_depth; ++i)
        QCLEAR(&stack[i].q, "stack[%lu].q", i);
    free(stack);
}

void report_depth(ulong depth) {
    ulong i;
    gmp_printf("%Qd is solved in %lu steps by [ ", r, depth);
    for (i = depth; i > 0; --i)
        printf("%lu ", stack[i].count);
    printf("]\n");
}

extern double timing(void);
void report_progress(ulong depth) {
    ulong i;
    ulong lim = (max_depth > 20) ? 20 : max_depth;
    gmp_printf("%Qd: tried %lu values (%.2fs) ", r, count, timing());
    for (i = 0; i < lim; ++i) {
        frame_t *f = &stack[i];
        printf("%lu%s", f->actual, (i < lim - 1) ? " ": "\n");
    }
}

void report_final(void) {
    ulong i;
    for (i = 0; i < max_depth; ++i) {
        frame_t *f = &stack[i];
        gmp_printf("%Qd - g > %lu: %lu, best_bits %lu/%lu\n",
                r, i + 1, f->actual, f->best_bits_num, f->best_bits_den);
    }
}

bool try_depth(ulong depth) {
    frame_t *cur = &(stack[depth]);
    frame_t *next = &(stack[depth + 1]);
    int final = (depth + 1 >= max_depth) ? 1 : 0;
    bool locally_solved = 0;
    ulong bits_num, bits_den;

    ++count;
    if ((count & 0xfffffff) == 0)
        report_progress(depth);

    mpq_inv(next->q, cur->q);

    /* if final, we only care whether we can reach 1 */
    if (final) {
        mpq_sub(next->q, next->q, rone);
        if (mpq_sgn(next->q) < 0)
            return 0;
        mpq_div(next->q, next->q, r);
        if (mpz_cmp_ui(mpq_denref(next->q), (ulong)1) != 0)
            return 0;
        next->count = mpz_strict_get_ui(mpq_numref(next->q));
        report_depth(depth + 1);
        return 1;
    }

    if (mpq_cmp(next->q, limit) <= 0) {
        next->count = (ulong)0;
    } else {
        /* count = floor((q - limit) / r) */
        mpq_sub(limit_calc, next->q, limit);
        mpq_div(limit_calc, limit_calc, r);
        mpz_fdiv_q(limit_div, mpq_numref(limit_calc), mpq_denref(limit_calc));
        next->count = mpz_strict_get_ui(limit_div);

        /* q -= count * r */
        mpq_set(limit_calc, r);
        mpz_mul_ui(mpq_numref(limit_calc), mpq_numref(limit_calc), next->count);
        mpq_canonicalize(limit_calc);
        mpq_sub(next->q, next->q, limit_calc);
    }

    /* now we've guaranteed curq - r <= limit, so all values we find are
     * useful.
     */
    /* while ((q -= r) > 0) { ... } */
    while (1) {
        mpq_sub(next->q, next->q, r);
        if (mpq_sgn(next->q) <= 0)
            return locally_solved;
        ++next->count;
        if (mpq_equal(next->q, rone)) {
            report_depth(depth + 1);
            max_depth = depth + 1;
            solved = 1;
            return 1;
        }
        ++next->actual;
        bits_num = mpz_bitsize(mpq_numref(next->q));
        if (!next->best_bits_num || next->best_bits_num >= bits_num) {
            next->best_bits_num = bits_num;
            bits_den = mpz_bitsize(mpq_denref(next->q));
            if (!next->best_bits_den || next->best_bits_den > bits_den)
                next->best_bits_den = bits_den;
        }

        /* and recurse */
        if (try_depth(depth + 1))
            locally_solved = 1;
    }
    /* not reached */
}

bool search_depth(void) {
    bool result;
    mpq_set(stack[0].q, rone);
    stack[0].count = 0;
    stack[0].actual = 1;
    result = try_depth(0);
    report_final();
    return result;
}
