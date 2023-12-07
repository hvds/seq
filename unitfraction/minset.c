#include "unit.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

mpq_t qr, qs;

void init(uint depth) {
    init_unit(depth + 1);
    QINIT(qr);
    QINIT(qs);
}

void done(void) {
    QCLEAR(qr);
    QCLEAR(qs);
    done_unit();
}

/* Find min(| S |) such that sum_{s in S}{1/s} = q.
 * With -m, finds min(| M |) such that sum_{m in M}{1/m} = q.
 */
int main(int argc, char** argv) {
    uint arg = 1;
    uint min_depth = 0;
    uint max_depth = 20;
    bool multi = 0;
    bool square = 0;

    while (arg < argc && argv[arg][0] == '-') {
        char* s = argv[arg++];
        if (strcmp(s, "--") == 0)
            break;
        if (s[1] == 'n') {
            min_depth = atoi(&s[2]);
            continue;
        }
        if (s[1] == 'x') {
            max_depth = atoi(&s[2]);
            continue;
        }
        if (s[1] == 'm') {
            multi = 1;
            continue;
        }
        if (s[1] == 'q') {
            square = 1;
            continue;
        }
        fprintf(stderr, "Unknown option '%s'\n", s);
        argc = -1;  /* force usage message */
        break;
    }

    init(max_depth);
    for (uint i = arg; i < argc; ++i) {
        mpq_set_ui(qr, 1, atoi(argv[i]));
        mpq_add(qs, qs, qr);
    }

    uint best = (square) ? find_square_set(qs, min_depth, max_depth)
            : (multi) ? find_multi(qs) : find_set(qs);
    if (best)
        printf("best %u\n", best);
    else
        printf("no solution found up to depth %d\n", max_depth);

    done();
    return 0;
}
