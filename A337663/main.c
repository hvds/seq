#include <stdlib.h>
#include <stdio.h>

#include "board.h"
#include "group.h"
#include "sym.h"

board_t *init(int n, int freq, char *start_hist) {
    init_sym();
    init_group();
    return init_board(n, freq, start_hist);
}

void finish(void) {
    finish_board();
    finish_group();
    finish_sym();
}

int main(int argc, char** argv) {
    board_t *b;
    int n = 2, freq = 100;

    setvbuf(stdout, (char *)NULL, _IOLBF, 0);

    if (argc > 1) {
        n = atoi(argv[1]);
        if (n >= 9) {
            fprintf(stderr, "Cannot yet calculate n >= 9, need to be able to"
                    " coalesce more than 2 groups simultaneously\n");
            exit(1);
        }
        if (n < 1) {
            fprintf(stderr, "Error, need 1 <= n <= 8\n");
            exit(1);
        }
    }
    if (argc > 2) {
        freq = atoi(argv[2]);
    }

    b = init(n, freq, (argc > 3) ? argv[3] : (char*)NULL);
    try_board(b);
    printf("a(%d) = %d (%lu)\n", n, best_k, board_count);

    finish();
    return 0;
}
