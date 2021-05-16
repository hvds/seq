#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

/* Type to use to store counts: any unsigned type up to int should work
 * as long as VAL_MAX is adjusted appropriately.
 */
typedef unsigned char val_t;
#define VAL_MAX 255

/* Number of cubes to sum */
#define COUNT 5

typedef unsigned int full_t;

/* The active function */
#define CUBE(x) x * x * x

full_t *p1;
val_t *px[COUNT - 1];
full_t max, cmax;

void init(void) {
    cmax = CUBE(max);
    p1 = malloc(max * sizeof(full_t));
    for (uint i = 0; i < COUNT - 1; ++i) {
        px[i] = calloc(cmax + 1, sizeof(val_t));
    }
}

void run(void) {
    for (full_t n = 1; n < max; ++n) {
        full_t cn = CUBE(n);
        p1[n - 1] = cn;
        for (full_t k = 0; k < n; ++k) {
            full_t off = p1[k] + cn;
            if (off > cmax)
                break;
            if (px[0][off] < VAL_MAX)
                ++px[0][off];
        }
        full_t tmax = cmax - cn;
        for (uint x = 0; x < COUNT - 2; ++x) {
            val_t *pf = px[x];
            val_t *pt = px[x + 1];
            for (full_t i = 1; i <= tmax; ++i) {
                full_t s = pt[cn + i] + pf[i];
                if (s > VAL_MAX)
                    s = VAL_MAX;
                pt[cn + i] = s;
            }
        }
    }
}

void report(void) {
    char file[5];
    for (uint x = 0; x < COUNT - 1; ++x) {
        snprintf(file, sizeof(file), "x%d", x + 2);
        int fd = open(file, O_CREAT | O_TRUNC | O_WRONLY, 0777);
        write(fd, px[x] + 1, cmax * sizeof(val_t));
        close(fd);
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <max>\n", argv[0]);
        exit(1);
    }
    max = (full_t)atoi(argv[1]);
    init();
    run();
    report();
    return 0;
}
