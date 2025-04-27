#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "int.h"

/* xdata/g3-3.15 */
char pfilepath[ 7 + 2 + 1 + 2 + 1 + 4 + 1 ];
char *filepath(uint ri) {
    snprintf(pfilepath, sizeof(pfilepath), "xdata/g%u-%u.%u",
            na, nb, ri);
    return &pfilepath[0];
}

int resolve_reader(uint ri) {
    if (ri == 0)
        return 0;
    int fdi = open(filepath(ri), O_RDONLY);
    if (fdi < 0)
        fail("Error opening %s for read: %s\n", filepath(ri), strerror(errno));
    return fdi;
}

int resolve_writer(uint ri) {
    int fdo = open(filepath(ri + 1), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fdo < 0)
        fail("Error opening %s for write: %s\n",
                filepath(ri + 1), strerror(errno));
    return fdo;
}

