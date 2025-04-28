#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "int.h"
#include "frag.h"

int fdi, fdor, fdop1, fdop2;
bool seen0;

int record_mark = 0x7265636d;
int file_mark = 0x66696c65;

/* xdata/gnn-nn.rnnnn.pnn */
char path[23 + 10];
char *resolve_path(uint ri) {
    snprintf(path, sizeof(path), "xdata/g%u-%u.r%u", na, nb, ri);
    return &path[0];
}

char *resolved_path(uint ri, uint pi) {
    snprintf(path, sizeof(path), "xdata/g%u-%u.r%u.p%u", na, nb, ri, pi);
    return &path[0];
}

void resolve_open(uint ri) {
    if (ri == 0) {
        fdi = -1;
        seen0 = 0;
    } else {
        fdi = open(resolve_path(ri), O_RDONLY);
        if (fdi < 0)
            fail("Error opening %s for read: %s\n",
                    resolve_path(ri), strerror(errno));
    }
    fdor = open(resolve_path(ri + 1),
            O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fdor < 0)
        fail("Error opening %s for write: %s\n",
                resolve_path(ri + 1), strerror(errno));
    fdop1 = open(resolved_path(ri, resolves[ri].pi),
            O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fdop1 < 0)
        fail("Error opening %s for write: %s\n",
                resolved_path(ri, resolves[ri].pi), strerror(errno));
    fdop2 = open(resolved_path(ri, resolves[ri].pj),
            O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fdop2 < 0)
        fail("Error opening %s for write: %s\n",
                resolved_path(ri, resolves[ri].pj), strerror(errno));
}

void resolve_close(uint ri) {
    ssize_t chars;

    if (fdi >= 0)
        close(fdi);
    chars = write(fdor, &file_mark, sizeof(file_mark));
    if (chars != sizeof(file_mark))
        fail("write error, wrote %d bytes of %d (%s)\n",
                chars, sizeof(file_mark), strerror(errno));
    close(fdor);
    chars = write(fdop1, &file_mark, sizeof(file_mark));
    if (chars != sizeof(file_mark))
        fail("write error, wrote %d bytes of %d (%s)\n",
                chars, sizeof(file_mark), strerror(errno));
    close(fdop1);
    chars = write(fdop2, &file_mark, sizeof(file_mark));
    if (chars != sizeof(file_mark))
        fail("write error, wrote %d bytes of %d (%s)\n",
                chars, sizeof(file_mark), strerror(errno));
    close(fdop2);
}

void integrate_open(uint ri, uint pi) {
    fdi = open(resolved_path(ri, pi), O_RDONLY);
    if (fdi < 0)
        fail("Error opening %s for read: %s\n",
                resolved_path(ri, pi), strerror(errno));
}

void integrate_close(uint ri, uint pi) {
    close(fdi);
}

bool read_frag(fid_t fi) {
    if (fdi < 0) {
        if (seen0)
            return 0;
        seen0 = 1;
        frag_p(fi)->ps = all_paths();
        for (uint vi = 1; vi <= nv; ++vi) {
            limitp_dup(range_low(frag_range(fi, vi)), LIM0);
            limitp_dup(range_high(frag_range(fi, vi)), LIM1);
        }
        return 1;
    } else {
        char buf[frag_size() + sizeof(record_mark)];
        ssize_t chars = read(fdi, buf, sizeof(buf));
        if (chars == sizeof(buf)) {
            if (memcmp(&buf[frag_size()], &record_mark, sizeof(record_mark)))
                fail("bad record_mark, possible corruption\n");
            memcpy(frag_p(fi), buf, frag_size());
            return 1;
        } else if (chars == sizeof(file_mark)
            && memcmp(&buf[0], &file_mark, sizeof(file_mark)) == 0
        ) {
            return 0;
        } else {
            fail("read error: read %d of %d bytes (%s)\n",
                    chars, frag_size() + sizeof(record_mark), strerror(errno));
            exit(1);
        }
    }
}

static inline bool unique_set(pathset_t ps) {
    return ((ps & (ps - 1)) == 0) ? 1 : 0;
}

void write_frag(fid_t fi, pathset_t ps) {
    int fdo = fdor;
    if (ps) {
        pathset_t fps = frag_ps(fi) & ~ps;
        frag_ps_set(fi, fps);
        if (unique_set(fps))
            fdo = (fps < ps) ? fdop1 : fdop2;
    }
    ssize_t chars = write(fdo, frag_p(fi), frag_size());
    if (chars != frag_size())
        fail("write error, wrote %d bytes of %d to %d (%s)\n",
                chars, frag_size(), fdo, strerror(errno));
    chars = write(fdo, &record_mark, sizeof(record_mark));
    if (chars != sizeof(record_mark))
        fail("write error, wrote %d bytes of %d (%s)\n",
                chars, sizeof(record_mark), strerror(errno));
    return;
}

