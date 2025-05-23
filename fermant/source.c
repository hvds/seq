#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "int.h"
#include "frag.h"
#include "calc.h"

int fdi, fdor, fdop1, fdop2;
bool seen0;
uint write_count = 0;

typedef struct {
    bool valid;
    uint rid;
    off_t ofdi, ofdor, ofdop1, ofdop2;
} resolve_cp_t;
typedef struct {
    bool valid;
    uint rid, pid;
    off_t ofdi;
} integrate_cp_t;
resolve_cp_t rcp;
integrate_cp_t icp;
uint rid, pid;

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

int ropen(char *path, bool recover, off_t off) {
    int fd = open(path, O_RDONLY);
    if (fd < 0)
        fail("Error opening %s for read: %s\n", path, strerror(errno));
    if (recover) {
        off_t got = lseek(fd, off, SEEK_SET);
        if (got != off)
            fail("Error recovering read at %ld (%ld) of %s: %s\n",
                    (long)off, (long)got, path, strerror(errno));
    }
    return fd;
}

int wopen(char *path, bool recover, off_t off) {
    int fd;
    if (debug_suppress_write)
        return -1;
    if (!recover || off == 0) {
        fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0666);
        if (fd < 0)
            fail("Error opening %s for write: %s\n", path, strerror(errno));
        return fd;
    }
    fd = open(path, O_RDONLY);
    off_t got = lseek(fd, off - sizeof(record_mark), SEEK_SET);
    if (got != off - sizeof(record_mark))
        fail("Error seeking to record mark in %s: %s\n", path, strerror(errno));
    char buf[sizeof(record_mark)];
    if (read(fd, buf, sizeof(record_mark)) != sizeof(record_mark))
        fail("Error reading record mark in %s: %s\n", path, strerror(errno));
    if (memcmp(buf, &record_mark, sizeof(record_mark)) != 0)
        fail("Record mark mismatch in %s\n", path);
    close(fd);
    fd = open(path, O_WRONLY);
    if (fd < 0)
        fail("panic: could not reopen %s: %s\n", path, strerror(errno));
    got = lseek(fd, off, SEEK_SET);
    if (got != off)
        fail("panic: could not re-seek to %s of %s: %s\n", path, strerror(errno));
    return fd;
}

bool resolve_open(uint ri) {
    if (icp.valid || (rcp.valid && ri != rcp.rid))
        return 0;
    rid = ri;
    if (ri == 0) {
        fdi = -1;
        seen0 = 0;
        if (rcp.valid && rcp.ofdi > 0)
            seen0 = 1;
    } else {
        fdi = ropen(resolve_path(ri), rcp.valid, rcp.ofdi);
    }
    fdor = wopen(resolve_path(ri + 1), rcp.valid, rcp.ofdor);
    fdop1 = wopen(resolved_path(ri, resolves[ri].pi), rcp.valid, rcp.ofdop1);
    fdop2 = wopen(resolved_path(ri, resolves[ri].pj), rcp.valid, rcp.ofdop2);
    rcp.valid = 0;
    return 1;
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

bool integrate_open(uint ri, uint pi) {
    if (icp.valid && (ri != icp.rid || pi != icp.pid))
        return 0;
    rid = ri;
    pid = pi;
    fdi = ropen(resolved_path(ri, pi), icp.valid, icp.ofdi);
    icp.valid = 0;
    return 1;
}

void integrate_close(uint ri, uint pi) {
    close(fdi);
}

bool read_frag(fid_t fi) {
    if (fdi < 0) {
        if (seen0)
            return 0;
        seen0 = 1;
        frag_p(fi)->id = nextfragid++;
        frag_p(fi)->parent = (fid_t)-1;
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

void do_sync(int fd) {
    if (fsync(fd) != 0)
        fail("error syncing fd %d: %s", fd, strerror(errno));
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
    if (!debug_suppress_write) {
        ssize_t chars = write(fdo, frag_p(fi), frag_size());
        if (chars != frag_size())
            fail("write error, wrote %d bytes of %d to %d (%s)\n",
                    chars, frag_size(), fdo, strerror(errno));
        chars = write(fdo, &record_mark, sizeof(record_mark));
        if (chars != sizeof(record_mark))
            fail("write error, wrote %d bytes of %d (%s)\n",
                    chars, sizeof(record_mark), strerror(errno));
        ++write_count;
        if (sync_count && (write_count % sync_count) == 0) {
            if (sync_stderr)
                do_sync(STDERR_FILENO);
            do_sync(fdor);
            do_sync(fdop1);
            do_sync(fdop2);
        }
    }
    return;
}

void resolve_checkpoint(void) {
    report("305 resolve %u %u %ld %ld %ld %ld (%.2fs)\n", rid, nextfragid,
        (long)lseek(fdi, 0, SEEK_CUR),
        (long)lseek(fdor, 0, SEEK_CUR),
        (long)lseek(fdop1, 0, SEEK_CUR),
        (long)lseek(fdop2, 0, SEEK_CUR),
        seconds((double)utime())
    );
}

void integrate_checkpoint(void) {
    char buf[4096];
    uint off = 0;
    for (uint i = 0; i < npaths; ++i)
        off += gmp_snprintf(&buf[off], sizeof(buf) - off, " %Qd",
                path_total[i]);
    /* if we don't fit, mark it so we don't silently recover wrong values */
    if (off >= sizeof(buf) - 1)
        buf[0] = 'x';
    report("305 integrate %u %u %u %ld%s (%.2fs)\n",
            rid, pid, nextfragid, (long)lseek(fdi, 0, SEEK_CUR), buf,
            seconds((double)utime()));
}

void set_resolve_cp(uint ri, off_t oi, off_t oor, off_t oop1, off_t oop2) {
    rcp.rid = ri;
    rcp.ofdi = oi;
    rcp.ofdor = oor;
    rcp.ofdop1 = oop1;
    rcp.ofdop2 = oop2;
    rcp.valid = 1;
}

void set_integrate_cp(uint ri, uint pi, off_t oi) {
    icp.rid = ri;
    icp.pid = pi;
    icp.ofdi = oi;
    icp.valid = 1;
}

long page_size; /* expecting 4096 */
#define TARGET_SIZE (64L * 1024L * 1024L)
void *mmfrag = NULL;
struct stat statbuf;
off_t chunk_size, file_size, indent, map_size;
uint records, this_record;

void mmfrag_init(void) {
    page_size = sysconf(_SC_PAGE_SIZE);
    if (TARGET_SIZE % page_size)
        fail("target size %lu is not a multiple of page size %ld\n",
                TARGET_SIZE, page_size);
    chunk_size = frag_size() + sizeof(record_mark);
    fdi = -1;
    mmfrag = NULL;
}

void mmfrag_done(void) {
    if (mmfrag)
        munmap(mmfrag, map_size);
    if (fdi >= 0)
        close(fdi);
}

uint mmfrag_open(bool resolved, uint ri, uint pi) {
    if (fdi >= 0)
        close(fdi);
    rid = ri;
    pid = pi;
    char *path = resolved ? resolved_path(rid, pid) : resolve_path(rid);
    fdi = ropen(path, 0, 0);
    if (fstat(fdi, &statbuf) != 0)
        fail("could not stat %s: %s\n", path, strerror(errno));
    file_size = statbuf.st_size - sizeof(file_mark);
    return file_size / chunk_size;
}

frag_t *mmfrag_get(uint rec) {
    off_t off = (off_t)rec * chunk_size;
    if (mmfrag && (off >= indent) && (off + chunk_size < indent + map_size))
        return (frag_t *)add_p(mmfrag, off - indent);
    if (mmfrag)
        munmap(mmfrag, map_size);
    indent = off - (off % page_size);
    map_size = file_size - indent;
    if (map_size > TARGET_SIZE)
        map_size = TARGET_SIZE;
    /* I assume that if there's any difference at all, the kernel can do
     * less work for MAP_SHARED on a read-only extent than for MAP_PRIVATE */
    mmfrag = mmap(NULL, map_size, PROT_READ, MAP_SHARED, fdi, indent);
    if (mmfrag == (void *)-1)
        fail("mmap error for size %ld from %ld for file(%u, %u): %s\n",
                map_size, indent, rid, pid, strerror(errno));
    return (frag_t *)add_p(mmfrag, off - indent);
}
