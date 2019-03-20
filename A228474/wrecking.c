/*
  A228474 - Wrecking ball, a Recaman-like sequence

  Starting at n, a(n) is the number of steps required to reach zero.
  On the k-th step (k=1,2,3,...), move a distance of k in the direction
  of zero. If the number landed on has been landed on before, move a distance
  of k away from zero instead. Set a(n) = -1 if n never reaches 0.
*/

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
/* See http://www.pixelbeat.org/programming/gcc/static_assert.html if this
   doesn't give you static_assert(). */
#define __USE_ISOC11
#include <assert.h>

/* We'll spill the bitvector to this file, make sure it is on a filesystem
   with plenty of room (probably a few TB).
*/
const char* map_path = "/tmp/wrecking_map";

/* Type used for the running value; must be signed */
typedef long long INT;

/* RANGE_BITS can probably go up to around 61 without reengineering, at the
   cost of more NCHUNKS. */
#define RANGE_BITS 48
#define RANGE_MAX ((INT)1 << (RANGE_BITS - 1))
#define RANGE_MIN (-RANGE_MAX)
static_assert(RANGE_BITS <= sizeof(INT) * 8,
        "RANGE_BITS too large for INT type");

/* Values seen are captured in a boolean bitvector notionally running from
   min_INT to max_INT (excluding zero), actually stored in NCHUNKS chunks
   of CHUNK_SIZE bytes. We keep at most LOCAL_CHUNKS of those chunks in
   memory via an LRU cache, spilling the rest to disk into the file map_path
   declared below.
*/
#define CHUNK_BITS 30
#define CHUNK_SIZE (1 << CHUNK_BITS)
#define BITS_PER_CHUNK (1 << (CHUNK_BITS + 3)
#define NCHUNKS (1 << (RANGE_BITS - CHUNK_BITS))
#define LOCAL_CHUNKS 16
static_assert(CHUNK_BITS < sizeof(int) * 8,
        "CHUNK_BITS too large to use 'int' to index chunks");
static_assert(RANGE_BITS - CHUNK_BITS < sizeof(int) * 8,
        "RANGE_BITS - CHUNK_BITS too large to index all_chunks[NCHUNKS]");


/* Each access increments generation, and sets lru_gen to it. */
INT generation = 0;
int max_mapped = -1;
int map_fd;

typedef struct {
    INT lru_gen;
    char* chunk;
} all_chunk_info;
all_chunk_info all_chunks[NCHUNKS];

typedef struct {
    int which_mapped;
    char* chunk;
} mapped_chunk_info;
mapped_chunk_info mapped[LOCAL_CHUNKS];


/* Called at start of run */
void init_chunks(void) {
    /* Note: you'd need to memfill all_chunks[] and mapped[] with zeros
       to reuse them within the same process. */
    map_fd = open(map_path, O_RDWR | O_CREAT | O_TRUNC, 0666);
    return;
}

/* Called at end of run */
void clear_chunks(void) {
    int i;

    /* Run has completed, so no point writing stuff to disk now. */
#if 0
    for (i = 0; i < LOCAL_CHUNKS; ++i) {
        if (mapped[i].chunk) {
            munmap((void*)mapped[i].chunk, CHUNK_SIZE);
        }
    }
#endif
    close(map_fd);
    return;
}

/* Load the specified chunk, return a pointer to it */
char* map_chunk(int index) {
    INT minru = -1, thisru;
    int i, cache_index;
    char* chunk;

    /* choose a local chunk to use */
    for (i = 0; i < LOCAL_CHUNKS; ++i) {
        if (!mapped[i].chunk) {
            cache_index = i;
            break;
        }
        thisru = all_chunks[ mapped[i].which_mapped ].lru_gen;
        if (minru == -1 || thisru < minru) {
            cache_index = i;
            minru = thisru;
        }
    }

    /* clean up any previous use of the local chunk */
    if (mapped[cache_index].chunk) {
        munmap((void*)mapped[cache_index].chunk, CHUNK_SIZE);
        all_chunks[ mapped[cache_index].which_mapped ].chunk = NULL;
    }

#ifdef DEBUG
    fprintf(stderr, "map chunk %d to local %d", index, cache_index);
    if (mapped[cache_index].chunk) {
        fprintf(stderr, " (was %d)\n", mapped[cache_index].which_mapped);
    } else {
        fprintf(stderr, " (was free)\n");
    }
#endif

    /* extend the file if needed */
    if (index > max_mapped) {
        off_t size = ((off_t)index + 1) << CHUNK_BITS;
        if (ftruncate(map_fd, size) != 0) {
            fprintf(stderr, "Could not extend file to size %lld", (INT)size);
            exit(1);
        }
        max_mapped = index;
    }

    /* now map the chunk */
    chunk = (char*)mmap(
        NULL, CHUNK_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED,
        map_fd, (off_t)index << CHUNK_BITS
    );
    if (chunk == MAP_FAILED) {
        fprintf(stderr, "Could not map chunk %d\n", index);
        exit(1);
    }

    mapped[cache_index].which_mapped = index;
    mapped[cache_index].chunk = chunk;
    all_chunks[index].chunk = chunk;
    return chunk;
}

/* Return a pointer to the chunk containing the specified value */
char* chunk_for(INT which) {
    /* chunks alternate +ve and -ve */
    int pair_index = llabs(which) >> (CHUNK_BITS + 3);
    int index = (pair_index << 1) + ((which < 0) ? 1 : 0);
    char* chunk = all_chunks[index].chunk;

    if (!chunk)
        chunk = map_chunk(index);
    all_chunks[index].lru_gen = ++generation;
    return chunk;
}

/* Return the offset within a chunk for the specified value */
int offset_for(INT which) {
    return (llabs(which) >> 3) & ((1 << CHUNK_BITS) - 1);
}

/* Return the bit value within the offset for the specified value */
int bit_for(INT which) {
    return 1 << (llabs(which) & 7);
}

/* Boolean, TRUE if the specified value is marked as seen in the bitvector */
int marked(INT which) {
    char* chunk = chunk_for(which);
    return (chunk[offset_for(which)] & bit_for(which)) ? 1 : 0;
}

/* Mark the specified value as seen in the bitvector */
void mark(INT which) {
    char* chunk = chunk_for(which);
    chunk[offset_for(which)] |= bit_for(which);
    return;
}

void A228474(int n0) {
    INT d = 0;
    INT n = (INT)n0;
    INT next;

    while (n) {
        mark(n);
        ++d;
        /* report progress from time to time */
        if ((d & 0x07ffffff) == 0)
            fprintf(stderr, "d=%lld, max_mapped=%d\n", d, max_mapped);

        if (n > 0) {
            next = n - d;
            if (marked(next)) {
                next = n + d;
                if (next < 0 || next >= RANGE_MAX) {
                    fprintf(stderr, "Overflow for %d at %lld: %lld -> %lld\n",
                            n0, d, n, next);
                    exit(1);
                }
            }
        } else {
            next = n + d;
            if (marked(next)) {
                next = n - d;
                if (next > 0 || next <= RANGE_MIN) {
                    fprintf(stderr, "Overflow for %d at %lld: %lld -> %lld\n",
                            n0, d, n, next);
                    exit(1);
                }
            }
        }
        n = next;
    }
    printf("a(%d) = %lld\n", n0, d);
    return;
}

int main(int argc, char** argv) {
    int n;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    n = atoi(argv[1]);

    init_chunks();
    A228474(n);
    clear_chunks();
    return 0;
}
