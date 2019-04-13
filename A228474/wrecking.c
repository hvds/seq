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
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
/* See http://www.pixelbeat.org/programming/gcc/static_assert.html if this
   doesn't give you static_assert(). */
#define __USE_ISOC11
#include <assert.h>

/* We'll spill the bitvector to this file, make sure it is on a filesystem
   with plenty of room (probably a few TB).
   Default path is '/tmp/wrecking_map-<n>-<k>'.
*/
char* map_path;

/* Type used for the running value; must be signed */
typedef long long INT;

/* RANGE_BITS can go up to 60 without reengineering; beyond that, we probably
   want to ditch all_chunks[], move lru_gen into mapped[] and just scan the
   mapped list each time.
*/
#define RANGE_BITS 48
#define RANGE_MAX ((INT)1 << (RANGE_BITS - 1))
#define RANGE_MIN (-RANGE_MAX)
static_assert(RANGE_BITS <= sizeof(INT) * 8,
        "RANGE_BITS too large for INT type");

/* Values seen are captured in a boolean bitvector notionally running from
   RANGE_MIN+1 to RANGE_MAX-1, actually stored in NCHUNKS chunks
   of CHUNK_SIZE bytes. We keep at most LOCAL_CHUNKS of those chunks in
   memory via an LRU cache, spilling the rest to disk into the file map_path
   declared above.
*/
#define CHUNK_BITS 30
#define CHUNK_SIZE (1 << CHUNK_BITS)
#define BITS_PER_CHUNK (1 << (CHUNK_BITS + 3)
#define NCHUNKS (1 << (RANGE_BITS - CHUNK_BITS))
#define LOCAL_CHUNKS 16
static_assert(RANGE_BITS - CHUNK_BITS <= 30,
        "RANGE_BITS too large, or CHUNK_BITS too small: NCHUNKS will overflow");
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
all_chunk_info *all_chunks;

typedef struct {
    int which_mapped;
    char* chunk;
} mapped_chunk_info;
mapped_chunk_info mapped[LOCAL_CHUNKS];


/* Called at start of run */
void init_chunks(void) {
    if (all_chunks == NULL) {
        all_chunks = (all_chunk_info *)calloc(NCHUNKS, sizeof(all_chunk_info));
        if (all_chunks == NULL) {
            fprintf(stderr,
                "Not enough memory, need room for %d chunks of %lu bytes\n",
                NCHUNKS, sizeof(all_chunk_info)
            );
            exit(1);
        }
    }
    memset((void*)mapped, 0, sizeof(mapped));
    map_fd = open(map_path, O_RDWR | O_CREAT | O_TRUNC, 0666);
    generation = 0;
    max_mapped = -1;
    return;
}

/* Called at end of run */
void clear_chunks(void) {
    int i;

    for (i = 0; i < LOCAL_CHUNKS; ++i) {
        if (mapped[i].chunk) {
            munmap((void*)mapped[i].chunk, CHUNK_SIZE);
        }
    }
    close(map_fd);
    memset((void*)all_chunks, 0, NCHUNKS * sizeof(all_chunk_info));
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
            fprintf(stderr, "Could not extend file to size %lld\n", (INT)size);
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

INT A228474(int n0, INT limit) {
    INT d = 0;
    INT n = (INT)n0;
    INT next;

    while (n) {
        mark(n);

        ++d;
        if (limit && d > limit) {
            return (INT)-1;
        }

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
    return d;
}

int main(int argc, char** argv) {
    int n, k = 1;
    char path[1+3+1+12+1+11+1+11+1];    /* /tmp/wrecking_map-nnn-kkk\0 */
    INT result, limit = 0;
    int arg = 1;

    while (arg < argc && argv[arg][0] == '-') {
        if (strcmp("--", argv[arg]) == 0) {
            ++arg;
            break;
        } else if (strcmp("-L", argv[arg]) == 0) {
            limit = atoll(argv[++arg]);
            ++arg;
        } else if (strcmp("-F", argv[arg]) == 0) {
            map_path = argv[++arg];
            ++arg;
        } else {
            fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
            return 1;
        }
    }

    if (argc - arg < 1 || argc - arg > 2) {
        fprintf(stderr,
            "Usage: %s [ -L <limit> ] [ -F <filepath> ] <n> [ <k> ]\n",
            argv[0]
        );
        return 1;
    }

    n = atoi(argv[arg++]);
    if (arg < argc) {
        k = atoi(argv[arg]);
    }

    if (map_path == (char*)NULL) {
        snprintf(path, sizeof(path), "/tmp/wrecking_map-%u-%u", n, k);
        map_path = path;
    }

    while (k) {
        init_chunks();
        result = A228474(n, limit);
        if (result >= 0) {
            printf("a(%d) = %lld\n", n, result);
        } else {
            printf("a(%d) > %lld\n", n, limit);
        }
        clear_chunks();
        ++n;
        --k;
    }

    if (unlink(map_path)) {
        fprintf(stderr, "Warning, could not remove map file '%s'\n", map_path);
    }
    if (all_chunks != NULL) {
        free(all_chunks);
    }
    return 0;
}
