#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "array.h"
#include "seen.h"

/*
 * Support for a "seen" pipeline for GMP big integers.
 *
 * Small numbers are dealt with by a bit vector; larger numbers are
 * recorded in a series of hash tables, grouped by the size of the GMP
 * integers in limbs.
 *
 * TODO: use 64-bit ints for arena size, arena offsets, buckets; reduce
 * memory footprint a little by not storing the hash value - it's quick
 * enough to recalculate, especially once it's 64 bits.
 */

typedef unsigned int uint;

/*
 * Each hash has a 2^x array of buckets, each pointing to a linked list
 * of values in that bucket terminated by NO_NEXT.
 * When the number of elements in the hash is the same as the number of
 * buckets, we double the number of buckets.
 */
typedef struct {
    uint *buckets;
    uint bucket_count;
    uint count;
} hash_t;

/*
 * Hash elements are stored in the arena, and referenced by an offset in
 * uints into the arena.
 */
typedef struct {
    uint next;
    uint hash;
    mp_limb_t limbs[0];
} entry_t;

/* TODO: use 0 as the marker instead */
#define NO_NEXT ((uint)~1)


/*
 * 512MB vector handles the first 2^32 numbers
 */
const uint vec_range = 0xffffffff;
const uint vec_size = 1 << 29;
char smallvec[1 << 29];


/*
 * Arena of seen values larger than vec_range, each holding the size and limbs
 * of the value, the hashed value, and a 'next' pointer for other values in
 * the same hash bucket
 */
array_t arena;

/*
 * Array of hashes, one for each limb size
 */
array_t hashes;


/*
 * Arena support
 */
inline entry_t *arena_off(uint u) {
    uint *arena_0 = (uint *)arena.array;
    return (entry_t *)&(arena_0[u]);
}

inline void resize_arena(uint size) {
    resize_array(&arena, size);
}


/*
 * Hash support
 */
inline hash_t *hash(uint limb_size) {
    return (hash_t *)array_element(&hashes, limb_size);
}

inline void init_hashes(uint limbs) {
    uint i, j;
    if (limbs >= hashes.count) {
        resize_array(&hashes, limbs + 1);
        for (i = hashes.count; i < limbs + 1; ++i) {
            hash_t *h = hash(i);
            h->bucket_count = 256;
            h->buckets = (uint *)malloc(sizeof(uint) * h->bucket_count);
            for (j = 0; j < h->bucket_count; ++j)
                h->buckets[j] = NO_NEXT;
            h->count = 0;
        }
        hashes.count = limbs + 1;
    }
}

void free_hash(uint limbs) {
    hash_t *h = hash(limbs);
    if (h->buckets)
        free(h->buckets);
}

/*
 * Double the number of buckets for this hash, reassign all the values.
 */
void split_hash(uint limbs) {
    hash_t *h = hash(limbs);
    uint bucket_count = h->bucket_count;
    uint i;
    /* keep some stats to check hash quality */
    ulong stat_n = 0;   /* number of non-empty buckets */
    ulong stat_n2 = 0;  /* sum of squares of length of linked lists */
    ulong len;

    h->bucket_count <<= 1;
    h->buckets = (uint *)realloc(h->buckets, sizeof(uint) * h->bucket_count);
    for (i = 0 ; i < bucket_count; ++i) {
        uint low = NO_NEXT, high = NO_NEXT;
        uint chain = h->buckets[i], next;
        len = 0;
        while (chain != NO_NEXT) {
            ++len;
            entry_t *he = arena_off(chain);
            next = he->next;
            if ((he->hash & (h->bucket_count - 1)) < bucket_count) {
                he->next = low;
                low = chain;
            } else {
                he->next = high;
                high = chain;
            }
            chain = next;
        }
        h->buckets[i] = low;
        h->buckets[i + bucket_count] = high;
        stat_n += len ? 1 : 0;
        stat_n2 += len * len;
    }
    printf("split hash[%d] to %d buckets (sum_n %ld, sum_n2 %ld)\n",
            limbs, h->bucket_count, stat_n, stat_n2);
}

/* Show stats on all hashes */
void dump_seen(void) {
    uint i, first = 1;
    printf("seen [");
    for (i = 1; i < hashes.count; ++i) {
        if (!hash(i)->buckets)
            continue;
        printf("%s%u:%u", first ? "" : " ", i, hash(i)->count);
        first = 0;
    }
    printf("]");
}


/*
 * Hash entry support
 */

/*
 * Add a new hash entry for the supplied integer at the top of the arena,
 * and return its offset.
 *
 * The last entry made can be discarded by setting arena.count to
 * the returned offset.
 */
uint make_entry(mpz_t z) {
    uint limbs = mpz_size(z);
    uint size = (sizeof(entry_t) + limbs * sizeof(mp_limb_t)) / sizeof(uint);
    uint offset = arena.count;
    uint i;
    mp_limb_t hv64 = 0; /* hash value */
    entry_t *he;
    mp_limb_t *hz;

    resize_arena(offset + size);
    he = arena_off(offset); /* Note that resize_arena() can invalidate this */
    arena.count += size;

    he->next = NO_NEXT;
    hz = (mp_limb_t *)&(he->limbs);
    for (i = 0; i < limbs; ++i) {
        hz[i] = mpz_getlimbn(z, (mp_size_t)(limbs - i - 1));
        hv64 ^= hz[i];
    }
    he->hash = (hv64 & 0xfffffffful) ^ (hv64 >> 32);
    return offset;
}

/*
 * Returns 1 if the specified hash entry has already been seen, else 0.
 */
int hash_exists(uint limbs, uint offset) {
    hash_t *h = hash(limbs);
    entry_t *he_s, *he_d;
    uint next_off;

    he_s = arena_off(offset);
    next_off = h->buckets[he_s->hash & (h->bucket_count - 1)];
    while (next_off != NO_NEXT) {
        he_d = arena_off(next_off);
        if (0 == memcmp((void *)&(he_s->limbs), (void *)&(he_d->limbs),
                limbs * sizeof(mp_limb_t))) {
            return 1;
        }
        next_off = he_d->next;
    }
    return 0;
}

/*
 * Insert the specified hash entry into the appropriate hash.
 */
void hash_insert(uint limbs, uint offset) {
    hash_t *h = hash(limbs);
    entry_t *he_s = arena_off(offset);
    uint bucket = he_s->hash & (h->bucket_count - 1);

    he_s->next = h->buckets[bucket];
    h->buckets[bucket] = offset;
    ++h->count;
    if (h->count >= h->bucket_count)
        split_hash(limbs);
}

void init_seen(void) {
    init_array(&arena, sizeof(uint), 1 << 24);
    init_array(&hashes, sizeof(hash_t), 10);
    memset(smallvec, 0, sizeof(smallvec));
}

void free_seen(void) {
    uint i;
    free_array(&arena);
    for (i = 0; i < hashes.count; ++i)
        free_hash(i);
    free_array(&hashes);
}

/*
 * Returns TRUE if the given integer has been seen before; else marks it as
 * seen and returns FALSE.
 */
int seen(mpz_t z) {
    size_t limbs = mpz_size(z);
    hash_t *h;
    uint offset;
    entry_t *he;

    /* Handle small integers directly by the bit-vector */
    if (limbs == 1 && mpz_cmp_ui(z, vec_range) < 0) {
        uint n = mpz_get_ui(z);
        uint byte = n >> 3;
        int offset = 1 << (n & 7);
        if (smallvec[byte] & offset)
            return 1;
        smallvec[byte] |= offset;
        return 0;
    }

    init_hashes(limbs);
    h = hash(limbs);
    offset = make_entry(z);
    if (hash_exists(limbs, offset)) {
        /* give back the space */
        arena.count = offset;
        return 1;
    }

    hash_insert(limbs, offset);
    return 0;
}

/*
 * Compare two integers, returning standard -1, 0 or +1.
 */
int limb_cmp(uint limbs, mp_limb_t *left, mp_limb_t *right) {
    uint i;

    for (i = 0; i < limbs; ++i) {
        if (left[i] < right[i])
            return -1;
        if (left[i] > right[i])
            return 1;
    }
    return 0;
}

/*
 * Return the highest number seen, writing it into the supplied mpz_t
 * (which the caller must have initialized).
 */
void max_seen(mpz_t z) {
    uint max_limbs = hashes.count;

    /* Find the hash for the largest limb size we've seen */
    while (max_limbs > 0) {
        --max_limbs;
        if (hash(max_limbs)->buckets)
            break;
    }
    if (max_limbs) {
        /* Find the largest value in this hash */
        uint i, off, first = 1;
        hash_t *h = hash(max_limbs);
        mp_limb_t best[max_limbs];

        for (i = 0; i < h->bucket_count; ++i) {
            off = h->buckets[i];
            while (off != NO_NEXT) {
                entry_t *he = arena_off(off);
                if (first || limb_cmp(max_limbs, best, he->limbs) < 0) {
                    memcpy(best, he->limbs, max_limbs * sizeof(mp_limb_t));
                    first = 0;
                }
                off = he->next;
            }
        }
        mpz_import(z, max_limbs, 1, sizeof(mp_limb_t), 0, 0, (void*)best);
    } else {
        /* nothing in the hashes, check smallvec */
        char* vec = smallvec + vec_size - 1;
        uint byte, offset;

        while (!*vec && vec > smallvec)
            --vec;
        byte = vec - smallvec;
        offset = 0;
        while (*vec > 1) {
            (*vec) >>= 1;
            ++offset;
        }
        mpz_set_ui(z, byte * 8 + offset);
    }
    return;
}
