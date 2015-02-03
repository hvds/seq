#define IS_PIECES_C
#include "pieces.h"
#include "symmetries.h"
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

uint piece_array[NODES + 2];
piece_t* pieces;
uint count_array[NODES + 1];
uint* counts;

void setup_pieces(void) {
}

void teardown_pieces(void) {
	free(pieces);
}

void load_pieces(FILE* f, uint size) {
	uint setsize = size - NODES * sizeof(uint);
	uint read, count[NODES], i, total;

	if (size < NODES * sizeof(uint)) {
		fprintf(stderr, "load_pieces: expected size %u > %u\n",
				size, NODES * sizeof(uint));
		fclose(f);
		exit(-1);
	}
	read = fread(&count, sizeof(uint), NODES, f);
	if (read < NODES) {
		fprintf(stderr,
			"load_pieces error reading counts (want %u, got %u): %s (%d)\n",
			NODES, read, strerror(errno), errno
		);
		fclose(f);
		exit(-1);
	}
	pieces = (piece_t*)malloc(setsize);
	read = fread(pieces, 1, setsize, f);
	if (read < setsize) {
		fprintf(stderr,
			"load_pieces error reading pieces (want %u, got %u): %s (%d)\n",
			setsize, read, strerror(errno), errno
		);
		fclose(f);
		exit(-1);
	}

	total = 0;
	for (i = 0; i < NODES; ++i) {
		piece_array[i + 1] = total;
		total += count[i];
	}
	piece_array[NODES + 1] = total;
	assert(total * sizeof(piece_t) == setsize);
	fprintf(stderr,
		"load_pieces: loaded %u pieces, size %u + %u\n",
		total, NODES * sizeof(uint), setsize
	);
	return;
}

void load_counts(FILE* f, uint size) {
	uint read, i;

	counts = (uint*)malloc(size);
	read = fread(counts, size, 1, f);
	if (read < 1) {
		fprintf(stderr,
			"load_counts error reading counts (want 1 * %u, got %u): %s (%d)\n",
			size, read, strerror(errno), errno
		);
		fclose(f);
		exit(-1);
	}

	count_array[0] = 0;
	for (i = 1; i <= NODES; ++i) {
		count_array[i] = count_array[i - 1] + i * (piece_array[i + 1] - piece_array[i]);
	}
	assert(count_array[NODES] * sizeof(uint) == size);
	fprintf(stderr,
		"load_counts: loaded %u counts, size %u\n",
		count_array[NODES], size
	);
	return;
}
