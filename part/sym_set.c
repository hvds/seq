#define IS_SYM_SET_C
#include "sym_set.h"
#include "symmetries.h"
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

char* sym_set_arena;
uint sym_set_count;

uint sym_set_size;
uint sym_set_count;
sym_set_t** sets;

void setup_sym_set(void) {
}

void teardown_sym_set(void) {
	free(sym_set_arena);
}

sym_set_t* symset_new(void) {
	return (sym_set_t*)malloc(sizeof(sym_set_t) + sym_count * sizeof(int));
}
sym_set_t* symset_resize(sym_set_t* ss) {
	return (sym_set_t*)realloc(ss, SSSIZE(ss));
}
void symset_delete(sym_set_t* ss) {
	free(ss);
}

void load_sym_sets(FILE* f, uint size) {
	uint setsize = size - sizeof(uint);
	uint read;

	if (size < sizeof(uint)) {
		fprintf(stderr, "load_sym_sets: expected size %u > %u\n",
				size, sizeof(uint));
		fclose(f);
		exit(-1);
	}
	read = fread(&sym_set_count, sizeof(uint), 1, f);
	if (read < 1) {
		fprintf(stderr,
			"load_sym_sets could not read sym_set_count: %s (%d)\n",
			strerror(errno), errno
		);
		fclose(f);
		exit(-1);
	}
	sym_set_arena = (char*)malloc(setsize);
	read = fread(sym_set_arena, sizeof(char), setsize, f);
	if (read < setsize) {
		fprintf(stderr,
			"load_sym_sets error reading sym_sets (want %u, got %u): %s (%d)\n",
			setsize, read, strerror(errno), errno
		);
		fclose(f);
		exit(-1);
	}
#ifndef NDEBUG
	{
		uint offset;
		uint count = 0;
		sym_set_t* ss;
		for (offset = 0; offset < setsize; offset += SSSIZE(ss)) {
			ss = SSO(offset);
			++count;
		}
		assert(count == sym_set_count);
		assert(offset == setsize);
	}
#endif
	fprintf(stderr, "load_sym_sets: loaded %u sym_sets, size %u\n",
			sym_set_count, setsize);
	return;
}
