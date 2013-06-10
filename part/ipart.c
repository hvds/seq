#include "ipart.h"
#include "part.h"

ipart_t* ip;

void setup_ipart(void) {
	ip = (ipart_t*)malloc(sizeof(ipart_t) + NODES * sizeof(unsigned int));
}

void teardown_ipart(void) {
	free(ip);
}

ipart_t* ipart_first(unsigned int major) {
	unsigned int remain = NODES;
	ip->count = 0;
	while (remain >= major) {
		ip->parts[ip->count++] = major;
		remain -= major;
	}
	if (remain > 0)
		ip->parts[ip->count++] = remain;
	return ip;
}

ipart_t* ipart_next(void) {
	unsigned int last = ip->count - 1;
	unsigned int remain = 0;
	unsigned int chunk;
	while (last > 0) {
		++remain;
		if (ip->parts[last] > 1)
			break;
		--last;
	}
	if (last == 0)
		return (ipart_t*)NULL;

	ip->count = last;
	chunk = --ip->parts[last++];
	while (remain >= chunk) {
		ip->parts[last++] = chunk;
		remain -= chunk;
	}
	if (remain > 0)
		ip->parts[last++] = remain;
	ip->count = last;
	return ip;
}

void fprint_ipart(FILE* stream, ipart_t* ipart) {
	unsigned int u;
	fprintf(stream, "[%u] %u", ipart->count, ipart->parts[0]);
	for (u = 1; u < ipart->count; ++u)
		fprintf(stream, ".%u", ipart->parts[u]);
}
