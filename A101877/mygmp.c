#include "mygmp.h"

inline void QINIT(mpq_t* q, char* format, ...) {
	mpq_init(*q);
	DEBUG_PRINT(q, "qinit", format);
}

inline void QCLEAR(mpq_t* q, char* format, ...) {
	mpq_clear(*q);
	DEBUG_PRINT(q, "qclear", format);
}

inline void ZINIT(mpz_t* z, char* format, ...) {
	mpz_init(*z);
	DEBUG_PRINT(z, "zinit", format);
}

inline void ZCLEAR(mpz_t* z, char* format, ...) {
	mpz_clear(*z);
	DEBUG_PRINT(z, "zclear", format);
}
