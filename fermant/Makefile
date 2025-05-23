CFILES = path.c source.c lincom.c limit.c frag.c calc.c diag.c num.c
HFILES = int.h path.h source.h lincom.h limit.h frag.h calc.h diag.h num.h

CC_OPT = -O2
CC_DOPT = -O0 -DDEBUG

int: Makefile int.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} int.c ${CFILES} -lgmp -lrt

dint: Makefile int.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_DOPT} ${DEFINES} int.c ${CFILES} -lgmp -lrt

check: Makefile check.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} check.c ${CFILES} -lgmp -lrt

dcheck: Makefile check.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_DOPT} ${DEFINES} check.c ${CFILES} -lgmp -lrt

showfrag: Makefile showfrag.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_DOPT} ${DEFINES} showfrag.c ${CFILES} -lgmp -lrt
