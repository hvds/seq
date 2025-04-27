CINT = int.c path.c source.c lincom.c limit.c frag.c calc.c diag.c num.c
HINT = int.h path.h source.h lincom.h limit.h frag.h calc.h diag.h num.h

CC_OPT = -O2
CC_DOPT = -O0 -DDEBUG

int: Makefile ${CINT} ${HINT}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} ${CINT} -lgmp -lrt

dint: Makefile ${CINT} ${HINT}
	gcc -o $@ -g ${CC_DOPT} ${DEFINES} ${CINT} -lgmp -lrt
