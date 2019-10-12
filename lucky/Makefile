#GCC=LD_LIBRARY_PATH="/opt/mpfr-2.4.1/lib:/opt/mpc-0.8.1/lib" /opt/gcc-4.5.0/bin/gcc
INCS=
LIBS=
LIBOPTS=
CC_OPT=-O6

PROG=cgen_lucky dcgen_lucky vsieve dvsieve

all: ${PROG}
	true

.SUFFIXES: .c .oo .do

.c.oo:
	${CC} -g ${CC_OPT} ${INCS} -c -o $*.oo $*.c

.c.do:
	${CC} -g -DDEBUG ${INCS} -c -o $*.do $*.c

${PROG}:
	${CC} -o $@ ${CC_OPT} $^ ${LIBOPTS} ${LIBS}

cgen_lucky.do cgen_lucky.oo: cgen_lucky.c
vsieve.do vsieve.oo: vsieve.c
cgen_lucky: cgen_lucky.oo
dcgen_lucky: cgen_lucky.do
vsieve: vsieve.oo
dvsieve: vsieve.do
