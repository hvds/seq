GCC=gcc
INCS=-I.
LIBS=-lgmp
LIBOPTS=
CC_OPT=-O6

PROG=g_egyptian dg_egyptian

all: ${PROG}
	true

.SUFFIXES: .c .oo .do

.c.oo:
	${GCC} -g ${CC_OPT} ${INCS} -c -o $*.oo $*.c

.c.do:
	${GCC} -g -DDEBUG ${INCS} -c -o $*.do $*.c

${PROG}:
	${GCC} -o $@ ${CC_OPT} $^ ${LIBOPTS} ${LIBS}

g_egyptian.do g_egyptian.oo: g_egyptian.c
g_egyptian: g_egyptian.oo
dg_egyptian: g_egyptian.do

clean:
	rm -f g_egyptian.do g_egyptian.oo ${PROG}
