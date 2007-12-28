LIBS=-lgmp
LIBOPTS=-L/opt/gmp/lib -Wl,--rpath -Wl,/opt/gmp/lib
INCS=-I/opt/gmp/include

all: gmseq ogmseq pg-ogmseq
	-@true

# unoptimised: for debugging
gmseq: gmseq.c
	gcc -g -O0 -Wall ${INCS} -o gmseq gmseq.c ${LIBOPTS} ${LIBS}

# optimised: for real use
ogmseq: gmseq.c
	gcc -g -O6 ${INCS} -o ogmseq gmseq.c ${LIBOPTS} ${LIBS}

# instrumented for profiling with gprof
pg-ogmseq: gmseq.c
	gcc -g -pg -O6 ${INCS} -o pg-ogmseq gmseq.c ${LIBOPTS} ${LIBS}
