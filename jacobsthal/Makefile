GCC=LD_LIBRARY_PATH="/opt/mpfr-2.4.1/lib:/opt/mpc-0.8.1/lib" /opt/gcc-4.5.0/bin/gcc
CC_OPT=-O6

all: cj gcj 0cj
	true

cj: cj.c
	${GCC} ${CC_OPT} -o cj cj.c

gcj: cj.c
	${GCC} ${CC_OPT} -g -o gcj cj.c

0cj: cj.c
	${GCC} -O0  -g -o 0cj cj.c

clean:
	rm -f cj gcj 0cj
