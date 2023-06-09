CC=gcc

INCLUDES=
LIBS=
DEBUG=-O0
OPT=-O3
# Any improvement from these extra optimization options is lost in noise
#OPT=-O3 -fomit-frame-pointer -Wunsafe-loop-optimizations -funroll-loops -freorder-blocks-and-partition -funswitch-loops -fwhole-program

A169858: A169858.c
	${CC} -o A169858 -g ${OPT} A169858.c -I${INCLUDES} -L${LIBS} -l gmp

#Not yet written
#A177834: A177834.c
#	${CC} -o A177834 -g ${OPT} A177834.c -I${INCLUDES} -L${LIBS} -l gmp

#A175516(n) is a variant of A169858(n) allowing only even numbers.
A175516: A169858.c
	${CC} -o A175516 -g ${OPT} -DA175516 A169858.c -I${INCLUDES} -L${LIBS} -l gmp
