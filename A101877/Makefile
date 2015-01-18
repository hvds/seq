LIBS=-lgmp
LIBOPTS=
INCS=

all: C_int dC_int mC_int
	-@true

clean:
	rm ${obj} ${obj_d} ${obj_m} C_int dC_int mC_int

.SUFFIXES: .c .o .do .mo

.c.o:
	gcc -g -O6 ${INCS} -c -o $*.o $*.c

.c.do:
	gcc -g -O0 -DDEBUG -Wall ${INCS} -c -o $*.do $*.c

.c.mo:
	gcc -g -O0 -DDEBUG -DDEBUG_GMP_LEAK -Wall ${INCS} -c -o $*.mo $*.c

obj = C_int.o vector.o match.o total.o
obj_d = C_int.do vector.do match.do total.do
obj_m = C_int.mo vector.mo match.mo total.mo

# unoptimised: for debugging
dC_int: $(obj_d)
	gcc -g -O0 -DDEBUG -Wall -o dC_int $(obj_d) ${LIBOPTS} ${LIBS}

# memory debugging (mpq leaks)
mC_int: $(obj_m)
	gcc -g -O0 -DDEBUG -DDEBUG_GMP_LEAK -Wall -o mC_int $(obj_m) ${LIBOPTS} ${LIBS}

# optimised: for real use
C_int: $(obj)
	gcc -g -O6 -o C_int $(obj) ${LIBOPTS} ${LIBS}

# instrumented for profiling with gprof
#pC_int: C_int.c
#	gcc -g -pg -O6 ${INCS} -o pC_int C_int.c ${LIBOPTS} ${LIBS}
