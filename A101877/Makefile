INCS=-I.
LIBS=-lgmp
LIBOPTS=

all: debug opt
	true

.SUFFIXES: .c .oo .do

.c.oo:
	gcc -g -O6 ${INCS} -c -o $*.oo $*.c

.c.do:
	gcc -g -DDEBUG ${INCS} -c -o $*.do $*.c
	# gcc -g -DDEBUG -DDEBUG_GMP_LEAK ${INCS} -c -o $*.do $*.c

TEST=t/inverse t/mbh t/prime t/walker t/pp
DTEST=t/d_inverse t/d_mbh t/d_prime t/d_walker t/d_pp

opt: C_int ${TEST}
debug: d_C_int ${DTEST}

C_int d_C_int ${TEST} ${DTEST}:
	gcc -o $@ $^ ${LIBOPTS} ${LIBS}

# Object dependencies
mygmp.do mygmp.oo: mygmp.c mygmp.h
inverse.do inverse.oo: inverse.c inverse.h
mbh.do mbh.oo: mbh.c mbh.h
prime.do prime.oo: prime.c prime.h mygmp.h
walker.do walker.oo: walker.c walker.h mbh.h mygmp.h
pp.do pp.oo: pp.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h
C_int.do C_int.oo: C_int.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h
t/inverse.do t/inverse.oo: t/inverse.c inverse.h
t/mbh.do t/mbh.oo: t/mbh.c mbh.h
t/prime.do t/prime.oo: t/prime.c prime.h mygmp.h
t/walker.do t/walker.oo: t/walker.c walker.h mbh.h mygmp.h
t/pp.do t/pp.oo: t/pp.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h

# Program dependencies
t/d_inverse: t/inverse.do inverse.do
t/inverse: t/inverse.oo inverse.oo
t/d_mbh: t/mbh.do mbh.do
t/mbh: t/mbh.oo mbh.oo
t/d_prime: t/prime.do prime.do mygmp.do
t/prime: t/prime.oo prime.oo mygmp.oo
t/d_walker: t/walker.do walker.do mbh.do mygmp.do
t/walker: t/walker.oo walker.oo mbh.oo mygmp.oo
t/d_pp: t/pp.do pp.do inverse.do mbh.do prime.do walker.do mygmp.do
t/pp: t/pp.oo pp.oo inverse.oo mbh.oo prime.oo walker.oo mygmp.oo
d_C_int: C_int.do pp.do inverse.do mbh.do prime.do walker.do mygmp.do
C_int: C_int.oo pp.oo inverse.oo mbh.oo prime.oo walker.oo mygmp.oo

dtest: debug
	echo -n "t/d_inverse: " ; t/d_inverse
	echo -n "t/d_mbh: " ; t/d_mbh
	echo -n "t/d_prime: " ; t/d_prime
	echo -n "t/d_walker: " ; t/d_walker
	echo -n "t/d_pp: " ; t/d_pp

test: opt
	echo -n "t/inverse: " ; t/inverse
	echo -n "t/mbh: " ; t/mbh
	echo -n "t/prime: " ; t/prime
	echo -n "t/walker: " ; t/walker
	echo -n "t/pp: " ; t/pp

clean:
	rm -f ${TEST} ${DTEST} C_int d_C_int *.do *.oo t/*.do t/*.oo core* vgcore*
