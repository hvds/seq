GCC=gcc
INCS=-I.
LIBS=-lgmp
LIBOPTS=
CC_OPT=-O6 -fgcse-sm -fgcse-las -fgcse-after-reload -Wunsafe-loop-optimizations -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
CC_ALL_OPT=-fwhole-program
# ... -fipa-struct-reorg -fvect-cost-model -fprefetch-loop-arrays
# -fwhole-program (with rearranged makefile, presumably)
# -flto -ffunction-sections -fdata-sections (but remove -g, or get experimental)
# -fprofile-use ...?
# increase params with <--param name=value>:
#   max-delay-slot-insn-search, max-delay-slot-live-search, max-inline-insns-single=300, max-inline-insns-auto=50, large-function-insns=2700, large-function-growth=100, ?large-unit-insns=10000, inline-unit-growth=30, ipcp-unit-growth=10, early-inlining-insns=8, max-early-inliner-iterations, max-unswitch-insns, max-unswitch-level, max-reload-search-insns=100, max-cselib-memory-locations=500

all: debug opt
	true

.SUFFIXES: .c .oo .do

.c.oo:
	${GCC} -g ${CC_OPT} ${INCS} -c -o $*.oo $*.c
	# ${GCC} -g -pg ${CC_OPT} ${INCS} -c -o $*.oo $*.c

.c.do:
	${GCC} -g -DDEBUG ${INCS} -c -o $*.do $*.c
	# ${GCC} -g -DDEBUG -DDEBUG_GMP_LEAK ${INCS} -c -o $*.do $*.c

TEST=t/mygmp t/inverse t/mbh t/prime t/walker t/pp
DTEST=t/d_mygmp t/d_inverse t/d_mbh t/d_prime t/d_walker t/d_pp
CFILE=pp.c walker.c mygmp.c mbh.c prime.c inverse.c
HFILE=pp.h walker.h mygmp.h mbh.h prime.h inverse.h
PROG=d_C_int C_int

opt: C_int ${TEST}
debug: d_C_int ${DTEST}
prog: ${PROG}

all.c: C_int.c ${CFILE} ${HFILE}
	cat C_int.c ${CFILE} >all.c

C_int: all.c
	${GCC} -o $@ -g ${CC_OPT} ${CC_ALL_OPT} ${INCS} -DALL_C all.c ${LIBOPTS} ${LIBS}

d_C_int: C_int.c ${CFILE} ${HFILE}
	${GCC} -o $@ -g -DDEBUG ${INCS} C_int.c ${CFILE} ${LIBOPTS} ${LIBS}

${TEST} ${DTEST}:
	${GCC} -o $@ ${CC_OPT} $^ ${LIBOPTS} ${LIBS}
	# ${GCC} -pg -o $@ ${CC_OPT} $^ ${LIBOPTS} ${LIBS}

# Object dependencies
mygmp.do mygmp.oo: mygmp.c mygmp.h
inverse.do inverse.oo: inverse.c inverse.h
mbh.do mbh.oo: mbh.c mbh.h
prime.do prime.oo: prime.c prime.h mygmp.h
walker.do walker.oo: walker.c walker.h mbh.h mygmp.h
pp.do pp.oo: pp.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h
C_int.do C_int.oo: C_int.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h
t/mygmp.do t/mygmp.oo: t/mygmp.c mygmp.h
t/inverse.do t/inverse.oo: t/inverse.c inverse.h
t/mbh.do t/mbh.oo: t/mbh.c mbh.h
t/prime.do t/prime.oo: t/prime.c prime.h mygmp.h
t/walker.do t/walker.oo: t/walker.c walker.h mbh.h mygmp.h pp.h
t/pp.do t/pp.oo: t/pp.c pp.h walker.h mbh.h inverse.h prime.h mygmp.h

# Program dependencies
t/d_mygmp: t/mygmp.do mygmp.do
t/mygmp: t/mygmp.oo mygmp.oo
t/d_inverse: t/inverse.do inverse.do
t/inverse: t/inverse.oo inverse.oo
t/d_mbh: t/mbh.do mbh.do
t/mbh: t/mbh.oo mbh.oo
t/d_prime: t/prime.do prime.do mygmp.do
t/prime: t/prime.oo prime.oo mygmp.oo
t/d_walker: t/walker.do pp.do inverse.do mbh.do prime.do walker.do mygmp.do
t/walker: t/walker.oo pp.oo inverse.oo mbh.oo prime.oo walker.oo mygmp.oo
t/d_pp: t/pp.do pp.do inverse.do mbh.do prime.do walker.do mygmp.do
t/pp: t/pp.oo pp.oo inverse.oo mbh.oo prime.oo walker.oo mygmp.oo

dtest: debug
	echo -n "t/d_mygmp: " ; t/d_mygmp
	echo -n "t/d_inverse: " ; t/d_inverse
	echo -n "t/d_mbh: " ; t/d_mbh
	echo -n "t/d_prime: " ; t/d_prime
	echo -n "t/d_walker: " ; t/d_walker
	echo -n "t/d_pp: " ; t/d_pp

test: opt
	echo -n "t/mygmp: " ; t/mygmp
	echo -n "t/inverse: " ; t/inverse
	echo -n "t/mbh: " ; t/mbh
	echo -n "t/prime: " ; t/prime
	echo -n "t/walker: " ; t/walker
	echo -n "t/pp: " ; t/pp

clean:
	rm -f ${TEST} ${DTEST} C_int d_C_int all.c *.do *.oo t/*.do t/*.oo core* vgcore*
