MPUGMP = /src/perl/Math-Prime-Util-GMP
PERL = /opt/maths-GMP/lib/perl5/5.34.0/x86_64-linux/CORE
COUL = coulfact.c diag.c rootmod.c coultau.c pell.c
HOUL = coulfact.h diag.h rootmod.h coultau.h pell.h
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
#CC_OPT = -Og
# CC_ALL_OPT = -fwhole-program

CFACTOR = ${MPUGMP}/factor.c ${MPUGMP}/prime_iterator.c ${MPUGMP}/ecm.c ${MPUGMP}/pbrent63.c ${MPUGMP}/isaac.c ${MPUGMP}/tinyqs.c ${MPUGMP}/squfof126.c ${MPUGMP}/simpqs.c ${MPUGMP}/primality.c ${MPUGMP}/utility.c ${MPUGMP}/gmp_main.c ${MPUGMP}/bls75.c ${MPUGMP}/real.c ${MPUGMP}/ecpp.c
HFACTOR = ${MPUGMP}/factor.h

all: coul pcoul

coul: Makefile coul.c ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o coul -g ${CC_OPT} -DSTANDALONE coul.c ${COUL} ${CFACTOR} -I${MPUGMP} -I${PERL} -lgmp -lm

pcoul: Makefile coul.c ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o pcoul -DPARALLEL -g ${CC_OPT} -DSTANDALONE coul.c ${COUL} ${CFACTOR} -I${MPUGMP} -I${PERL} -lgmp -lm

test_pell: Makefile test_pell.c ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o test_pell -g ${CC_OPT} -DSTANDALONE test_pell.c ${COUL} ${CFACTOR} -I${MPUGMP} -I${PERL} -lgmp -lm
