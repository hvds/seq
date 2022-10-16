MPUGMP = /src/perl/Math-Prime-Util-GMP
COUL = coulfact.c diag.c rootmod.c coultau.c pell.c
HOUL = coulfact.h diag.h rootmod.h coultau.h pell.h coul.h

CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
dcoul dpcoul: CC_OPT = -O0

CFACTOR = ${MPUGMP}/factor.c ${MPUGMP}/prime_iterator.c ${MPUGMP}/ecm.c ${MPUGMP}/pbrent63.c ${MPUGMP}/isaac.c ${MPUGMP}/tinyqs.c ${MPUGMP}/squfof126.c ${MPUGMP}/simpqs.c ${MPUGMP}/primality.c ${MPUGMP}/utility.c ${MPUGMP}/gmp_main.c ${MPUGMP}/bls75.c ${MPUGMP}/real.c ${MPUGMP}/ecpp.c
HFACTOR = ${MPUGMP}/factor.h

DEFINES := -DSTANDALONE
pcoul dpcoul ftest: DEFINES += -DPARALLEL
ifdef SQONLY
    DEFINES += -DSQONLY
endif
ifdef CHECK_OVERFLOW
    DEFINES += -DCHECK_OVERFLOW
endif
ifdef VERBOSE
    DEFINES += -DVERBOSE
endif
ifdef LARGE_MIN
    DEFINES += -DLARGE_MIN
endif
ifdef TRY_HARDER
    DEFINES += -DTRY_HARDER
endif

all: pcoul

coul pcoul dcoul dpcoul: Makefile coul.c ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} coul.c ${COUL} ${CFACTOR} -I${MPUGMP} -lgmp -lm

test_pell: Makefile test_pell.c ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o test_pell -g ${CC_OPT} ${DEFINES} test_pell.c ${COUL} ${CFACTOR} -I${MPUGMP} -lgmp -lm

ftest: Makefile ftest.c coultau.c ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o ftest -g ${CC_OPT} ${DEFINES} ftest.c coultau.c ${CFACTOR} -I${MPUGMP} -lgmp -lm
