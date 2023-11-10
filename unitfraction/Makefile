MPUGMP = /src/perl/Math-Prime-Util-GMP

GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}
dfind_primitive dminset: CC_OPT := -O0

CFACTOR = ${MPUGMP}/factor.c ${MPUGMP}/ecm.c ${MPUGMP}/pbrent63.c ${MPUGMP}/isaac.c ${MPUGMP}/tinyqs.c ${MPUGMP}/squfof126.c ${MPUGMP}/simpqs.c ${MPUGMP}/primality.c ${MPUGMP}/utility.c ${MPUGMP}/gmp_main.c ${MPUGMP}/bls75.c ${MPUGMP}/real.c ${MPUGMP}/ecpp.c ${MPUGMP}/prime_iterator.c
HFACTOR = ${MPUGMP}/factor.h
DEFINES := -DSTANDALONE

dfind_primitive find_primitive: C_FILES := find_primitive.c
dminset minset: C_FILES := minset.c

default: find_primitive minset
all: dfind_primitive find_primitive dminset minset

dfind_primitive find_primitive dminset minset: Makefile ${C_FILES} unit.c unit.h ${CFACTOR} ${HFACTOR}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} ${C_FILES} unit.c ${CFACTOR} -I${MPUGMP} -lgmp -lm -lrt
