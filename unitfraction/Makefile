GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}

C_FILES := unit.c diag.c sieve.c
H_FILES := unit.h diag.h sieve.h

dfind_primitive dminset: CC_OPT := -O0

default: find_primitive minset
all: dfind_primitive find_primitive dminset minset

dfind_primitive find_primitive: Makefile find_primitive.c ${C_FILES} ${H_FILES}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} find_primitive.c ${C_FILES} -lgmp -lm -lrt

dminset minset: Makefile minset.c ${C_FILES} ${H_FILES}
	gcc -o $@ -g ${CC_OPT} ${DEFINES} minset.c ${C_FILES} -lgmp -lm -lrt
