GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}
dfind_primitive dminset: CC_OPT := -O0
dfind_primitive find_primitive: C_FILES := find_primitive.c
dminset minset: C_FILES := minset.c

default: find_primitive minset
all: dfind_primitive find_primitive dminset minset

dfind_primitive find_primitive dminset minset: Makefile ${C_FILES} unit.c unit.h
	gcc -o $@ -g ${CC_OPT} ${C_FILES} unit.c -lgmp -lm -lrt
