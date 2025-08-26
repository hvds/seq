CFILES = cf.c
HFILES =

GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}
dcf: CC_OPT = -O0

default: dcf cf

all: dcf cf

cf dcf: Makefile ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_OPT} ${CFILES}
