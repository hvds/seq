CFILES = allsub.c
HFILES = cf.h allsub.h
DEFINES =

GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}

dcf dA051013: CC_OPT = -O0
dcf dA051013: DEFINES += -DDEBUG
dA051013 A051013: DEFINES += -DMAX_ALLSUB=64

default: dcf cf

all: dcf cf dA051013 A051013

cf dcf: Makefile cf.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_OPT} cf.c ${CFILES} ${DEFINES}

dA051013 A051013: Makefile A051013.c ${CFILES} ${HFILES}
	gcc -o $@ -g ${CC_OPT} A051013.c ${CFILES} ${DEFINES}
