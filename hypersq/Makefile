C_FILES = hypersq.c
H_FILES =

# We will be using gcc builtins that compile directly to single CPU
# instructions where available.
CC_ARCH = -march=native

GCC_MAJOR := $(shell gcc -dumpversion)
ifeq "${GCC_MAJOR}" "7"
  CC_EXTRA_OPT = -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im
endif
CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops ${CC_EXTRA_OPT}

DEFINES :=
dhypersq: CC_OPT = -O0
mhypersq: DEFINES += -DMULDIM

default: hypersq
all: dhypersq hypersq mhypersq

dhypersq hypersq mhypersq: Makefile ${C_FILES} ${H_FILES}
	gcc -o $@ -g ${CC_ARCH} ${DEFINES} ${CC_OPT} ${C_FILES}
