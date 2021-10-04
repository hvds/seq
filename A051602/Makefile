CC_OPT=-O6 -fgcse-sm -fgcse-las -fgcse-after-reload -Wunsafe-loop-optimizations -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
CC_ALL_OPT=-fwhole-program

all: try gtry

try: try.c loc.c loc.h sym.c sym.c Makefile
	cat try.c loc.c sym.c >.whole.c
	gcc -o try ${CC_OPT} ${CC_ALL_OPT} -g .whole.c

gtry: try.c loc.c loc.h sym.c sym.c Makefile
	gcc -DDEBUG -o gtry -O0 -g try.c loc.c sym.c
