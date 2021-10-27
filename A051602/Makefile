CC_OPT=-O6 -fgcse-sm -fgcse-las -fgcse-after-reload -Wunsafe-loop-optimizations -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
CC_ALL_OPT=-fwhole-program

all: recurse grecurse

recurse: recurse.c loc.c loc.h sym.c sym.c Makefile
	cat recurse.c loc.c sym.c >.recurse.c
	gcc -o recurse ${CC_OPT} ${CC_ALL_OPT} -g .recurse.c

grecurse: recurse.c loc.c loc.h sym.c sym.c Makefile
	gcc -DDEBUG -o grecurse -O0 -g recurse.c loc.c sym.c
