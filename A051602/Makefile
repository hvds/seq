CC_OPT=-O6 -fgcse-sm -fgcse-las -fgcse-after-reload -Wunsafe-loop-optimizations -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops
CC_ALL_OPT=-fwhole-program

all: recurse grecurse advance gadvance

recurse: recurse.c loc.c loc.h sym.c sym.c Makefile
	cat recurse.c loc.c sym.c >.recurse.c
	gcc -o recurse ${CC_OPT} ${CC_ALL_OPT} -g .recurse.c

advance: advance.c loc.c loc.h sym.c sym.c grid.h Makefile
	cat advance.c loc.c sym.c >.advance.c
	gcc -o advance ${CC_OPT} ${CC_ALL_OPT} -g .advance.c

grecurse: recurse.c loc.c loc.h sym.c sym.c Makefile
	gcc -DDEBUG -o grecurse -O0 -g recurse.c loc.c sym.c

gadvance: advance.c loc.c loc.h sym.c sym.c grid.h Makefile
	gcc -DDEBUG -o gadvance -O0 -g advance.c loc.c sym.c

madvance: advance.c loc.c loc.h sym.c sym.c grid.h Makefile
	clang -DDEBUG -o madvance -O0 -g -fsanitize=address advance.c loc.c sym.c
# ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./madvance ...
