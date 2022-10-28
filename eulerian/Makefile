CC_OPT = -O6 -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftracer -funroll-loops -fvariable-expansion-in-unroller -freorder-blocks-and-partition -funswitch-loops -ftree-loop-linear -ftree-loop-distribution -ftree-loop-im

eulerian: Makefile eulerian.c diag.c diag.h
	gcc -g -o eulerian ${CC_OPT} eulerian.c diag.c

debug: Makefile eulerian.c diag.c diag.h
	gcc -g -o debug -O0 eulerian.c diag.c
