cA337663: main.c board.c board.h group.c group.h sym.c sym.h
	gcc -O6 -o cA337663 main.c board.c group.c sym.c
dA337663: main.c board.c board.h group.c group.h sym.c sym.h
	gcc -O0 -g -o dA337663 main.c board.c group.c sym.c
uA337663: main.c board.c board.h group.c group.h sym.c sym.h
	clang -g -o uA337663 -fsanitize=address main.c board.c group.c sym.c
# ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./uA337663
