# main program, aggressively optimized for sped
cA337663: whole_file.c main.c board.c board.h group.c group.h sym.c sym.h loc.h
	gcc -O3 -o cA337663 -fwhole-program whole_file.c

# debug version, no optimization
dA337663: main.c board.c board.h group.c group.h sym.c sym.h loc.h
	gcc -O0 -g -o dA337663 main.c board.c group.c sym.c

# debug version, for finding bounds errors and memory leaks
uA337663: main.c board.c board.h group.c group.h sym.c sym.h loc.h
	clang -g -o uA337663 -fsanitize=address main.c board.c group.c sym.c
# ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./uA337663

# tests
t/test_sym: t/test_sym.c t/test_sym.h t/common.c t/common.h sym.c sym.h loc.h
	gcc -O3 -g -o t/test_sym t/test_sym.c t/common.c sym.c
t/test_group: t/test_group.c t/test_group.h t/common.c t/common.h group.c group.h sym.c sym.h loc.h
	gcc -O3 -g -o t/test_group t/test_group.c t/common.c group.c sym.c
