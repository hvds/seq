all: try gtry

try: try.c grid.c grid.h
	cat try.c grid.c >.whole.c
	gcc -o try -O6 -fwhole-program -g .whole.c

gtry: try.c grid.c grid.h
	gcc -o gtry -O0 -g try.c grid.c
