all: try gtry

try: try.c loc.c loc.h
	cat try.c loc.c >.whole.c
	gcc -o try -O6 -fwhole-program -g .whole.c

gtry: try.c loc.c loc.h
	gcc -o gtry -O0 -g try.c loc.c
