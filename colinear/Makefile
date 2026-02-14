all: dbreadth ddepth obreadth odepth

obreadth: Makefile breadth.c poly.c poly.h
	gcc -g -O6 -o obreadth breadth.c poly.c

dbreadth: Makefile breadth.c poly.c poly.h
	gcc -g -O0 -o dbreadth breadth.c poly.c

odepth: Makefile depth.c poly.c poly.h
	gcc -g -O6 -o odepth depth.c poly.c ../inc/diag.c -lgmp -lrt

ddepth: Makefile depth.c poly.c poly.h
	gcc -g -O0 -o ddepth depth.c poly.c ../inc/diag.c -lgmp -lrt
