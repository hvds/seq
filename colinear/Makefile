all: dbreadth obreadth

obreadth: Makefile breadth.c poly.c poly.h
	gcc -g -O6 -o obreadth breadth.c poly.c

dbreadth: Makefile breadth.c poly.c poly.h
	gcc -g -O0 -o dbreadth breadth.c poly.c
