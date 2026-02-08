all: dcolinear scolinear

scolinear: Makefile scolinear.c
	gcc -g -O6 -o scolinear scolinear.c

dcolinear: Makefile scolinear.c
	gcc -g -O0 -o dcolinear scolinear.c
