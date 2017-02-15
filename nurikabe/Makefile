all: cframe cframe-d

cframe: cframe.c Makefile
	gcc -o cframe -g -O6 cframe.c

cframe-d: cframe.c Makefile
	gcc -o cframe-d -g -O0 -DDEBUG cframe.c
