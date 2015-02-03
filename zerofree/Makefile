all: zerofree debug

zerofree: zerofree.c seen.h seen.c array.h array.c
	gcc -O6 -o zerofree zerofree.c seen.c array.c -lgmp

debug: zerofree.c seen.h seen.c array.h array.c
	gcc -g -o debug zerofree.c seen.c array.c -lgmp

clean:
	rm -f zerofree debug
