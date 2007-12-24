gmseq: gmseq.o
	gcc -o gmseq gmseq.o -L/opt/gmp/lib -Wl,--rpath -Wl,/opt/gmp/lib -lgmp

gmseq.o: gmseq.c
	gcc -c -g -O0 -I/opt/gmp/include -o gmseq.o gmseq.c

ogmseq: ogmseq.o
	gcc -o ogmseq ogmseq.o -L/opt/gmp/lib -Wl,--rpath -Wl,/opt/gmp/lib -lgmp

ogmseq.o: gmseq.c
	gcc -c -g -O6 -I/opt/gmp/include -o ogmseq.o gmseq.c
