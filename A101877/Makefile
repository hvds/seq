all: debug opt
	true

debug: t/d_inverse t/d_mbh t/d_prime t/d_usable t/d_walker t/d_pp
	true

opt: t/inverse t/mbh t/prime t/usable t/walker t/pp
	true

dC_int: C_int.c pp.c pp.h inverse.c inverse.h mbh.c mbh.h prime.c prime.h usable.c usable.h walker.c walker.h
	gcc -g -DDEBUG -I. -o dC_int C_int.c pp.c inverse.c mbh.c prime.c usable.c walker.c

C_int: C_int.c pp.c pp.h inverse.c inverse.h mbh.c mbh.h prime.c prime.h usable.c usable.h walker.c walker.h
	gcc -g -O6 -I. -o C_int C_int.c pp.c inverse.c mbh.c prime.c usable.c walker.c

dtest: debug
	echo -n "t/d_inverse: " ; t/d_inverse
	echo -n "t/d_mbh: " ; t/d_mbh
	echo -n "t/d_prime: " ; t/d_prime
	echo -n "t/d_usable: " ; t/d_usable
	echo -n "t/d_walker: " ; t/d_walker
	echo -n "t/d_pp: " ; t/d_pp

test: opt
	echo -n "t/inverse: " ; t/inverse
	echo -n "t/mbh: " ; t/mbh
	echo -n "t/prime: " ; t/prime
	echo -n "t/usable: " ; t/usable
	echo -n "t/walker: " ; t/walker
	echo -n "t/pp: " ; t/pp

t/d_inverse: t/inverse.c inverse.c inverse.h
	gcc -g -DDEBUG -I. -o t/d_inverse t/inverse.c inverse.c

t/inverse: t/inverse.c inverse.c inverse.h
	gcc -g -O6 -I. -o t/inverse t/inverse.c inverse.c

t/d_mbh: t/mbh.c mbh.c mbh.h
	gcc -g -DDEBUG -I. -o t/d_mbh t/mbh.c mbh.c

t/mbh: t/mbh.c mbh.c mbh.h
	gcc -g -O6 -I. -o t/mbh t/mbh.c mbh.c

t/d_prime: t/prime.c prime.c prime.h
	gcc -g -DDEBUG -I. -o t/d_prime t/prime.c prime.c

t/prime: t/prime.c prime.c prime.h
	gcc -g -O6 -I. -o t/prime t/prime.c prime.c

t/d_usable: t/usable.c usable.c usable.h prime.c prime.h
	gcc -g -DDEBUG -I. -o t/d_usable t/usable.c usable.c prime.c

t/usable: t/usable.c usable.c usable.h prime.c prime.h
	gcc -g -O6 -I. -o t/usable t/usable.c usable.c prime.c

t/d_walker: t/walker.c walker.c walker.h mbh.c mbh.h
	gcc -g -DDEBUG -I. -o t/d_walker t/walker.c walker.c mbh.c

t/walker: t/walker.c walker.c walker.h mbh.c mbh.h
	gcc -g -O6 -I. -o t/walker t/walker.c walker.c mbh.c

t/d_pp: t/pp.c pp.c pp.h inverse.c inverse.h mbh.c mbh.h prime.c prime.h usable.c usable.h walker.c walker.h
	gcc -g -DDEBUG -I. -o t/d_pp t/pp.c pp.c inverse.c mbh.c prime.c usable.c walker.c

t/pp: t/pp.c pp.c pp.h inverse.c inverse.h mbh.c mbh.h prime.c prime.h usable.c usable.h walker.c walker.h
	gcc -g -O6 -I. -o t/pp t/pp.c pp.c inverse.c mbh.c prime.c usable.c walker.c

clean:
	rm -f t/d_inverse t/inverse t/d_mbh t/mbh t/d_prime t/prime t/d_usable t/usable t/d_walker t/walker t/d_pp t/pp C_int dC_int core* vgcore*
