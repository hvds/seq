CFLAGS = -g -O6
CFILES = part.c vec.c set.c clock.c symmetries.c pieces.c sym_set.c
HFILES = part.h vec.h set.h clock.h symmetries.h pieces.h sym_set.h

part1: $(CFILES) $(HFILES)
	gcc -DNBASE=1 -o part1 $(CFLAGS) $(CFILES)

part2: $(CFILES) $(HFILES)
	gcc -DNBASE=2 -o part2 $(CFLAGS) $(CFILES)

part3: $(CFILES) $(HFILES)
	gcc -DNBASE=3 -o part3 $(CFLAGS) $(CFILES)

part4: $(CFILES) $(HFILES)
	gcc -DNBASE=4 -o part4 $(CFLAGS) $(CFILES)

#part5p: $(CFILES) $(HFILES)
#	gcc -DNBASE=5 -DPIECE_ONLY -o part5 $(CFLAGS) $(CFILES)

part5r: $(CFILES) $(HFILES)
	gcc -DNBASE=5 -DREVERSE -o part5r $(CFLAGS) $(CFILES)

part5: $(CFILES) $(HFILES)
	gcc -DNBASE=5 -DNOPREP -o part5 $(CFLAGS) $(CFILES)

all: part1 part2 part3 part4 part5 part5r
