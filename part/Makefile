CFLAGS = -g -O6
COMMONCFILES = vec.c set.c clock.c symmetries.c
PARTCFILES = part.c pieces.c sym_set.c $(COMMONCFILES)
BUILDCFILES = build_pieces.c $(COMMONCFILES)
HFILES = part.h vec.h set.h clock.h symmetries.h pieces.h sym_set.h

PART = part1 part2 part3 part4 part5 part5r
BUILD = build1 build2 build3 build4 build5

part1: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=1 -o part1 $(CFLAGS) $(PARTCFILES)

part2: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=2 -o part2 $(CFLAGS) $(PARTCFILES)

part3: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=3 -o part3 $(CFLAGS) $(PARTCFILES)

part4: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=4 -DQUIET -o part4 $(CFLAGS) $(PARTCFILES)

part5r: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=5 -DQUIET -DREVERSE -o part5r $(CFLAGS) $(PARTCFILES)

part5: $(PARTCFILES) $(HFILES)
	gcc -DNBASE=5 -DQUIET -o part5 $(CFLAGS) $(PARTCFILES)

build1: $(BUILDCFILES) $(HFILES)
	gcc -DNBASE=1 -o build1 $(CFLAGS) $(BUILDCFILES)

build2: $(BUILDCFILES) $(HFILES)
	gcc -DNBASE=2 -o build2 $(CFLAGS) $(BUILDCFILES)

build3: $(BUILDCFILES) $(HFILES)
	gcc -DNBASE=3 -o build3 $(CFLAGS) $(BUILDCFILES)

build4: $(BUILDCFILES) $(HFILES)
	gcc -DNBASE=4 -o build4 $(CFLAGS) $(BUILDCFILES)

build5: $(BUILDCFILES) $(HFILES)
	gcc -DNBASE=5 -o build5 $(CFLAGS) $(BUILDCFILES)

allpart: $(PART)

allbuild: $(BUILD)

all: allpart allbuild
