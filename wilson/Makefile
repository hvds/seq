all: search-breadth search-depth

search-breadth: search-breadth.c breadth.c breadth.h mygmp.c mygmp.h
	gcc -g -O6 -Wall -o search-breadth search-breadth.c breadth.c mygmp.c -lgmp

search-depth: search-depth.c depth.c depth.h mygmp.c mygmp.h
	gcc -g -O6 -Wall -o search-depth search-depth.c depth.c mygmp.c -lgmp
