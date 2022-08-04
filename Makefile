all:
	gcc -c cylindricalMicelle.c -lm -Wall
	gcc -c packing.c -lm -Wall
	gcc -c readInputs.c -lm -Wall
	gcc cylindricalMicelle.o packing.o readInputs.o -o cylindricalMicelle -lm -Wall
