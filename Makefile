all:
	gcc -c cylindricalMicelle.c -lm -Wall
	gcc -c packing.c -lm -Wall
	gcc -c readInputs.c -lm -Wall
	gcc -c printOutput.c -lm -Wall
	gcc cylindricalMicelle.o packing.o readInputs.o printOutput.o -o cylindricalMicelle -lm -Wall
