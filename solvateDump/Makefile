all:
	gcc -c inputs.c -lm -Wall
	gcc -c solvate.c -lm -Wall
	gcc solvate.o inputs.o -o solvateDump -lm -Wall
run:
	./solvate finalStructure_c25d1.lammps05 lastframe.lammpstrj solvate.data
