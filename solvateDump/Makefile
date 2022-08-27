all:
	gcc -o solvate solvate.c -lm -Wall
	./solvate finalStructure_c25d1.lammps05 lastframe.lammpstrj solvate.data
compile:
	gcc -o solvate solvate.c -lm -Wall
run:
	./solvate finalStructure_c25d1.lammps05 lastframe.lammpstrj solvate.data
