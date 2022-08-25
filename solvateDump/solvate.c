#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "structs.h"
#include "inputs.h"

/*

This code takes dump file as input, adds water and creates a data file.

Number of arguments:
~~~~~~~~~~~~~~~~~~~~

argv[0] = ./program
argv[1] = input data filename
argv[2] = input dump filename
argv[3] = output data filename

*/

int main(int argc, char const *argv[])
{
	FILE *inputData, *inputDump, *output;

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafile, datafile_raw;

	datafile = readData (argv[1], &atoms, &bonds, &angles, &dihedrals, &impropers);

	int nAtoms = getNatoms (argv[2]);
	DUMP *traj;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	traj = readLastFrame (pipeString, nAtoms, dumpDimension);

	return 0;
}