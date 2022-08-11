#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "structs.h"
#include "readInputs.h"
#include "packing.h"

/*
Get the number of input surfactant molecules from user
Read and save all the coordinates and bond information from all the surfactant PDB molecule
Orient the surfactant molecule along X axis
Using the packing factor for each surfactant, calculate the angle between the surfactants
Using the angle between surfactant molecules in a single layer, find the total number of surfactant molecules in a layer
Using the number of surfactant molecules in a layer, calculate the number of layers to add
According to the packing angle, calculate the number of surfactant molecules to fit in one layer
Replicate the layer along the axis of the cylinder
Use bond information from the input files to create bond information for the assembled cylindrical micelle

IMPORTANT: Take all the input values from a config file
*/

int main(int argc, char const *argv[])
{
	checkArguments (argc);

	FILE *readConfig;
	readConfig = fopen (argv[1], "r");

	int nSurfactants = checkNSurfactants (readConfig);

	SURFACTANT *inputStructures;
	inputStructures = (SURFACTANT *) malloc (nSurfactants * sizeof (SURFACTANT));

	inputStructures = storeSurfactantInformation (inputStructures, nSurfactants, readConfig);
	inputStructures = getNAtoms (inputStructures, nSurfactants);
	inputStructures = getNBonds (inputStructures, nSurfactants);

	// Create variables to store atom coordinates and bond connectivity information
	COORDINATES **inputCoordinates;
	BONDS **inputBonds;

	inputCoordinates = readCoordinates (inputCoordinates, nSurfactants, inputStructures);
	inputBonds = readBonds (inputBonds, nSurfactants, inputStructures);

	FARTHESTPOINTS *inputStructures_farPoints;
	inputStructures_farPoints = (FARTHESTPOINTS *) malloc (nSurfactants * sizeof (FARTHESTPOINTS));

	inputStructures_farPoints = calculateFarPoints (inputStructures_farPoints, inputCoordinates, nSurfactants, inputStructures);

	inputCoordinates = orientSurfactants (inputCoordinates, nSurfactants, inputStructures, inputStructures_farPoints);

	// Find the longest dimension on X, Y, and Z
	CARTESIAN *loDimension, *hiDimension, globalSurfactantlo, globalSurfactanthi;
	loDimension = (CARTESIAN *) malloc (nSurfactants * sizeof (CARTESIAN));
	hiDimension = (CARTESIAN *) malloc (nSurfactants * sizeof (CARTESIAN));

	computeLongestDimension (&loDimension, &hiDimension, inputCoordinates, nSurfactants, inputStructures);
	calculateGlobalMinMax (&globalSurfactanthi, &globalSurfactantlo, loDimension, hiDimension, nSurfactants);

	// Replicate the molecule (coordinates and bonds)
	COORDINATES *outputCoordinates;
	BONDS *outputBonds;

	outputCoordinates = replicateSurfactants (inputCoordinates, inputBonds, globalSurfactantlo, globalSurfactanthi, nSurfactants, inputStructures);

	// Testing results from replicateSurfactants function
	// int totalAtoms = countTotalAtoms (inputStructures, nSurfactants);
	// for (int i = 0; i < totalAtoms; ++i)
	// {
	// 	printf("%f %f %f %s %s %s\n", outputCoordinates[i].x, outputCoordinates[i].y, outputCoordinates[i].z, outputCoordinates[i].atomName1, outputCoordinates[i].atomName2, outputCoordinates[i].molName);
	// 	fflush (stdout);
	// 	usleep (100000);
	// }

	outputBonds = addBonds (outputCoordinates, inputCoordinates, inputBonds, inputStructures, nSurfactants);

	// Testing results from addBonds function
	// int totalAtoms = countTotalAtoms (inputStructures, nSurfactants);
	// for (int i = 0; i < totalAtoms; ++i)
	// {
	// 	printf("%3d %3d %3d %3d %3d %3d\n", outputBonds[i].atom1, outputBonds[i].atom2, outputBonds[i].atom3, outputBonds[i].atom4, outputBonds[i].atom5, outputBonds[i].atom6);
	// 	fflush (stdout);
	// 	usleep (100000);
	// }

	// Save the above information as *.car and *.mdf files
	writeCar ();
	writeMdf ();

	free (inputStructures);
	fclose (readConfig);
	return 0;
}