#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include "structs.h"
#include "readInputs.h"
#include "packing.h"
#include "printOutput.h"

int main(int argc, char const *argv[])
{
	checkArguments (argc);

	FILE *readConfig;
	readConfig = fopen (argv[1], "r");

	int nSurfactants = checkNSurfactants (readConfig);
	printf("Number of surfactant molecules: %d\n", nSurfactants);

	SURFACTANT *inputStructures;
	inputStructures = (SURFACTANT *) malloc (nSurfactants * sizeof (SURFACTANT));

	inputStructures = storeSurfactantInformation (inputStructures, nSurfactants, readConfig);
	inputStructures = getNAtoms (inputStructures, nSurfactants);
	inputStructures = getNBonds (inputStructures, nSurfactants);

	/* Printing the number of atoms for each surfactant */
	// for (int i = 0; i < nSurfactants; ++i)
	// {
	// 	printf("Number of atoms in surfactant %d: %d\n", i + 1, inputStructures[i].nAtoms);
	// }
	// exit (1);

	// Create variables to store atom coordinates and bond connectivity information
	COORDINATES **inputCoordinates;
	BONDS **inputBonds;

	inputCoordinates = readCoordinates (inputCoordinates, nSurfactants, inputStructures);
	inputBonds = readBonds (inputBonds, nSurfactants, inputStructures);

	FARTHESTPOINTS *inputStructures_farPoints;
	inputStructures_farPoints = (FARTHESTPOINTS *) malloc (nSurfactants * sizeof (FARTHESTPOINTS));

	inputStructures_farPoints = calculateFarPoints (inputStructures_farPoints, inputCoordinates, nSurfactants, inputStructures);

	printf("\nFar points for surfactant 2:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nsino: %d (%.2f, %.2f, %.2f)\nsino: %d (%.2f, %.2f, %.2f)\n\n", inputStructures_farPoints[1].sino1, inputStructures_farPoints[1].x1, inputStructures_farPoints[1].y1, inputStructures_farPoints[1].z1, inputStructures_farPoints[1].sino2, inputStructures_farPoints[1].x2, inputStructures_farPoints[1].y2, inputStructures_farPoints[1].z2);

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
	int totalAtoms = countTotalAtoms (inputStructures, nSurfactants);

	// Testing results from replicateSurfactants function
	// for (int i = 0; i < totalAtoms; ++i)
	// {
	// 	printf("%f %f %f %s %s %s\n", outputCoordinates[i].x, outputCoordinates[i].y, outputCoordinates[i].z, outputCoordinates[i].atomName1, outputCoordinates[i].atomName2, outputCoordinates[i].molName);
	// 	// fflush (stdout);
	// 	// usleep (10000);
	// }
	// exit (1);

	outputBonds = addBonds (outputCoordinates, inputCoordinates, inputBonds, inputStructures, nSurfactants);

	// Testing results from addBonds function
	// for (int i = 0; i < totalAtoms; ++i)
	// {
	// 	printf("%3d %3d %3d %3d %3d %3d\n", outputBonds[i].atom1, outputBonds[i].atom2, outputBonds[i].atom3, outputBonds[i].atom4, outputBonds[i].atom5, outputBonds[i].atom6);
	// 	fflush (stdout);
	// 	usleep (100000);
	// }

	// Save the above information as *.car and *.mdf files
	writeCar (outputCoordinates, totalAtoms, inputStructures, nSurfactants);
	writeMdf (outputCoordinates, outputBonds, totalAtoms, inputStructures, nSurfactants);

	printf("Creating lammps data file...\n");
	system ("./msi2lmp_gcc32.exe finalStructure");

	// printf("Adding polarizable gold substrate...\n");
	// CARTESIAN lowerBounds, upperBounds;
	// computeOutputBounds (&lowerBounds, &upperBounds, outputCoordinates, totalAtoms);
	// // Add a gold surface in XY plane.
	// // Surfactants are filled first in X axis, then in Y axis. Increment in Z axis was given at the end.
	// // So, dimensions along X and Y will be greater than in X direction.

	// printf("Filling the remaining area with water molecules...\n");

	free (inputStructures);
	fclose (readConfig);
	return 0;
}