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

void writeCar (COORDINATES *outputCoordinates, int totalAtoms, SURFACTANT *inputStructures, int nSurfactants)
{
	FILE *outputCAR;
	outputCAR = fopen ("finalStructure.car", "w");

	for (int i = 0; i < totalAtoms; ++i)
	{
		printf("%s%-5d%14.9f%14.9f%14.9f XXXX AF1%6s*%7s%8.3f\n", outputCoordinates[i].atomName1, i, outputCoordinates[i].x, outputCoordinates[i].y, outputCoordinates[i].z, outputCoordinates[i].atomName1, outputCoordinates[i].atomName1, outputCoordinates[i].col10);
		// printf("%-5d%6f\n", 23, 56);
		fflush (stdout);
		usleep (100000);
	}

	fclose (outputCAR);
}

void writeMdf (BONDS *outputBonds, int totalAtoms, SURFACTANT *inputStructures, int nSurfactants)
{
	FILE *outputMDF;
	outputMDF = fopen ("finalStructure.mdf", "w");

	fclose (outputMDF);
}