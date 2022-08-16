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

	for (int i = 0; i < totalAtoms; ++i) {
		fprintf(outputCAR, "%s%-5d%14.9f%14.9f%14.9f XXXX AF1%6s*%7s%8.3f\n", outputCoordinates[i].atomName1, i, outputCoordinates[i].x, outputCoordinates[i].y, outputCoordinates[i].z, outputCoordinates[i].atomName1, outputCoordinates[i].atomName1, outputCoordinates[i].col10); }

	fclose (outputCAR);
}

void writeMdf (COORDINATES *outputCoordinates, BONDS *outputBonds, int totalAtoms, SURFACTANT *inputStructures, int nSurfactants)
{
	FILE *outputMDF;
	outputMDF = fopen ("finalStructure.mdf", "w");

	for (int i = 0; i < totalAtoms; ++i)
	{
		fprintf(outputMDF, "XXXX_AF1:%s%-5d%5s%4s*%7s%6d%3d%11.4f%2d%2d%2d%7.4f%8.4f\n", outputCoordinates[i].atomName1, i, outputCoordinates[i].atomName1, outputCoordinates[i].atomName1, "?", 0, 0, outputCoordinates[i].col10, 0, 0, 8, 1.0, 0.0);
	}

	fclose (outputMDF);
}