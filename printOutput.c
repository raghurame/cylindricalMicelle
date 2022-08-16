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

	int atomID1, atomID2, atomID3, atomID4, atomID5;

	for (int i = 0; i < totalAtoms; ++i)
	{
		atomID1 = outputBonds[i].atom1;
		atomID2 = outputBonds[i].atom2;
		atomID3 = outputBonds[i].atom3;
		atomID4 = outputBonds[i].atom4;
		atomID5 = outputBonds[i].atom5;

		printf("%d %d %d %d\n", outputBonds[i].atom1, outputBonds[i].atom2, outputBonds[i].atom3, outputBonds[i].atom4, outputBonds[i].atom5);
		usleep (100000);

		fprintf(outputMDF, "XXXX_AF1:%s%-5d%5s%4s*%7s%6d%3d%11.4f%2d%2d%2d%7.4f%8.4f", outputCoordinates[i].atomName1, i, outputCoordinates[i].atomName1, outputCoordinates[i].atomName1, "?", 0, 0, outputCoordinates[i].col10, 0, 0, 8, 1.0, 0.0);

		// Printing information about bonded atoms
		if (atomID1 > 0) {
			fprintf(outputMDF, " %s%d ", outputCoordinates[atomID1 + 1].atomName1, atomID1 + 1); }
		if (atomID2 > 0) {
			fprintf(outputMDF, " %s%d ", outputCoordinates[atomID2 + 1].atomName1, atomID2 + 1); }
		if (atomID3 > 0) {
			fprintf(outputMDF, " %s%d ", outputCoordinates[atomID3 + 1].atomName1, atomID3 + 1); }
		if (atomID4 > 0) {
			fprintf(outputMDF, " %s%d ", outputCoordinates[atomID4 + 1].atomName1, atomID4 + 1); }
		if (atomID5 > 0) {
			fprintf(outputMDF, " %s%d ", outputCoordinates[atomID5 + 1].atomName1, atomID5 + 1); }

		fprintf(outputMDF, "\n");
	}

	fclose (outputMDF);
}