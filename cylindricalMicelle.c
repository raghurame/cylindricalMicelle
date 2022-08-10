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

float calculateDistance (float x1, float y1, float z1, float x2, float y2, float z2)
{
	float distance;

	distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

	return distance;
}

FARTHESTPOINTS *calculateFarPoints (FARTHESTPOINTS *inputStructures_farPoints, COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{
	float distance_local, distance_max = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		distance_max = 0;
		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			for (int k = 0; k < inputStructures[i].nAtoms; ++k)
			{
				distance_local = calculateDistance (inputCoordinates[i][j].x, inputCoordinates[i][j].y, inputCoordinates[i][j].z, inputCoordinates[i][k].x, inputCoordinates[i][k].y, inputCoordinates[i][k].z);

				if (distance_local > distance_max)
				{
					distance_max = distance_local;

					// Here, [i] corresponds to the surfactant molecule
					inputStructures_farPoints[i].sino1 = j + 1;
					inputStructures_farPoints[i].sino2 = k + 1;
					inputStructures_farPoints[i].x1 = inputCoordinates[i][j].x;
					inputStructures_farPoints[i].y1 = inputCoordinates[i][j].y;
					inputStructures_farPoints[i].z1 = inputCoordinates[i][j].z;
					inputStructures_farPoints[i].x2 = inputCoordinates[i][k].x;
					inputStructures_farPoints[i].y2 = inputCoordinates[i][k].y;
					inputStructures_farPoints[i].z2 = inputCoordinates[i][k].z;
					inputStructures_farPoints[i].distance = distance_max;
				}
			}
		}
	}

	return inputStructures_farPoints;
}

CENTEROFMASS *computeCOM (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{
	CENTEROFMASS *comForMolecules;
	comForMolecules = (CENTEROFMASS *) malloc (nSurfactants * sizeof (CENTEROFMASS));

	for (int i = 0; i < nSurfactants; ++i)
	{
		// Resetting the x, y, z values, which will be incremented later
		comForMolecules[i].x = 0;
		comForMolecules[i].y = 0;
		comForMolecules[i].z = 0;

		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			comForMolecules[i].x += inputCoordinates[i][j].x;
			comForMolecules[i].y += inputCoordinates[i][j].y;
			comForMolecules[i].z += inputCoordinates[i][j].z;
		}

		comForMolecules[i].x /= inputStructures[i].nAtoms;
		comForMolecules[i].y /= inputStructures[i].nAtoms;
		comForMolecules[i].z /= inputStructures[i].nAtoms;
	}

	return comForMolecules;
}

COORDINATES **translateCoordinates (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, CENTEROFMASS *comForMolecules, FARTHESTPOINTS *inputStructures_farPoints)
{
	CARTESIAN location_farPoint1;
	int sino1_local;

	for (int i = 0; i < nSurfactants; ++i)
	{
		// Calculate the distance between far point-1 and the center of mass
		sino1_local = inputStructures_farPoints[i].sino1;
		location_farPoint1.x = inputCoordinates[i][sino1_local].x;
		location_farPoint1.y = inputCoordinates[i][sino1_local].y;
		location_farPoint1.z = inputCoordinates[i][sino1_local].z;
		// distanceFromCOM.x = inputCoordinates[i][sino1_local].x - comForMolecules[i].x;
		// distanceFromCOM.y = inputCoordinates[i][sino1_local].y - comForMolecules[i].y;
		// distanceFromCOM.z = inputCoordinates[i][sino1_local].z - comForMolecules[i].z;

		// Translate all atoms based on the above calculated distance

		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			inputCoordinates[i][j].x -= location_farPoint1.x;
			inputCoordinates[i][j].y -= location_farPoint1.y;
			inputCoordinates[i][j].z -= location_farPoint1.z;
		}
	}

	return inputCoordinates;
}

COORDINATES **alignMolecule (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, FARTHESTPOINTS *inputStructures_farPoints)
{
	// vector 1 is from sino1 and sino2 of the farthest points.
	// vector 2 is vectorX, which is a unit vector along X axis.
	// Find the angle between the two vectors.
	// Then rotate the molecule by the same angle, but in negative side.
	float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, dotProduct, magnitude1, magnitude2, cosTheta, theta;

	for (int i = 0; i < nSurfactants; ++i)
	{
		x1 = inputStructures_farPoints[i].x1;
		y1 = inputStructures_farPoints[i].y1;
		z1 = inputStructures_farPoints[i].z1;

		x2 = inputStructures_farPoints[i].x2;
		y2 = inputStructures_farPoints[i].y2;
		z2 = inputStructures_farPoints[i].z2;

		x3 = 0;
		y3 = 0;
		z3 = 0;

		x4 = 1;
		y4 = 0;
		z4 = 0;

		// Finding angle in XY plane
		dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3));
		magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)); 
		magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)); 

		cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2));
		theta = acosf (cosTheta);

		// Aligning the molecule along XY plane
		// Not modifying the Z values
		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			// printf("mol %d: %.2f, %.2f --> ", i, inputCoordinates[i][j].x, inputCoordinates[i][j].y);
			inputCoordinates[i][j].x = (inputCoordinates[i][j].x * cosf (theta)) - (inputCoordinates[i][j].y * sinf (theta));
			inputCoordinates[i][j].y = (inputCoordinates[i][j].x * sinf (theta)) - (inputCoordinates[i][j].y * cosf (theta));
			// printf("%.2f, %.2f\n", inputCoordinates[i][j].x, inputCoordinates[i][j].y);
		}

		// Finding angle in XZ plane
		dotProduct = ((x2 - x1) * (x4 - x3)) + ((z2 - z1) * (z4 - z3));
		magnitude1 = ((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)); 
		magnitude2 = ((x4 - x3) * (x4 - x3)) + ((z4 - z3) * (z4 - z3)); 

		cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2));
		theta = acosf (cosTheta);

		// Aligning the molecule along XZ plane
		// Here Y values are not modified
		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			inputCoordinates[i][j].x = (inputCoordinates[i][j].x * cosf (theta)) - (inputCoordinates[i][j].z * sinf (theta));
			inputCoordinates[i][j].z = (inputCoordinates[i][j].x * sinf (theta)) - (inputCoordinates[i][j].z * cosf (theta));
		}
	}

	return inputCoordinates;
}

COORDINATES **orientSurfactants (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, FARTHESTPOINTS *inputStructures_farPoints)
{
	CENTEROFMASS *comForMolecules;
	comForMolecules = computeCOM (inputCoordinates, nSurfactants, inputStructures); // assuming all atoms have same mass

	// Printing the center of mass for each molecule
	for (int i = 0; i < nSurfactants; ++i)
	{
		printf("center of mass for mol: %d => %f %f %f\n", i + 1, comForMolecules[i].x, comForMolecules[i].y, comForMolecules[i].z);
	}

	inputCoordinates = translateCoordinates (inputCoordinates, nSurfactants, inputStructures, comForMolecules, inputStructures_farPoints);

	// Aligning the molecule
	inputCoordinates = alignMolecule (inputCoordinates, nSurfactants, inputStructures, inputStructures_farPoints);

	return inputCoordinates;
}

void computeLongestDimension (CARTESIAN **loDimension, CARTESIAN **hiDimension, COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{	
	for (int i = 0; i < nSurfactants; ++i)
	{
		for (int j = 0; j < inputStructures[i].nAtoms; ++j)
		{
			if (j == 0)
			{
				(*loDimension)[i].x = inputCoordinates[i][j].x;
				(*hiDimension)[i].x = inputCoordinates[i][j].x;
				(*loDimension)[i].y = inputCoordinates[i][j].y;
				(*hiDimension)[i].y = inputCoordinates[i][j].y;
				(*loDimension)[i].z = inputCoordinates[i][j].z;
				(*hiDimension)[i].z = inputCoordinates[i][j].z;
			}

			if (inputCoordinates[i][j].x < (*loDimension)[i].x) (*loDimension)[i].x = inputCoordinates[i][j].x;
			if (inputCoordinates[i][j].x > (*hiDimension)[i].x) (*hiDimension)[i].x = inputCoordinates[i][j].x;
			if (inputCoordinates[i][j].y < (*loDimension)[i].y) (*loDimension)[i].y = inputCoordinates[i][j].y;
			if (inputCoordinates[i][j].y > (*hiDimension)[i].y) (*hiDimension)[i].y = inputCoordinates[i][j].y;
			if (inputCoordinates[i][j].z < (*loDimension)[i].z) (*loDimension)[i].z = inputCoordinates[i][j].z;
			if (inputCoordinates[i][j].z > (*hiDimension)[i].z) (*hiDimension)[i].z = inputCoordinates[i][j].z;
		}
	}	
}

int countTotalAtoms (SURFACTANT *inputStructures, int nSurfactants)
{
	int totalAtoms = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		totalAtoms += (inputStructures[i].nAtoms * inputStructures[i].nMolecules);
	}

	return totalAtoms;
}

int countTotalMolecules (SURFACTANT *inputStructures, int nSurfactants)
{
	int totalMolecules = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		totalMolecules += inputStructures[i].nMolecules;
	}

	return totalMolecules;
}

MOLECULELOG *initializeAllMoleculesLog (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants)
{
	int currentMolecule = 0;
	for (int i = 0; i < nSurfactants; ++i)
	{
		for (int j = 0; j < inputStructures[i].nMolecules; ++j)
		{
			allMoleculeLog[currentMolecule].surfactantID = i;
			allMoleculeLog[currentMolecule].fillStatus = 1; // Initializing with '1's
			currentMolecule++;
		}
	}

	return allMoleculeLog;
}

int findRemainingMolecules (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants)
{
	int remainingMolecules = 0, currentMolecule = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		for (int j = 0; j < inputStructures[i].nMolecules; ++j)
		{
			remainingMolecules += allMoleculeLog[currentMolecule].fillStatus;
			currentMolecule++;
		}
	}

	return remainingMolecules;
}

MOLECULELOG *pickMolecule (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants, double randomFlip, int *pickedSurfactantMolID)
{
	int currentMolecule = 0, currentUnselectedMolecule = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		for (int j = 0; j < inputStructures[i].nMolecules; ++j)
		{
			if (allMoleculeLog[currentMolecule].fillStatus == 1)
			{
				if (currentUnselectedMolecule == (int) randomFlip)
				{
					allMoleculeLog[currentMolecule].fillStatus = 0;
					*pickedSurfactantMolID = allMoleculeLog[currentMolecule].surfactantID;
				}
				currentUnselectedMolecule++;
			}
			currentMolecule++;
		}
	}

	return allMoleculeLog;
}

CARTESIAN computeMaxSurfactantLength (CARTESIAN maxSurfactantLenght, CARTESIAN globalSurfactantlo, CARTESIAN globalSurfactanthi, float tolerance)
{
	maxSurfactantLenght.x = (globalSurfactanthi.x - globalSurfactantlo.x) + tolerance;
	maxSurfactantLenght.y = (globalSurfactanthi.y - globalSurfactantlo.y) + tolerance;
	maxSurfactantLenght.z = (globalSurfactanthi.z - globalSurfactantlo.z) + tolerance;

	return maxSurfactantLenght;
}

COORDINATES *replicateSurfactants (COORDINATES **inputCoordinates, BONDS **inputBonds, CARTESIAN globalSurfactantlo, CARTESIAN globalSurfactanthi, int nSurfactants, SURFACTANT *inputStructures)
{
	COORDINATES *outputCoordinates;

	// Finding the max surfactant length, to calculate the lattice positions for translations
	CARTESIAN maxSurfactantLenght; float tolerance = 2.0;
	maxSurfactantLenght = computeMaxSurfactantLength (maxSurfactantLenght, globalSurfactantlo, globalSurfactanthi, tolerance);

	// Random picking of surfactant molecules
	srand (time (NULL));
	double randomFlip;
	int totalAtoms = countTotalAtoms (inputStructures, nSurfactants), totalMolecules = countTotalMolecules (inputStructures, nSurfactants), remainingMolecules = totalMolecules, pickedSurfactantMolID;

	outputCoordinates = (COORDINATES *) malloc (totalAtoms * sizeof (COORDINATES));

	MOLECULELOG *allMoleculeLog;
	allMoleculeLog = (MOLECULELOG *) malloc (totalMolecules * sizeof (MOLECULELOG));

	// Variables for creating surfactant lattice for translations
	int nMoleculesPerSide = ceil (cbrt (totalMolecules)), currentLatticeInX = 0, currentLatticeInY = 0, currentLatticeInZ = 0, xIncrement = 0, yIncrement = 0, zIncrement = 0, currentAtom = 0;

	allMoleculeLog = initializeAllMoleculesLog (allMoleculeLog, inputStructures, nSurfactants);

	for (int i = 0; i < totalMolecules; ++i)
	{
		randomFlip = rand ()/(double) RAND_MAX;
		remainingMolecules = findRemainingMolecules (allMoleculeLog, inputStructures, nSurfactants);
		randomFlip *= (double) remainingMolecules;

		allMoleculeLog = pickMolecule (allMoleculeLog, inputStructures, nSurfactants, randomFlip, &pickedSurfactantMolID);
		// printf("%d\n", pickedSurfactantMolID);
		// usleep (10000);

		for (int j = 0; j < inputStructures[pickedSurfactantMolID].nAtoms; ++j)
		{
			outputCoordinates[currentAtom].x = inputCoordinates[pickedSurfactantMolID][j].x + (maxSurfactantLenght.x * xIncrement);
			outputCoordinates[currentAtom].y = inputCoordinates[pickedSurfactantMolID][j].y + (maxSurfactantLenght.y * yIncrement);
			outputCoordinates[currentAtom].z = inputCoordinates[pickedSurfactantMolID][j].z + (maxSurfactantLenght.z * zIncrement);
			outputCoordinates[currentAtom].sino = inputCoordinates[pickedSurfactantMolID][j].sino;
			outputCoordinates[currentAtom].col5 = inputCoordinates[pickedSurfactantMolID][j].col5;
			outputCoordinates[currentAtom].col9 = inputCoordinates[pickedSurfactantMolID][j].col9;
			outputCoordinates[currentAtom].col10 = inputCoordinates[pickedSurfactantMolID][j].col10;

			strncpy (outputCoordinates[currentAtom].atomName1, inputCoordinates[pickedSurfactantMolID][j].atomName1, 5);
			strncpy (outputCoordinates[currentAtom].atomName2, inputCoordinates[pickedSurfactantMolID][j].atomName2, 5);
			strncpy (outputCoordinates[currentAtom].molName, inputCoordinates[pickedSurfactantMolID][j].molName, 5);

			if (currentAtom < totalAtoms) {
				currentAtom++; }
		}

		// printf("%d %d %d %d\n", nMoleculesPerSide, xIncrement, yIncrement, zIncrement);
		// fflush (stdout);
		// usleep (100000);

		if (xIncrement == (nMoleculesPerSide - 1)) {
			xIncrement = 0;
			yIncrement++; }
		else if (xIncrement < (nMoleculesPerSide - 1)) {
			xIncrement++; }

		if (yIncrement == nMoleculesPerSide) {
			xIncrement = 0;
			yIncrement = 0;
			zIncrement++; }
	}

	return outputCoordinates;
}

void calculateGlobalMinMax (CARTESIAN *globalSurfactanthi, CARTESIAN *globalSurfactantlo, CARTESIAN *loDimension, CARTESIAN *hiDimension, int nSurfactants)
{
	for (int i = 0; i < nSurfactants; ++i)
	{
		if (i == 0)
		{
			(*globalSurfactanthi).x = hiDimension[i].x;
			(*globalSurfactanthi).y = hiDimension[i].y;
			(*globalSurfactanthi).z = hiDimension[i].z;

			(*globalSurfactantlo).x = loDimension[i].x;
			(*globalSurfactantlo).y = loDimension[i].y;
			(*globalSurfactantlo).z = loDimension[i].z;
		}
		else
		{
			if (hiDimension[i].x > (*globalSurfactanthi).x) {
				(*globalSurfactanthi).x = hiDimension[i].x; }
			if (hiDimension[i].y > (*globalSurfactanthi).y) {
				(*globalSurfactanthi).y = hiDimension[i].y; }
			if (hiDimension[i].z > (*globalSurfactanthi).z) {
				(*globalSurfactanthi).z = hiDimension[i].z; }

			if (loDimension[i].x < (*globalSurfactantlo).x) {
				(*globalSurfactantlo).x = loDimension[i].x; }
			if (loDimension[i].y < (*globalSurfactantlo).y) {
				(*globalSurfactantlo).y = loDimension[i].y; }
			if (loDimension[i].z < (*globalSurfactantlo).z) {
				(*globalSurfactantlo).z = loDimension[i].z; }
		}
	}

	printf("\n\nMax dimensions (for all surfactants)\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nxlo: %.3f; xhi: %.3f\nylo: %.3f; yhi: %.3f\nzlo: %.3f; zhi: %.3f\n\n", 
		(*globalSurfactantlo).x, 
		(*globalSurfactanthi).x, 
		(*globalSurfactantlo).y, 
		(*globalSurfactanthi).y, 
		(*globalSurfactantlo).z, 
		(*globalSurfactanthi).z);
}

BONDS *addBonds (COORDINATES *outputCoordinates, COORDINATES **inputCoordinates, BONDS **inputBonds, SURFACTANT *inputStructures, int nSurfactants)
{
	int totalAtoms = countTotalAtoms (inputStructures, nSurfactants);

	BONDS *outputBonds;
	outputBonds = (BONDS *) malloc (totalAtoms * sizeof (BONDS));

	// Compare the molName between outputCoordinates and inputCoordinates. If the molName matches, then add bonds corresponding to the number of atoms in inputCoordinates. After that point, again compare the molName between outputCoordinates and inputCoordinates, then repeat the process multiple times.

	int currentSurfactant, surfactantBeginning = -1, surfactantEnd = -1;

	for (int i = 0; i < totalAtoms; ++i)
	{
		if (i > surfactantEnd)
		{
			for (int j = 0; j < nSurfactants; ++j)
			{
				if (strstr (outputCoordinates[i].molName, inputCoordinates[j][0].molName))
				{
					currentSurfactant = j;

					surfactantBeginning = i;
					surfactantEnd = surfactantBeginning + inputStructures[j].nAtoms - 1;
				}
			}
		}


		if (i <= surfactantEnd && i >= surfactantBeginning)
		{
			// printf("i = %d belongs to surfactant %d (%s)                   \r", i, currentSurfactant + 1, outputCoordinates[i].molName);
			// printf("current atom in this surfactant: %d                       \r", i - surfactantBeginning);
			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom1 != 0) {
				outputBonds[i].atom1 = inputBonds[currentSurfactant][i - surfactantBeginning].atom1 + surfactantBeginning; }
			else {
				outputBonds[i].atom1 = 0; }

			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom2 != 0) {
				outputBonds[i].atom2 = inputBonds[currentSurfactant][i - surfactantBeginning].atom2 + surfactantBeginning; }
			else {
				outputBonds[i].atom2 = 0; }

			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom3 != 0) {
				outputBonds[i].atom3 = inputBonds[currentSurfactant][i - surfactantBeginning].atom3 + surfactantBeginning; }
			else {
				outputBonds[i].atom3 = 0; }

			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom4 != 0) {
				outputBonds[i].atom4 = inputBonds[currentSurfactant][i - surfactantBeginning].atom4 + surfactantBeginning; }
			else {
				outputBonds[i].atom4 = 0; }

			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom5 != 0) {
				outputBonds[i].atom5 = inputBonds[currentSurfactant][i - surfactantBeginning].atom5 + surfactantBeginning; }
			else {
				outputBonds[i].atom5 = 0; }

			if (inputBonds[currentSurfactant][i - surfactantBeginning].atom6 != 0) {
				outputBonds[i].atom6 = inputBonds[currentSurfactant][i - surfactantBeginning].atom6 + surfactantBeginning; }
			else {
				outputBonds[i].atom6 = 0; }

			// printf("%2d => %3d %3d %3d %3d %3d %3d\n", i+1, outputBonds[i].atom1, outputBonds[i].atom2, outputBonds[i].atom3, outputBonds[i].atom4, outputBonds[i].atom5, outputBonds[i].atom6);

			// fflush (stdout);
			// usleep (100000);
		}
	}

	return outputBonds;
}

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
	int totalAtoms = countTotalAtoms (inputStructures, nSurfactants);
	for (int i = 0; i < totalAtoms; ++i)
	{
		printf("%3d %3d %3d %3d %3d %3d\n", outputBonds[i].atom1, outputBonds[i].atom2, outputBonds[i].atom3, outputBonds[i].atom4, outputBonds[i].atom5, outputBonds[i].atom6);
		usleep (100000);
		
	}

	// Save the above information as *.car and *.mdf files

	free (inputStructures);
	fclose (readConfig);
	return 0;
}