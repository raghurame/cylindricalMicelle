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
			printf("mol %d: %.2f, %.2f --> ", i, inputCoordinates[i][j].x, inputCoordinates[i][j].y);
			inputCoordinates[i][j].x = (inputCoordinates[i][j].x * cosf (theta)) - (inputCoordinates[i][j].y * sinf (theta));
			inputCoordinates[i][j].y = (inputCoordinates[i][j].x * sinf (theta)) - (inputCoordinates[i][j].y * cosf (theta));
			printf("%.2f, %.2f\n", inputCoordinates[i][j].x, inputCoordinates[i][j].y);
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

	/*
	Pack them in some lattice for now (later, pack them in micelle structure)
	While packing them in lattice structure, don't worry about packing factor. Just maintain some tolerance
	Later, while packing the moleules in a micelle structure, think about implementing soft repulsive potential and energy minimization.
	*/

	FARTHESTPOINTS *inputStructures_farPoints;
	inputStructures_farPoints = (FARTHESTPOINTS *) malloc (nSurfactants * sizeof (FARTHESTPOINTS));

	inputStructures_farPoints = calculateFarPoints (inputStructures_farPoints, inputCoordinates, nSurfactants, inputStructures);

	inputCoordinates = orientSurfactants (inputCoordinates, nSurfactants, inputStructures, inputStructures_farPoints);

	// Find the longest dimension on X, Y, and Z

	// Replicate the molecule (coordinates and bonds)

	// Save the above information as *.car and *.mdf files

	free (inputStructures);
	fclose (readConfig);
	return 0;
}