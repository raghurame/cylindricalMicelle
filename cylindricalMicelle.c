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

CENTEROFMASS *computeCOM (COORDINATES *inputCoordinates, int nSurfactants, SURFACTANT inputStructures)
{
	CENTEROFMASS *comForMolecules;
	comForMolecules = (CENTEROFMASS *) malloc (nSurfactants * sizeof (CENTEROFMASS));

	

	return comForMolecules;
}

COORDINATES **orientSurfactants (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, FARTHESTPOINTS *inputStructures_farPoints)
{
	// This function should do the following
	// Center the surfactant molecule. Keep one point as (0, 0, 0)
	// Find the angle between the vector and X axis in XY/XZ planes
	// Separately rotate the vector on XY/XZ planes

	CENTEROFMASS *comForMolecules;
	comForMolecules = computeCOM (inputCoordinates, nSurfactants, inputStructures); // assuming all atoms have same mass

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
	Find two farthest points in every molecule
	Orient the chain along X axis using these two points.
	Pack them in some lattice for now (later, pack them in micelle structure)
	While packing them in lattice structure, don't worry about packing factor. Just maintain some tolerance
	Later, while packing the moleules in a micelle structure, think about implementing soft repulsive potential and energy minimization.
	*/

	FARTHESTPOINTS *inputStructures_farPoints;
	inputStructures_farPoints = (FARTHESTPOINTS *) malloc (nSurfactants * sizeof (FARTHESTPOINTS));

	inputStructures_farPoints = calculateFarPoints (inputStructures_farPoints, inputCoordinates, nSurfactants, inputStructures);

	inputCoordinates = orientSurfactants (inputCoordinates, nSurfactants, inputStructures, inputStructures_farPoints);

	free (inputStructures);
	fclose (readConfig);
	return 0;
}