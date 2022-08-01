#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

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

typedef struct surfactant
{
	int nMolecules, nAtoms, nBonds;
	float packingFactor;
	char filename[100];
} SURFACTANT;

typedef struct coordinates
{
	int sino, col5;
	float x, y, z, col9, col10;
	char atomName1[5], atomName2[5], molName[5];
} COORDINATES;

typedef struct bonds
{
	int atom1, atom2, atom3, atom4, atom5, atom6;
} BONDS;

int checkNSurfactants (FILE *readConfig)
{
	char lineString[1000];
	int nSurfactants = 0;

	while (fgets (lineString, 1000, readConfig) != NULL)
	{
		if (strstr (lineString, "nsurfactants"))
			sscanf (lineString, "%*s %d\n", &nSurfactants);
	}

	if (nSurfactants == 0)
	{
		printf("Number of surfactant molecules not initialized !\n");
		exit(1);
	}

	rewind (readConfig);
	return nSurfactants;
}

void checkArguments (int argc)
{
	if (argc == 1)
	{
		printf("\nARGUMENTS REQUIRED:\n~~~~~~~~~~~~~~~~~~~\n\n  [~] argv[0] = ./program\n  [~] argv[1] = input config file.\n\n");
		exit (1);
	}
}

SURFACTANT *storeSurfactantInformation (SURFACTANT *inputStructures, int nSurfactants, FILE *readConfig)
{
	char lineString[1000];
	int currentSurfactant = 0, MOLECULE_BEGIN = 0, MOLECULE_END = 0;

	while (fgets (lineString, 1000, readConfig) != NULL)
	{
		if (strstr (lineString, "begin surfactant"))
		{
			MOLECULE_END = 0;
			MOLECULE_BEGIN = 1;
			currentSurfactant++;
		}

		if (strstr (lineString, "end surfactant"))
		{
			MOLECULE_END = 1;
			MOLECULE_BEGIN = 0;
		}

		if (currentSurfactant > nSurfactants)
			break;

		if (MOLECULE_BEGIN == 1)
		{
			if (strstr (lineString, "filename"))
				sscanf (lineString, "%*s %s", &inputStructures[currentSurfactant - 1].filename);

			if (strstr (lineString, "nmolecules"))
				sscanf (lineString, "%*s %d", &inputStructures[currentSurfactant - 1].nMolecules);

			if (strstr (lineString, "packingfactor"))
				sscanf (lineString, "%*s %f", &inputStructures[currentSurfactant - 1].packingFactor);
		}
	}

	return inputStructures;
}

SURFACTANT *getNAtoms (SURFACTANT *inputStructures, int nSurfactants)
{
	char lineString[1000];
	int nAtoms_local = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		nAtoms_local = 0;
		FILE *inputFile;

		inputFile = fopen (inputStructures[i].filename, "r");

		while (fgets (lineString, 1000, inputFile) != NULL)
		{
			if (strstr (lineString, "ATOM") || strstr (lineString, "HETATM"))
				nAtoms_local++;
		}

		inputStructures[i].nAtoms = nAtoms_local;

		fclose (inputFile);
	}

	return inputStructures;
}

SURFACTANT *getNBonds (SURFACTANT *inputStructures, int nSurfactants)
{
	char lineString[1000];
	int nBonds_local = 0;

	for (int i = 0; i < nSurfactants; ++i)
	{
		nBonds_local = 0;
		FILE *inputFile;

		inputFile = fopen (inputStructures[i].filename, "r");

		while (fgets (lineString, 1000, inputFile) != NULL)
		{
			if (strstr (lineString, "CONECT"))
				nBonds_local++;
		}

		inputStructures[i].nBonds = nBonds_local;
		fclose (inputFile);
	}

	return inputStructures;
}

COORDINATES **assigningMemory_atoms (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{
	inputCoordinates = (COORDINATES **) malloc (nSurfactants * sizeof (COORDINATES *));

	for (int i = 0; i < nSurfactants; ++i)
	{
		printf("Assigning %d mem for surfactant %d\n", inputStructures[i].nAtoms, nSurfactants);
		inputCoordinates[i] = (COORDINATES *) malloc (inputStructures[i].nAtoms * sizeof (COORDINATES));
	}

	return inputCoordinates;
}

BONDS **assigningMemory_bonds (BONDS **inputBonds, int nSurfactants, SURFACTANT *inputStructures)
{
	inputBonds = (BONDS **) malloc (nSurfactants * sizeof (BONDS **));

	for (int i = 0; i < nSurfactants; ++i)
		inputBonds[i] = (BONDS *) malloc (inputStructures[i].nBonds * sizeof (BONDS));

	return inputBonds;
}

COORDINATES **readCoordinates (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{
	char lineString[2000];
	int sino_local;

	for (int i = 0; i < nSurfactants; ++i)
	{
		FILE *inputFile;
		inputFile = fopen (inputStructures[i].filename, "r");

		printf("Opening %s\n", inputStructures[i].filename);

		while (fgets (lineString, 2000, inputFile) != NULL)
		{
			// This is not counting for all the lines in DDAB. Check the follwing using printf statements
			if (strstr (lineString, "ATOM") || strstr (lineString, "HETATM"))
			{
				sscanf (lineString, "%*s %d %*s %*s %*d %*f %*f %*f %*f %*f %*s\n", &sino_local);
				sscanf (lineString, "%*s %d %s %s %d %f %f %f %f %f %s\n", 
					&inputCoordinates[i][sino_local].sino, 
					&inputCoordinates[i][sino_local].atomName1, 
					&inputCoordinates[i][sino_local].molName, 
					&inputCoordinates[i][sino_local].col5, 
					&inputCoordinates[i][sino_local].x, 
					&inputCoordinates[i][sino_local].y, 
					&inputCoordinates[i][sino_local].z, 
					&inputCoordinates[i][sino_local].col9, 
					&inputCoordinates[i][sino_local].col10, 
					&inputCoordinates[i][sino_local].atomName2);

				printf("%d\n", inputCoordinates[i][sino_local].sino);
			}
		}

		fclose (inputFile);
	}

	return inputCoordinates;
}

BONDS **readBonds (BONDS **inputBonds, int nSurfactants, SURFACTANT *inputStructures)
{
	char lineString[2000];
	int sino_local;

	for (int i = 0; i < nSurfactants; ++i)
	{
		FILE *inputFile;
		inputFile = fopen (inputStructures[i].filename, "r");

		// Initializing the values to '0'
		printf("%d\n", inputStructures[i].nBonds);
		fflush (stdout);

		for (int j = 0; j < inputStructures[i].nBonds; ++j)
		{
			printf("%d\n", inputBonds[i][j].atom1);
			// inputBonds[i][j].atom1 = 0;
			// inputBonds[i][j].atom2 = 0;
			// inputBonds[i][j].atom3 = 0;
			// inputBonds[i][j].atom4 = 0;
			// inputBonds[i][j].atom5 = 0;
			// inputBonds[i][j].atom6 = 0;
		}

		while (fgets (lineString, 2000, inputFile) != NULL)
		{
			if (strstr (lineString, "CONECT"))
			{
				sscanf (lineString, "%*s %d %*d %*d %*d %*d %*d %*d\n", &sino_local);
				sscanf (lineString, "%*s %*d %d %d %d %d %d %d\n", 
					&inputBonds[i][sino_local].atom1, 
					&inputBonds[i][sino_local].atom2, 
					&inputBonds[i][sino_local].atom3, 
					&inputBonds[i][sino_local].atom4, 
					&inputBonds[i][sino_local].atom5, 
					&inputBonds[i][sino_local].atom6);
			}
		}

		fclose (inputFile);
	}

	return inputBonds;
}

int main(int argc, char const *argv[])
{

	checkArguments (argc);

	FILE *readConfig;
	readConfig = fopen (argv[1], "r");

	int nSurfactants = checkNSurfactants (readConfig); // This value is coming from input config file.
	/*	
		SURFACTANT:
		~~~~~~~~~~
		int nMolecules, nAtoms, nBonds;
		float packingFactor;
		char filename[100];
	*/
	SURFACTANT *inputStructures;
	inputStructures = (SURFACTANT *) malloc (nSurfactants * sizeof (SURFACTANT));

	inputStructures = storeSurfactantInformation (inputStructures, nSurfactants, readConfig);
	inputStructures = getNAtoms (inputStructures, nSurfactants);
	inputStructures = getNBonds (inputStructures, nSurfactants);

	// Create variables to store atom coordinates and bond connectivity information
	COORDINATES **inputCoordinates;
	BONDS **inputBonds;

	inputCoordinates = assigningMemory_atoms (inputCoordinates, nSurfactants, inputStructures);
	inputBonds = assigningMemory_bonds (inputBonds, nSurfactants, inputStructures);

	/*
		COORDINATES:
		~~~~~~~~~~~
		int sino, col5;
		float x, y, z, col9, col10;
		char atomName1[5], atomName2[5], molName[5];

		BONDS:
		~~~~~
		int atom1, atom2, atom3, atom4, atom5, atom6;
	*/

	inputCoordinates = readCoordinates (inputCoordinates, nSurfactants, inputStructures);
	inputBonds = readBonds (inputBonds, nSurfactants, inputStructures);

	free (inputStructures);
	fclose (readConfig);
	return 0;
}