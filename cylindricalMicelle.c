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


int main(int argc, char const *argv[])
{

	checkArguments (argc);

	FILE *readConfig;
	readConfig = fopen (argv[1], "r");

	/*
		These are the parameters to take from input config file:
			
			1. Total number of surfactant molecules
			2. Input surfactant file name along with the total number of molecules
			3. Packing factor of each molecules

			NOTE: Any line starting with '#' can be considered as commented line
	*/

	int nSurfactants = checkNSurfactants (readConfig); // This value is coming from input config file.
	SURFACTANT *inputStructures;
	inputStructures = (SURFACTANT *) malloc (nSurfactants * sizeof (SURFACTANT));

	inputStructures = storeSurfactantInformation (inputStructures, nSurfactants, readConfig);
	inputStructures = getNAtoms (inputStructures, nSurfactants);
	inputStructures = getNBonds (inputStructures, nSurfactants);

	for (int i = 0; i < nSurfactants; ++i)
	{
		printf("%s %d %f %d %d\n", inputStructures[i].filename, inputStructures[i].nMolecules, inputStructures[i].packingFactor, inputStructures[i].nAtoms, inputStructures[i].nBonds);
	}

	free (inputStructures);
	fclose (readConfig);
	return 0;
}