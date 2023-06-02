#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "structs.h"
#include "readInputs.h"

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

	printf("Number of surfactants: %d\n", nSurfactants);

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
			printf("\n");
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

			printf("Input filename: %s\n", inputStructures[currentSurfactant - 1].filename);
			printf("Number of molecules: %d\n", inputStructures[currentSurfactant - 1].nMolecules);
		}
	}

	return inputStructures;
}

SURFACTANT *getNAtoms (SURFACTANT *inputStructures, int nSurfactants)
{
	char lineString[1000];
	int nAtoms_local = 0;
	printf("Reading atoms from input pdb file...\n");

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
	printf("Reading bonds from input pdb file...\n");
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

COORDINATES **readCoordinates (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures)
{
	char lineString[2000];
	int sino_local = 0;

	inputCoordinates = (COORDINATES **) malloc (nSurfactants * sizeof (COORDINATES *));

	for (int i = 0; i < nSurfactants; ++i)
	{
		printf("Assigning %d mem (COORDINATES) for %s\n", inputStructures[i].nAtoms, inputStructures[i].filename);
		inputCoordinates[i] = (COORDINATES *) malloc (inputStructures[i].nAtoms * sizeof (COORDINATES));
	}

	for (int i = 0; i < nSurfactants; ++i)
	{
		sino_local = 0;
		FILE *inputFile;
		inputFile = fopen (inputStructures[i].filename, "r");

		char atomName1_temp[5], molName_temp[5], atomName2_temp[5];
		int col5_temp;
		float x_temp, y_temp, z_temp, col9_temp, col10_temp;

		printf("Opening %s\n", inputStructures[i].filename);

		while (fgets (lineString, 2000, inputFile) != NULL)
		{
			if (strstr (lineString, "ATOM") || strstr (lineString, "HETATM"))
			{
				sscanf (lineString, "%*s %*d %s %s %d %f %f %f %f %f %s\n", 
					&atomName1_temp, 
					&molName_temp, 
					&col5_temp, 
					&x_temp, 
					&y_temp, 
					&z_temp, 
					&col9_temp, 
					&col10_temp, 
					&atomName2_temp);

				strcpy (inputCoordinates[i][sino_local].atomName1, atomName1_temp);
				strcpy (inputCoordinates[i][sino_local].molName, molName_temp);
				inputCoordinates[i][sino_local].col5  = col5_temp;
				inputCoordinates[i][sino_local].x  = x_temp;
				inputCoordinates[i][sino_local].y  = y_temp;
				inputCoordinates[i][sino_local].z  = z_temp;
				inputCoordinates[i][sino_local].col9  = col9_temp;
				inputCoordinates[i][sino_local].col10  = col10_temp;
				strcpy (inputCoordinates[i][sino_local].atomName2, atomName2_temp);

				sino_local++;
			}
		}

		fclose (inputFile);
	}

	return inputCoordinates;
}

BONDS **readBonds (BONDS **inputBonds, int nSurfactants, SURFACTANT *inputStructures)
{
	char lineString[2000];
	int sino_local = 0;

	inputBonds = (BONDS **) malloc (nSurfactants * sizeof (BONDS *));

	for (int i = 0; i < nSurfactants; ++i)
	{
		printf("Assigning %d mem (BONDS) for %s\n", inputStructures[i].nBonds, inputStructures[i].filename);
		inputBonds[i] = (BONDS *) malloc (inputStructures[i].nBonds * sizeof (BONDS));
	}

	for (int i = 0; i < nSurfactants; ++i)
	{
		for (int j = 0; j < inputStructures[i].nBonds; ++j)
		{
			inputBonds[i][j].atom1 = 0;
			inputBonds[i][j].atom2 = 0;
			inputBonds[i][j].atom3 = 0;
			inputBonds[i][j].atom4 = 0;
			inputBonds[i][j].atom5 = 0;
			inputBonds[i][j].atom6 = 0;
		}
	}

	for (int i = 0; i < nSurfactants; ++i)
	{
		sino_local = 0;
		FILE *inputFile;
		inputFile = fopen (inputStructures[i].filename, "r");

		char atomName1_temp[5], molName_temp[5], atomName2_temp[5];
		int col5_temp;
		float x_temp, y_temp, z_temp, col9_temp, col10_temp;

		printf("Opening %s\n", inputStructures[i].filename);

		while (fgets (lineString, 2000, inputFile) != NULL)
		{
			if (strstr (lineString, "CONECT"))
			{
				sscanf (lineString, "%*s %d %d %d %d %d %d\n", 
					&inputBonds[i][sino_local].atom1,
					&inputBonds[i][sino_local].atom2,
					&inputBonds[i][sino_local].atom3,
					&inputBonds[i][sino_local].atom4,
					&inputBonds[i][sino_local].atom5,
					&inputBonds[i][sino_local].atom6);

				sino_local++;
			}
		}

		fclose (inputFile);
	}

	return inputBonds;
}