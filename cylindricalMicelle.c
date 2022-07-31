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
	int nMolecules;
	float packingFactor;
	char filename[100];
} SURFACTANT;

int checkNSurfactants (FILE *readConfig)
{
	
}

int main(int argc, char const *argv[])
{
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

	fclose (readConfig);
	return 0;
}