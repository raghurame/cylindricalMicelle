#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct dump
{
	int id, type, ix, iy, iz;
	float x, y, z, xs, ys, zs;
} DUMP;

typedef struct bounds
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
} BOUNDS;

typedef struct mass
{
	int atomType;
	float mass;
} ATOMIC_MASS;

int getNatoms (FILE *inputDump)
{
	char lineString[2000];
	int lineNumber = 0, natoms;

	while (fgets (lineString, 2000, inputDump) != NULL)
	{
		lineNumber++;

		if (lineNumber == 4)
		{
			sscanf (lineString, "%d", &natoms);
			rewind (inputDump);
			return natoms;
		}
	}
	return 0;
}

DUMP *readLastDumpFrame_full (char *pipeString, int nAtoms)
{
	FILE *input;
	input = popen (pipeString, "r");

	DUMP *traj, *traj_temp, com, max;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	traj_temp = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	int lineCount = 0;

	char lineString[1000];

	while (fgets (lineString, 1000, input) != NULL)
	{
		sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", 
			&traj_temp[lineCount].id, 
			&traj_temp[lineCount].type, 
			&traj_temp[lineCount].x, 
			&traj_temp[lineCount].y, 
			&traj_temp[lineCount].z, 
			&traj_temp[lineCount].xs, 
			&traj_temp[lineCount].ys, 
			&traj_temp[lineCount].zs, 
			&traj_temp[lineCount].ix, 
			&traj_temp[lineCount].iy, 
			&traj_temp[lineCount].iz);
		lineCount++;
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nAtoms; ++j)
		{
			if (traj_temp[j].id == i + 1)
			{
				traj[i].id = traj_temp[j].id;
				traj[i].type = traj_temp[j].type;
				traj[i].x = traj_temp[j].x;
				traj[i].y = traj_temp[j].y;
				traj[i].z = traj_temp[j].z;
				traj[i].xs = traj_temp[j].xs;
				traj[i].ys = traj_temp[j].ys;
				traj[i].zs = traj_temp[j].zs;
				traj[i].ix = traj_temp[j].ix;
				traj[i].iy = traj_temp[j].iy;
				traj[i].iz = traj_temp[j].iz;
			}
		}
	}

	pclose (input);
	return traj;
}

DUMP *readLastDumpFrame_minimal (char *pipeString, int nAtoms)
{
	FILE *input;
	input = popen (pipeString, "r");

	DUMP *traj, *traj_temp, com, max;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	traj_temp = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	int lineCount = 0;

	char lineString[1000];

	while (fgets (lineString, 1000, input) != NULL)
	{
		sscanf (lineString, "%d %d %f %f %f\n", 
			&traj_temp[lineCount].id, 
			&traj_temp[lineCount].type, 
			&traj_temp[lineCount].x, 
			&traj_temp[lineCount].y, 
			&traj_temp[lineCount].z);
		lineCount++;
	}

	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nAtoms; ++j)
		{
			if (traj_temp[j].id == i + 1)
			{
				traj[i].id = traj_temp[j].id;
				traj[i].type = traj_temp[j].type;
				traj[i].x = traj_temp[j].x;
				traj[i].y = traj_temp[j].y;
				traj[i].z = traj_temp[j].z;
			}
		}
	}

	pclose (input);
	return traj;
}

BOUNDS calculateDumpBoundary (DUMP *goldCoords, int nGoldAtoms, BOUNDS goldBoundary)
{
	goldBoundary.xlo = 0;
	goldBoundary.xhi = 0;
	goldBoundary.ylo = 0;
	goldBoundary.yhi = 0;
	goldBoundary.zlo = 0;
	goldBoundary.zhi = 0;

	for (int i = 0; i < nGoldAtoms; ++i)
	{
		if (i == 0)
		{
			goldBoundary.xlo = goldCoords[i].x; goldBoundary.xhi = goldCoords[i].x;
			goldBoundary.ylo = goldCoords[i].y; goldBoundary.yhi = goldCoords[i].y;
			goldBoundary.zlo = goldCoords[i].z; goldBoundary.zhi = goldCoords[i].z;
		}
		else
		{
			if (goldCoords[i].x < goldBoundary.xlo) {
				goldBoundary.xlo = goldCoords[i].x; }
			else if (goldCoords[i].x > goldBoundary.xhi) {
				goldBoundary.xhi = goldCoords[i].x; }

			if (goldCoords[i].y < goldBoundary.ylo) {
				goldBoundary.ylo = goldCoords[i].y; }
			else if (goldCoords[i].y > goldBoundary.yhi) {
				goldBoundary.yhi = goldCoords[i].y; }

			if (goldCoords[i].z < goldBoundary.zlo) {
				goldBoundary.zlo = goldCoords[i].z; }
			else if (goldCoords[i].z > goldBoundary.zhi) {
				goldBoundary.zhi = goldCoords[i].z; }
		}
	}

	return goldBoundary;
}

DUMP computeCenterDump (DUMP center, DUMP *coords, int nAtoms)
{
	center.x = 0; center.y = 0; center.z = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		center.x += coords[i].x;
		center.y += coords[i].y;
		center.z += coords[i].z;
	}

	center.x /= nAtoms;
	center.y /= nAtoms;
	center.z /= nAtoms;

	return center;
}

int checkNAtoms (DUMP *coords, int nAtoms, int atomType)
{
	int nAtoms_ofThatAtomType = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (coords[i].type == atomType)
		{
			nAtoms_ofThatAtomType++;
		}
	}

	return nAtoms_ofThatAtomType;
}

int main(int argc, char const *argv[])
{
	if (argc == 1) 	{
		(void)printf("\nARGUMENTS:\n~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = main lammps dump filename\n{~} argv[2] = atom type to be replaced\n{~} argv[3] = substitute lammps dump filename\n{~} argv[4] = substitute atom type\n{~} argv[5] = output filename.\n\n");
		exit (1); }

	FILE *inputDump_file, *substituteCoords_file, *output_file;
	inputDump_file = fopen (argv[1], "r");
	int atomTypeMain = atoi (argv[2]);
	substituteCoords_file = fopen (argv[3], "r");
	int atomTypeSubstitute = atoi (argv[4]);
	output_file = fopen (argv[5], "w");

	int nMainAtoms = getNatoms (inputDump_file), nSubstituteAtoms = getNatoms (substituteCoords_file);

	DUMP *inputDump_traj, *substituteAtoms_traj;
	inputDump_traj = (DUMP *) malloc (nMainAtoms * sizeof (DUMP));
	substituteAtoms_traj = (DUMP *) malloc (nSubstituteAtoms * sizeof (DUMP));

	// Reading the last timeframe of main dump file
	char *pipe_lastframe;
	pipe_lastframe = (char *) malloc (100 * sizeof (char));
	snprintf (pipe_lastframe, 100, "tail -%d %s", nMainAtoms, argv[1]);
	inputDump_traj = readLastDumpFrame_full (pipe_lastframe, nMainAtoms);

	// Reading the last timeframe of substitute dump file
	snprintf (pipe_lastframe, 100, "tail -%d %s", nSubstituteAtoms, argv[3]);
	substituteAtoms_traj = readLastDumpFrame_minimal (pipe_lastframe, nSubstituteAtoms);

	// Replace a particular atom type in inputDump_file with that present in substituteCoords_file
	// The number of atoms in both cases must be same, so it can be a one-to-one substitute
	int nMainAtoms_ofThatAtomType = checkNAtoms (inputDump_traj, nMainAtoms, atomTypeMain), nSubstituteAtoms_ofThatAtomType = checkNAtoms (substituteAtoms_traj, nSubstituteAtoms, atomTypeSubstitute);

	printf("%d %d\n", nMainAtoms_ofThatAtomType, nSubstituteAtoms_ofThatAtomType);

	fclose (inputDump_file);
	fclose (substituteCoords_file);
	fclose (output_file);
	return 0;
}