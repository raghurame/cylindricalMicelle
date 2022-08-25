#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

/*

This code takes dump file as input, adds water and creates a data file.

Number of arguments:
~~~~~~~~~~~~~~~~~~~~

argv[0] = ./program
argv[1] = input data filename
argv[2] = input dump filename
argv[3] = output data filename

*/

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct dump
{
	int id, type, ix, iy, iz;
	float x, y, z, xs, ys, zs;
} DUMP;

typedef struct bounds
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
} BOUNDS;

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers, DATAFILE_INFO *datafile, BOUNDS *datafileBoundary)
{
	printf("Reading LAMMPS data file...\n");
	FILE *input;
	input = fopen (dataFileName, "r");

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	// DATAFILE_INFO datafile;
	(*datafile).nAtoms = -1;
	(*datafile).nBonds = -1;
	(*datafile).nAngles = -1;
	(*datafile).nDihedrals = -1;
	(*datafile).nImpropers = -1;

	char lineString[1000];

	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &(*datafile).nAtoms);
			(*atoms) = (DATA_ATOMS *) malloc ((*datafile).nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &(*datafile).nBonds);
			(*bonds) = (DATA_BONDS *) malloc ((*datafile).nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &(*datafile).nAngles);
			(*angles) = (DATA_ANGLES *) malloc ((*datafile).nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &(*datafile).nDihedrals);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc ((*datafile).nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &(*datafile).nImpropers);
			(*impropers) = (DATA_IMPROPERS *) malloc ((*datafile).nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
			sscanf (lineString, "%d \n", &(*datafile).nAtomTypes);

		if (strstr (lineString, "bond types"))
			sscanf (lineString, "%d \n", &(*datafile).nBondTypes);

		if (strstr (lineString, "angle types"))
			sscanf (lineString, "%d \n", &(*datafile).nAngleTypes);

		if (strstr (lineString, "dihedral types"))
			sscanf (lineString, "%d \n", &(*datafile).nDihedralTypes);

		if (strstr (lineString, "improper types"))
			sscanf (lineString, "%d \n", &(*datafile).nImproperTypes);

		if (((*datafile).nAtoms >= 0) && ((*datafile).nBonds >= 0) && ((*datafile).nAngles >= 0) && ((*datafile).nDihedrals >= 0) && ((*datafile).nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "xlo") && strstr (lineString, "xhi")) {
			sscanf (lineString, "%f %f \n", &(*datafileBoundary).xlo, &(*datafileBoundary).xhi); }

		if (strstr (lineString, "ylo") && strstr (lineString, "yhi")) {
			sscanf (lineString, "%f %f \n", &(*datafileBoundary).ylo, &(*datafileBoundary).yhi); }

		if (strstr (lineString, "zlo") && strstr (lineString, "zhi")) {
			sscanf (lineString, "%f %f \n", &(*datafileBoundary).zlo, &(*datafileBoundary).zhi); }

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == (*datafile).nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == (*datafile).nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == (*datafile).nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == (*datafile).nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == (*datafile).nImpropers)
				isImproperLine = 0;
		}
	}

	// printf("checking\n");
	// printf("Printing boundary information from data file:\nxlo: %f; xhi: %f\nylo: %f; yhi: %f\nzlo: %f; zhi: %f\n", (*datafileBoundary).xlo, (*datafileBoundary).xhi, (*datafileBoundary).ylo, (*datafileBoundary).yhi, (*datafileBoundary).zlo, (*datafileBoundary).zhi);
	// printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", (*datafile).nAtoms, (*datafile).nBonds, (*datafile).nAngles, (*datafile).nDihedrals, (*datafile).nImpropers);
	// fflush (stdout);
}

DUMP *readLastDumpFrame (const char *pipeString, int nAtoms)
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
		sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", &traj_temp[lineCount].id, &traj_temp[lineCount].type, &traj_temp[lineCount].x, &traj_temp[lineCount].y, &traj_temp[lineCount].z, &traj_temp[lineCount].xs, &traj_temp[lineCount].ys, &traj_temp[lineCount].zs, &traj_temp[lineCount].ix, &traj_temp[lineCount].iy, &traj_temp[lineCount].iz);
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

int getNatoms (const char *inputFileName)
{
	FILE *read = fopen (inputFileName, "r");
	char lineString[2000];
	int lineNumber = 0, natoms;

	while (fgets (lineString, 2000, read) != NULL)
	{
		lineNumber++;

		if (lineNumber == 4)
		{
			sscanf (lineString, "%d", &natoms);
			return natoms;
		}
	}
	return 0;
}

int main(int argc, char const *argv[])
{
	FILE *inputData, *inputDump, *output;
	int nAtoms = 5;// getNatoms (argv[2]);

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafileInfo;

	BOUNDS datafileBoundary, dumpfileBoundary;

	printf("%s\n", argv[1]);
	printf("%s\n", argv[2]);
	// readData (argv[1], &atoms, &bonds, &angles, &dihedrals, &impropers, &datafileInfo, &datafileBoundary);

	DUMP *traj;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	// traj = readLastDumpFrame (argv[2], nAtoms);

	return 0;
}