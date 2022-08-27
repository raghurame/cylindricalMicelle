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

void readData (FILE *input, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers, DATAFILE_INFO *datafile, BOUNDS *datafileBoundary)
{
	printf("Reading LAMMPS data file...\n");

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

	printf("\nPrinting boundary information from data file:\n\n  xlo: %f; xhi: %f\n  ylo: %f; yhi: %f\n  zlo: %f; zhi: %f\n", (*datafileBoundary).xlo, (*datafileBoundary).xhi, (*datafileBoundary).ylo, (*datafileBoundary).yhi, (*datafileBoundary).zlo, (*datafileBoundary).zhi);
	printf("\nFrom input data file:\n\n  nAtoms: %d\n  nBonds: %d\n  nAngles: %d\n  nDihedrals: %d\n  nImpropers: %d\n\n", (*datafile).nAtoms, (*datafile).nBonds, (*datafile).nAngles, (*datafile).nDihedrals, (*datafile).nImpropers);
	rewind (input);
}

DUMP *readLastDumpFrame (char *pipeString, int nAtoms)
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

BOUNDS readDumpfileBoundary (const char *pipe_dumpBoundary)
{
	FILE *readingBoundary;
	readingBoundary = popen (pipe_dumpBoundary, "r");

	char lineString[2000];

	BOUNDS dumpfileBoundary;

	fgets (lineString, 2000, readingBoundary);
	sscanf (lineString, "%f %f\n", &dumpfileBoundary.xlo, &dumpfileBoundary.xhi);
	fgets (lineString, 2000, readingBoundary);
	sscanf (lineString, "%f %f\n", &dumpfileBoundary.ylo, &dumpfileBoundary.yhi);
	fgets (lineString, 2000, readingBoundary);
	sscanf (lineString, "%f %f\n", &dumpfileBoundary.zlo, &dumpfileBoundary.zhi);

	printf("Printing Boundary information from dump file:\n\nxlo: %f;xhi: %f;\nylo: %f; yhi: %f;\nzlo: %f; zhi: %f;\n\n", dumpfileBoundary.xlo, dumpfileBoundary.xhi, dumpfileBoundary.ylo, dumpfileBoundary.yhi, dumpfileBoundary.zlo, dumpfileBoundary.zhi);

	pclose (readingBoundary);
	return dumpfileBoundary;
}

int calculateNWater (BOUNDS dumpfileBoundary)
{
	float xLength = (dumpfileBoundary.xhi - dumpfileBoundary.xlo), yLength = (dumpfileBoundary.yhi - dumpfileBoundary.ylo), zLength = (dumpfileBoundary.zhi - dumpfileBoundary.zlo), avogadroNumber = 6.023, waterDensity = 1.0, molarMass = 18;
	float nWater_max = (waterDensity * xLength * yLength * zLength * avogadroNumber * 0.1 / molarMass);
	return nWater_max;
}

float deg2rad (float degrees) {
	return degrees * M_PI / 180.0; }

void findAttachedHydrogens (float ox, float oy, float oz, float *h1x, float *h1y, float *h1z, float *h2x, float *h2y, float *h2z)
{
	// HOH angle is 113.24. Consider O at the origin and the two Hs placed along the XY plane (facing along the positive X axis)
	// Each OH bond makes 56.62 degrees with the X axis. 
	// So, for H1, the X coordinate is cos (56.62) and the Y coordinate is sin (56.62)
	// For H2, the X coordinate is cos (56.62) and the Y coordinate  is -sin (56.62)
	(*h1x) = ox + cosf (deg2rad (56.62));
	(*h1y) = oy + sinf (deg2rad (56.62));
	(*h1z) = oz;
	(*h2x) = ox + cosf (deg2rad (56.62));
	(*h2y) = oy -sinf (deg2rad (56.62));
	(*h2z) = oz;
}

DATA_ATOMS *populateWater (DATA_ATOMS *atomsWater, int nWater, BOUNDS dumpfileBoundary, DATA_ATOMS *atoms, DATAFILE_INFO datafileInfo)
{
	int nBins_x = (int) floor (cbrt (nWater)), nBins_y = (int) floor (cbrt (nWater)), nBins_z = (int) floor (cbrt (nWater));
	float distSeparation_x = (dumpfileBoundary.xhi - dumpfileBoundary.xlo) / nBins_x, distSeparation_y = (dumpfileBoundary.yhi - dumpfileBoundary.ylo) / nBins_y, distSeparation_z = (dumpfileBoundary.zhi - dumpfileBoundary.zlo) / nBins_z;
	float distance_O_mol, distance_H1_mol, distance_H2_mol;
	int currentWaterAtom = 0, isOverlap = 0;

	// Distributing water evenly in cartesian space
	for (int i = 0; i < nBins_x; ++i)
	{
		for (int j = 0; j < nBins_y; ++j)
		{
			for (int k = 0; k < nBins_z; ++k)
			{
				atomsWater[currentWaterAtom].x = dumpfileBoundary.xlo + ((i + 1) * distSeparation_x) - (distSeparation_x / 2);
				atomsWater[currentWaterAtom].y = dumpfileBoundary.ylo + ((j + 1) * distSeparation_y) - (distSeparation_y / 2);
				atomsWater[currentWaterAtom].z = dumpfileBoundary.zlo + ((k + 1) * distSeparation_z) - (distSeparation_z / 2);

				findAttachedHydrogens (atomsWater[currentWaterAtom].x, atomsWater[currentWaterAtom].y, atomsWater[currentWaterAtom].z, &atomsWater[currentWaterAtom + 1].x, &atomsWater[currentWaterAtom + 1].y, &atomsWater[currentWaterAtom + 1].z, &atomsWater[currentWaterAtom + 2].x, &atomsWater[currentWaterAtom + 2].y, &atomsWater[currentWaterAtom + 2].z);

				// Resetting the isOverlap variable before checking the distances
				isOverlap = 0;

				for (int i = 0; i < datafileInfo.nAtoms; ++i)
				{
					distance_O_mol = sqrt (pow (atomsWater[currentWaterAtom].x - atoms[i].x, 2) + pow (atomsWater[currentWaterAtom].y - atoms[i].y, 2) + pow (atomsWater[currentWaterAtom].z - atoms[i].z, 2));
					distance_H1_mol = sqrt (pow (atomsWater[currentWaterAtom + 1].x - atoms[i].x, 2) + pow (atomsWater[currentWaterAtom + 1].y - atoms[i].y, 2) + pow (atomsWater[currentWaterAtom + 1].z - atoms[i].z, 2));
					distance_H2_mol = sqrt (pow (atomsWater[currentWaterAtom + 2].x - atoms[i].x, 2) + pow (atomsWater[currentWaterAtom + 2].y - atoms[i].y, 2) + pow (atomsWater[currentWaterAtom + 2].z - atoms[i].z, 2));

					if ((distance_O_mol < 3) || (distance_H1_mol < 3) || (distance_H2_mol < 3)) {
						isOverlap = 1; }
				}

				if (isOverlap == 0) {
					currentWaterAtom += 3;
					printf("Adding water... %d/%d                           \r", (int) (currentWaterAtom / 3), nWater);
					fflush (stdout); }
			}
		}
	}

	printf("\n\nMax number of water that can be added (based on overall simulation volume): %d\nNumber of water molecules added after checking for overlaps: %d\n\n", nWater, (int) floor (currentWaterAtom / 3));

	return atomsWater;
}

int main(int argc, char const *argv[])
{
	FILE *inputData, *inputDump, *output;
	inputData = fopen (argv[1], "r");
	inputDump = fopen (argv[2], "r");
	output = fopen (argv[3], "w");

	int nAtoms = getNatoms (inputDump);

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafileInfo;

	BOUNDS datafileBoundary, dumpfileBoundary;

	readData (inputData, &atoms, &bonds, &angles, &dihedrals, &impropers, &datafileInfo, &datafileBoundary);

	DUMP *traj;
	traj = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	char *pipe_lastframe;
	pipe_lastframe = (char *) malloc (50 * sizeof (char));
	snprintf (pipe_lastframe, 50, "tail -%d %s", nAtoms, argv[2]);
	traj = readLastDumpFrame (pipe_lastframe, nAtoms);

	char *pipe_dumpBoundary;
	pipe_dumpBoundary = (char *) malloc (50 * sizeof (char));
	snprintf (pipe_dumpBoundary, 50, "tail -%d %s | head -3", (nAtoms + 4), argv[2]);
	dumpfileBoundary = readDumpfileBoundary (pipe_dumpBoundary);

	// Replacing the coordinates in data file with the coordinates obtained from trajectory file
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].x = traj[i].x;
		atoms[i].y = traj[i].y;
		atoms[i].z = traj[i].z;
	}

	// Add SPC/Fw water molecules within the dumpfileBoundary
	// Adding new atoms, bonds, and angles for new water molecules
	DATA_ATOMS *atomsWater;
	DATA_BONDS *bondsWater;
	DATA_ANGLES *anglesWater;
	DATA_DIHEDRALS *dihedralsWater;
	DATA_IMPROPERS *impropersWater;

	int nWater = calculateNWater (dumpfileBoundary);

	// There are three atoms, 2 bonds and 1 angle for every water molecule
	atomsWater = (DATA_ATOMS *) malloc (nWater * 3 * sizeof (DATA_ATOMS));
	bondsWater = (DATA_BONDS *) malloc (nWater * 2 * sizeof (DATA_BONDS));
	anglesWater = (DATA_ANGLES *) malloc (nWater * sizeof (DATA_ANGLES));

	atomsWater = populateWater (atomsWater, nWater, dumpfileBoundary, atoms, datafileInfo);

	// Recalculate the simulation box size based on the added water molecules
	// because some Hs can protrude sligtly outside the simulation box

	// Then add bonds and angles

	// Print the final data file (with masses for water)

	return 0;
}