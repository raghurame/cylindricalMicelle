#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

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
	int maxMolType;
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

typedef struct mass
{
	int atomType;
	float mass;
} ATOMIC_MASS;

void readData (FILE *input, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers, DATAFILE_INFO *datafile, BOUNDS *datafileBoundary, ATOMIC_MASS **mass)
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

	(*datafile).maxMolType = 0;

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

		if (strstr (lineString, "atom types")) {
			sscanf (lineString, "%d \n", &(*datafile).nAtomTypes);
			(*mass) = (ATOMIC_MASS *) malloc ((*datafile).nAtomTypes * sizeof (ATOMIC_MASS)); }

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

		if (strstr (lineString, "Masses"))
		{
			fgets (lineString, 1000, input);
			for (int i = 0; i < (*datafile).nAtomTypes; ++i)
			{
				fgets (lineString, 1000, input);
				sscanf (lineString, "%d %f", &(*mass)[i].atomType, &(*mass)[i].mass);
			}
		}

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
			if ((*atoms)[nAtomLine].molType > (*datafile).maxMolType) {
				(*datafile).maxMolType = (*atoms)[nAtomLine].molType; }
			nAtomLine++;
			if (nAtomLine == (*datafile).nAtoms) {
				isAtomLine = 0; }
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
	printf("\nMasses:\n\n");
	for (int i = 0; i < (*datafile).nAtomTypes; ++i) {
		printf("%d %f\n", (*mass)[i].atomType, (*mass)[i].mass); }
	printf("\n");
	printf("Max mol type present in the data file: %d\n\n", (*datafile).maxMolType);
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

void printDatafile (DATA_ATOMS *atoms, DATA_BONDS *bonds, DATA_ANGLES *angles, DATA_DIHEDRALS *dihedrals, DATA_IMPROPERS *impropers, DATAFILE_INFO datafileInfo, BOUNDS datafileBoundary, ATOMIC_MASS *mass, DUMP *goldCoords, int nGoldAtoms, BOUNDS overallBoundary, FILE *outputData, FILE *outputXYZ)
{
	// Printing the header information
	fprintf(outputData, "%s\n\n", "LAMMPS data file containing the original structures and the polarizable gold");
	fprintf(outputData, "%d atoms\n%d bonds\n%d angles\n%d dihedrals\n%d impropers\n\n%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nMasses\n\n", datafileInfo.nAtoms + (nGoldAtoms * 2), datafileInfo.nBonds + nGoldAtoms, datafileInfo.nAngles, datafileInfo.nDihedrals, datafileInfo.nImpropers, datafileInfo.nAtomTypes + 2, datafileInfo.nBondTypes + 1, datafileInfo.nAngleTypes, datafileInfo.nDihedralTypes, overallBoundary.xlo, overallBoundary.xhi, overallBoundary.ylo, overallBoundary.yhi, overallBoundary.zlo, overallBoundary.zhi);

	fprintf(outputXYZ, "%d\n", datafileInfo.nAtoms + (nGoldAtoms * 2));
	fprintf(outputXYZ, "%s\n", "Dummy comment line");

	// Printing masses
	for (int i = 0; i < datafileInfo.nAtomTypes; ++i) {
		fprintf(outputData, "%d %f\n", mass[i].atomType, mass[i].mass); }
	fprintf(outputData, "%d %f\n", datafileInfo.nAtomTypes + 1, 195.96); // For gold atom
	fprintf(outputData, "%d %f\n", datafileInfo.nAtomTypes + 2, 1.0); // For electron cloud

	// Printing orignal atoms
	fprintf(outputData, "\nAtoms\n\n");
	for (int i = 0; i < datafileInfo.nAtoms; ++i) {
		fprintf(outputXYZ, "C %f %f %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
		fprintf(outputData, "%d %d %d %f %f %f %f\n", atoms[i].id, atoms[i].molType, atoms[i].atomType, atoms[i].charge, atoms[i].x, atoms[i].y, atoms[i].z); }

	// Printing gold atoms
	for (int i = 0; i < nGoldAtoms; ++i) {
		fprintf(outputXYZ, "C %f %f %f\n", goldCoords[i].x, goldCoords[i].y, goldCoords[i].z);
		fprintf(outputData, "%d %d %d 1.0 %f %f %f\n", (datafileInfo.nAtoms + i + 1), (datafileInfo.maxMolType + 1), (datafileInfo.nAtomTypes + 1), goldCoords[i].x, goldCoords[i].y, goldCoords[i].z); }

	// Printing the electron cloud
	for (int i = 0; i < nGoldAtoms; ++i) {
		fprintf(outputXYZ, "C %f %f %f\n", goldCoords[i].x, goldCoords[i].y, goldCoords[i].z);
		fprintf(outputData, "%d %d %d -1.0 %f %f %f\n", (datafileInfo.nAtoms + nGoldAtoms + i + 1), (datafileInfo.maxMolType + 2), (datafileInfo.nAtomTypes + 2), goldCoords[i].x, goldCoords[i].y, goldCoords[i].z); }

	// Printing the original bonds
	fprintf(outputData, "\nBonds\n\n");
	for (int i = 0; i < datafileInfo.nBonds; ++i) {
		fprintf(outputData, "%d %d %d %d\n", bonds[i].id, bonds[i].bondType, bonds[i].atom1, bonds[i].atom2); }

	// Printing the bonds in polarizable gold
	for (int i = 0; i < nGoldAtoms; ++i) {
		fprintf(outputData, "%d %d %d %d\n", (datafileInfo.nBonds + 1 + i), (datafileInfo.nBondTypes + 1), (datafileInfo.nAtoms + i + 1), (datafileInfo.nAtoms + nGoldAtoms + i + 1)); }

	// Printing the original angles
	fprintf(outputData, "\nAngles\n\n");
	for (int i = 0; i < datafileInfo.nAngles; ++i) {
		fprintf(outputData, "%d %d %d %d %d\n", angles[i].id, angles[i].angleType, angles[i].atom1, angles[i].atom2, angles[i].atom3); }

	// Printing the original dihedrals
	fprintf(outputData, "\nDihedrals\n\n");
	for (int i = 0; i < datafileInfo.nDihedrals; ++i) {
		fprintf(outputData, "%d %d %d %d %d %d\n", dihedrals[i].id, dihedrals[i].dihedralType, dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4); }

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

BOUNDS calculateDataBoundary (DATA_ATOMS *atoms, DATAFILE_INFO datafileInfo, BOUNDS originalBoundary)
{
	for (int i = 0; i < datafileInfo.nAtoms; ++i)
	{
		if (i == 0)
		{
			originalBoundary.xlo = atoms[i].x; originalBoundary.xhi = atoms[i].x;
			originalBoundary.ylo = atoms[i].y; originalBoundary.yhi = atoms[i].y;
			originalBoundary.zlo = atoms[i].z; originalBoundary.zhi = atoms[i].z;
		}
		else
		{
			if (atoms[i].x < originalBoundary.xlo) {
				originalBoundary.xlo = atoms[i].x; }
			else if (atoms[i].x > originalBoundary.xhi) {
				originalBoundary.xhi = atoms[i].x; }

			if (atoms[i].y < originalBoundary.ylo) {
				originalBoundary.ylo = atoms[i].y; }
			else if (atoms[i].y > originalBoundary.yhi) 	{
				originalBoundary.yhi = atoms[i].y; }

			if (atoms[i].z < originalBoundary.zlo) {
				originalBoundary.zlo = atoms[i].z; }
			else if (atoms[i].z > originalBoundary.zhi) {
				originalBoundary.zhi = atoms[i].z; }
		}
	}
	return originalBoundary;
}

DUMP computeCenterData (DUMP center, DATA_ATOMS *atoms, int nAtoms)
{
	center.x = 0; center.y = 0; center.z = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		center.x += atoms[i].x;
		center.y += atoms[i].y;
		center.z += atoms[i].z;
	}

	center.x /= nAtoms;
	center.y /= nAtoms;
	center.z /= nAtoms;

	return center;
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

DATA_ATOMS *translateDataAtoms (DATA_ATOMS *atoms, DATAFILE_INFO datafileInfo, DUMP *goldCoords, int nGoldAtoms, BOUNDS goldBoundary, BOUNDS originalBoundary)
{
	// Calculating the center for datafile atoms
	DUMP centerDatafile, centerGold;
	centerDatafile = computeCenterData (centerDatafile, atoms, datafileInfo.nAtoms);
	centerGold = computeCenterDump (centerGold, goldCoords, nGoldAtoms);

	printf("Center of mass for gold:\n\nx: %f; y: %f; z: %f\n\nCenter of mass for datafile atoms:\n\nx: %f; y: %f; z: %f\n\n", centerGold.x, centerGold.y, centerGold.z, centerDatafile.x, centerDatafile.y, centerDatafile.z);

	float thresholdDistance = 2.0;
	float translate_x = (centerGold.x - centerDatafile.x), translate_y = (centerGold.y - centerDatafile.y), translate_z = (goldBoundary.zhi - originalBoundary.zlo) + thresholdDistance;

	for (int i = 0; i < datafileInfo.nAtoms; ++i)
	{
		atoms[i].x += translate_x;
		atoms[i].y += translate_y;
		atoms[i].z += translate_z;
	}

	return atoms;
}

BOUNDS calculateOverallBoundary (BOUNDS newBoundary, BOUNDS goldBoundary, BOUNDS overallBoundary)
{
	if (newBoundary.xhi > goldBoundary.xhi)
		overallBoundary.xhi = newBoundary.xhi;
	else
		overallBoundary.xhi = goldBoundary.xhi;

	if (newBoundary.yhi > goldBoundary.yhi)
		overallBoundary.yhi = newBoundary.yhi;
	else
		overallBoundary.yhi = goldBoundary.yhi;

	if (newBoundary.zhi > goldBoundary.zhi)
		overallBoundary.zhi = newBoundary.zhi;
	else
		overallBoundary.zhi = goldBoundary.zhi;

	if (newBoundary.xlo < goldBoundary.xlo)
		overallBoundary.xlo = newBoundary.xlo;
	else
		overallBoundary.xlo = goldBoundary.xlo;

	if (newBoundary.ylo < goldBoundary.ylo)
		overallBoundary.ylo = newBoundary.ylo;
	else
		overallBoundary.ylo = goldBoundary.ylo;

	if (newBoundary.zlo < goldBoundary.zlo)
		overallBoundary.zlo = newBoundary.zlo;
	else
		overallBoundary.zlo = goldBoundary.zlo;

	return overallBoundary;
}

int main(int argc, char const *argv[])
{
	FILE *inputData, *inputGold, *outputData, *outputXYZ;
	inputData = fopen (argv[1], "r");
	inputGold = fopen (argv[2], "r");
	outputData = fopen (argv[3], "w");
	outputXYZ = fopen (argv[4], "w");

	int nGoldAtoms = getNatoms (inputGold);

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafileInfo;

	BOUNDS datafileBoundary, goldBoundary, originalBoundary, newBoundary, overallBoundary;
	ATOMIC_MASS *mass;

	// Reading the input data file
	readData (inputData, &atoms, &bonds, &angles, &dihedrals, &impropers, &datafileInfo, &datafileBoundary, &mass);

	DUMP *goldCoords;
	goldCoords = (DUMP *) malloc (nGoldAtoms * sizeof (DUMP));

	// Reading the last timeframe of gold dump file
	char *pipe_lastframe;
	pipe_lastframe = (char *) malloc (50 * sizeof (char));
	snprintf (pipe_lastframe, 50, "tail -%d %s", nGoldAtoms, argv[2]);
	goldCoords = readLastDumpFrame (pipe_lastframe, nGoldAtoms);

	goldBoundary = calculateDumpBoundary (goldCoords, nGoldAtoms, goldBoundary);
	originalBoundary = calculateDataBoundary (atoms, datafileInfo, originalBoundary);

	atoms = translateDataAtoms (atoms, datafileInfo, goldCoords, nGoldAtoms, goldBoundary, originalBoundary);
	newBoundary = calculateDataBoundary (atoms, datafileInfo, newBoundary);
	overallBoundary = calculateOverallBoundary (newBoundary, goldBoundary, overallBoundary);

	printf("Overall bounds after adding polarizable gold:\n\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\n", overallBoundary.xlo, overallBoundary.xhi, overallBoundary.ylo, overallBoundary.yhi, overallBoundary.zlo, overallBoundary.zhi);

	printDatafile (atoms, bonds, angles, dihedrals, impropers, datafileInfo, datafileBoundary, mass, goldCoords, nGoldAtoms, overallBoundary, outputData, outputXYZ);

	// Create electron cloud and bond information for the gold surface

	// Calculate the combined bounds for gold and surfactant molecules

	// Print the new data file with gold and electron positions

	return 0;
}