#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

/*
Read the data file (argv[1])
Define the range of atoms in the data file to modify (argv[2], argv[3])
Check the number of atoms
Read the pdb file (argv[4])
Define the range of atoms in the pdb file to substitute (argv[5], argv[6])
Check the number of atoms

If both the numbers match, replace the coordinates
Print the new datafile
*/

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct pdb_atoms
{
	int sino, molType;
	float x, y, z, float1, charge;
	char *atomName, *molName, *char1;
} PDB_ATOMS;

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
	rewind (input);

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

	// printf("\nPrinting boundary information from data file:\n\n  xlo: %f; xhi: %f\n  ylo: %f; yhi: %f\n  zlo: %f; zhi: %f\n", (*datafileBoundary).xlo, (*datafileBoundary).xhi, (*datafileBoundary).ylo, (*datafileBoundary).yhi, (*datafileBoundary).zlo, (*datafileBoundary).zhi);
	// printf("\nFrom input data file:\n\n  nAtoms: %d\n  nBonds: %d\n  nAngles: %d\n  nDihedrals: %d\n  nImpropers: %d\n\n", (*datafile).nAtoms, (*datafile).nBonds, (*datafile).nAngles, (*datafile).nDihedrals, (*datafile).nImpropers);
	// printf("\nMasses:\n\n");
	// for (int i = 0; i < (*datafile).nAtomTypes; ++i) {
	// 	printf("%d %f\n", (*mass)[i].atomType, (*mass)[i].mass); }
	// printf("\n");
	// printf("Max mol type present in the data file: %d\n\n", (*datafile).maxMolType);
	rewind (input);
}

PDB_ATOMS *readPDB (FILE *inputPDB, PDB_ATOMS *pdbCoordinates, int nPDBAtomsToAdd)
{
	char lineString[2000];
	int lineNumber = 0;
	rewind (inputPDB);

	while (fgets (lineString, 2000, inputPDB) != NULL)
	{
		sscanf (lineString, "%d\n", &lineNumber);

		if (lineNumber >= 0 && lineNumber < nPDBAtomsToAdd)
		{
			sscanf (lineString, "%d %s %s %s %d %f %f %f %f %f\n", 
				&pdbCoordinates[lineNumber - 1].sino, 
				pdbCoordinates[lineNumber - 1].atomName, 
				pdbCoordinates[lineNumber - 1].molName, 
				pdbCoordinates[lineNumber - 1].char1, 
				&pdbCoordinates[lineNumber - 1].molType, 
				&pdbCoordinates[lineNumber - 1].x, 
				&pdbCoordinates[lineNumber - 1].y, 
				&pdbCoordinates[lineNumber - 1].z, 
				&pdbCoordinates[lineNumber - 1].float1, 
				&pdbCoordinates[lineNumber - 1].charge);
		}
	}

	rewind (inputPDB);
	return pdbCoordinates;
}

PDB_ATOMS *initializePDBchars (PDB_ATOMS *pdbCoordinates, int nPDBAtomsToAdd)
{
	for (int i = 0; i < nPDBAtomsToAdd; ++i)
	{
		pdbCoordinates[i].atomName = (char *) malloc (10 * sizeof (char));
		pdbCoordinates[i].molName = (char *) malloc (10 * sizeof (char));
		pdbCoordinates[i].char1 = (char *) malloc (10 * sizeof (char));
	}

	return pdbCoordinates;
}

typedef struct atomIndices
{
	int start, end;
} ATOM_INDICES;

void printNewData (DATA_ATOMS *atoms, DATA_BONDS *bonds, DATA_ANGLES *angles, DATA_DIHEDRALS *dihedrals, DATA_IMPROPERS *impropers, DATAFILE_INFO datafileInfo, BOUNDS overallBoundary, ATOMIC_MASS *mass, PDB_ATOMS *pdbCoordinates, ATOM_INDICES pdb, ATOM_INDICES datafile, FILE *outputData, FILE *outputXYZ)
{
	int nDataAtomsToReplace = (datafile.end - datafile.start + 1), nPDBAtomsToAdd = (pdb.end - pdb.start + 1);
	// Printing the header information
	fprintf(outputData, "%s\n\n", "LAMMPS data file containing the original structures and the polarizable gold");
	fprintf(outputData, "%d atoms\n%d bonds\n%d angles\n%d dihedrals\n%d impropers\n\n%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nMasses\n\n", 
		datafileInfo.nAtoms, 
		datafileInfo.nBonds, 
		datafileInfo.nAngles, 
		datafileInfo.nDihedrals, 
		datafileInfo.nImpropers, 
		datafileInfo.nAtomTypes, 
		datafileInfo.nBondTypes, 
		datafileInfo.nAngleTypes, 
		datafileInfo.nDihedralTypes, 
		overallBoundary.xlo, 
		overallBoundary.xhi, 
		overallBoundary.ylo, 
		overallBoundary.yhi, 
		overallBoundary.zlo, 
		overallBoundary.zhi);

	fprintf(outputXYZ, "%d\n", datafileInfo.nAtoms);
	fprintf(outputXYZ, "%s\n", "Dummy comment line");

	// Printing masses
	for (int i = 0; i < datafileInfo.nAtomTypes; ++i) {
		fprintf(outputData, "%d %f\n", mass[i].atomType, mass[i].mass); }

	// Printing atoms from data file
	fprintf(outputData, "\nAtoms\n\n");

	int pdbIndex = 0;

	for (int i = 0; i < datafileInfo.nAtoms; ++i) 
	{
		fprintf(outputXYZ, "C %f %f %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
		fprintf(outputData, "%d %d %d %f %f %f %f\n", atoms[i].id, atoms[i].molType, atoms[i].atomType, atoms[i].charge, atoms[i].x, atoms[i].y, atoms[i].z); 
	}

	// Printing bonds
	fprintf(outputData, "\nBonds\n\n");
	for (int i = 0; i < datafileInfo.nBonds; ++i) {
		fprintf(outputData, "%d %d %d %d\n", bonds[i].id, bonds[i].bondType, bonds[i].atom1, bonds[i].atom2); }

	// Printing angles
	fprintf(outputData, "\nAngles\n\n");
	for (int i = 0; i < datafileInfo.nAngles; ++i) {
		fprintf(outputData, "%d %d %d %d %d\n", angles[i].id, angles[i].angleType, angles[i].atom1, angles[i].atom2, angles[i].atom3); }

	// Printing dihedrals
	fprintf(outputData, "\nDihedrals\n\n");
	for (int i = 0; i < datafileInfo.nDihedrals; ++i) {
		fprintf(outputData, "%d %d %d %d %d %d\n", dihedrals[i].id, dihedrals[i].dihedralType, dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4); }

}

BOUNDS findPDBboundary (PDB_ATOMS *pdbCoordinates, int nPDBAtomsToAdd, BOUNDS pdbBoundary)
{
	pdbBoundary.xlo = 0;
	pdbBoundary.xhi = 0;
	pdbBoundary.ylo = 0;
	pdbBoundary.yhi = 0;

	for (int i = 0; i < nPDBAtomsToAdd; ++i)
	{
		if (i == 0)
		{
			pdbBoundary.xlo = pdbCoordinates[i].x; pdbBoundary.xhi = pdbCoordinates[i].x;
			pdbBoundary.ylo = pdbCoordinates[i].y; pdbBoundary.yhi = pdbCoordinates[i].y;
			pdbBoundary.zlo = pdbCoordinates[i].z; pdbBoundary.zhi = pdbCoordinates[i].z;
		}
		else
		{
			if (pdbCoordinates[i].x < pdbBoundary.xlo) {
				pdbBoundary.xlo = pdbCoordinates[i].x; }
			else if (pdbCoordinates[i].x > pdbBoundary.xhi) {
				pdbBoundary.xhi = pdbCoordinates[i].x; }

			if (pdbCoordinates[i].y < pdbBoundary.ylo) {
				pdbBoundary.ylo = pdbCoordinates[i].y; }
			else if (pdbCoordinates[i].y > pdbBoundary.yhi) {
				pdbBoundary.yhi = pdbCoordinates[i].y; }

			if (pdbCoordinates[i].z < pdbBoundary.zlo) {
				pdbBoundary.zlo = pdbCoordinates[i].z; }
			else if (pdbCoordinates[i].z > pdbBoundary.zhi) {
				pdbBoundary.zhi = pdbCoordinates[i].z; }
		}
	}

	return pdbBoundary;
}

BOUNDS getOverallBoundary (BOUNDS overallBoundary, BOUNDS datafileBoundary, BOUNDS pdbBoundary)
{
	// printf("\nBoundary of datafile:\n\nxlo: %f; xhi: %f;\nylo: %f; yhi: %f;\nzlo: %f; zhi: %f;\n", 
	// 	datafileBoundary.xlo, 
	// 	datafileBoundary.xhi, 
	// 	datafileBoundary.ylo, 
	// 	datafileBoundary.yhi, 
	// 	datafileBoundary.zlo, 
	// 	datafileBoundary.zhi);

	// printf("\nBoundary of PDB file:\n\nxlo: %f; xhi: %f;\nylo: %f; yhi: %f;\nzlo: %f; zhi: %f;\n", 
	// 	pdbBoundary.xlo, 
	// 	pdbBoundary.xhi, 
	// 	pdbBoundary.ylo, 
	// 	pdbBoundary.yhi, 
	// 	pdbBoundary.zlo, 
	// 	pdbBoundary.zhi);

	if (datafileBoundary.xlo < pdbBoundary.xlo) {
		overallBoundary.xlo = datafileBoundary.xlo; }
	else {
		overallBoundary.xlo = pdbBoundary.xlo; }

	if (datafileBoundary.ylo < pdbBoundary.ylo) {
		overallBoundary.ylo = datafileBoundary.ylo; }
	else {
		overallBoundary.ylo = pdbBoundary.ylo; }

	if (datafileBoundary.zlo < pdbBoundary.zlo) {
		overallBoundary.zlo = datafileBoundary.zlo; }
	else {
		overallBoundary.zlo = pdbBoundary.zlo; }

	if (datafileBoundary.xhi < pdbBoundary.xhi) {
		overallBoundary.xhi = datafileBoundary.xhi; }
	else {
		overallBoundary.xhi = pdbBoundary.xhi; }

	if (datafileBoundary.yhi < pdbBoundary.yhi) {
		overallBoundary.yhi = datafileBoundary.yhi; }
	else {
		overallBoundary.yhi = pdbBoundary.yhi; }

	if (datafileBoundary.zhi < pdbBoundary.zhi) {
		overallBoundary.zhi = datafileBoundary.zhi; }
	else {
		overallBoundary.zhi = pdbBoundary.zhi; }

	// printf("\nOverall boundary:\n\nxlo: %f; xhi: %f;\nylo: %f; yhi: %f;\nzlo: %f; zhi: %f;\n", 
	// 	overallBoundary.xlo, 
	// 	overallBoundary.xhi, 
	// 	overallBoundary.ylo, 
	// 	overallBoundary.yhi, 
	// 	overallBoundary.zlo, 
	// 	overallBoundary.zhi);

	return overallBoundary;
}

DATA_ATOMS *assignMolType (DATA_ATOMS *atoms, DATAFILE_INFO datafileInfo)
{
	int previousBromine = 0, currentBromine, molLength;
	int atomTypeBr = 5;

	// Reassign molType for gold and electron cloud
	// molType = 3 for gold core 
	// molType = 4 for electron cloud
	for (int i = 0; i < datafileInfo.nAtoms; ++i)
	{
		if (atoms[i].atomType == 6)
		{
			atoms[i].molType = 3;
		}
		else if (atoms[i].atomType == 7)
		{
			atoms[i].molType = 4;
		}
	}

	for (int i = 0; i < datafileInfo.nAtoms; ++i)
	{
		if (atoms[i].atomType == atomTypeBr)
		{
			currentBromine = i + 1;

			// Assigning molType
			for (int j = previousBromine; j < currentBromine; ++j)
			{
				if ((currentBromine - previousBromine) == 63)
				{
					atoms[j].molType = 1;
				}
				else if ((currentBromine - previousBromine) == 84)
				{
					atoms[j].molType = 2;
				}
			}

			previousBromine = currentBromine;
		}
	}

	return atoms;
}

DATA_ATOMS *replaceAtoms (DATA_ATOMS *atoms, DATAFILE_INFO datafileInfo, PDB_ATOMS *pdbCoordinates, int nPDBAtoms)
{
	int datafileIndex = 0, pdbIndex = 0;
	int nAtoms_CTAB = 63, nAtoms_DDAB = 84;

	// printf("nPDBAtoms: %d\n", nPDBAtoms);
	// printf("datafileInfo.nAtoms: %d\n", datafileInfo.nAtoms);
	// fflush (stdout);

	while (1)
	{
		if ((datafileIndex < datafileInfo.nAtoms))
		{
			if (pdbIndex >= nPDBAtoms) {
				pdbIndex = 0; }
			if (datafileIndex >= datafileInfo.nAtoms) {
				goto leaveThisWhileLoop; }

			if (pdbIndex < nPDBAtoms)
			{
				printf("datafileindex: %d ==> %d <--> %s (pdbIndex: %d)                 \r", datafileIndex, atoms[datafileIndex].molType, pdbCoordinates[pdbIndex].molName, pdbIndex);
				// fflush (stdout);

				if (atoms[datafileIndex].molType == 1 && strstr (pdbCoordinates[pdbIndex].molName, "CTA"))
				{
					printf("\nMatched with CTAB and molType = 1. Iterating for the next %d times\n", nAtoms_CTAB);
					fflush (stdout);
					for (int i = 0; i < nAtoms_CTAB; ++i)
					{
						atoms[datafileIndex].x = pdbCoordinates[pdbIndex].x;
						atoms[datafileIndex].y = pdbCoordinates[pdbIndex].y;
						atoms[datafileIndex].z = pdbCoordinates[pdbIndex].z;

						datafileIndex++; pdbIndex++;

						// printf("%d\n", datafileIndex);
						// fflush (stdout);

						if (pdbIndex >= nPDBAtoms) {
							pdbIndex = 0; }
						if (datafileIndex >= datafileInfo.nAtoms) {
							goto leaveThisWhileLoop; }
					}

				}
				else if ((atoms[datafileIndex].molType == 2) && strstr (pdbCoordinates[pdbIndex].molName, "DDA"))
				{
					printf("\nMatched with DDAB and molType = 2. Iterating for the next %d times\n", nAtoms_DDAB);
					fflush (stdout);

					for (int i = 0; i < nAtoms_DDAB; ++i)
					{
						atoms[datafileIndex].x = pdbCoordinates[pdbIndex].x;
						atoms[datafileIndex].y = pdbCoordinates[pdbIndex].y;
						atoms[datafileIndex].z = pdbCoordinates[pdbIndex].z;

						datafileIndex++; pdbIndex++;

						// printf("%d\n", datafileIndex);
						// fflush (stdout);

						if (pdbIndex >= nPDBAtoms) {
							pdbIndex = 0; }
						if (datafileIndex >= datafileInfo.nAtoms) {
							goto leaveThisWhileLoop; }
					}
				}
				else if (atoms[datafileIndex].molType > 2)
				{
					datafileIndex++;
				}

				if (pdbIndex < nPDBAtoms)
				{
					pdbIndex++;
				}
			}
		}
		else if (datafileIndex >= datafileInfo.nAtoms)
		{
			goto leaveThisWhileLoop;
		}
	}

	leaveThisWhileLoop: ;

	return atoms;
}

int findNPDBAtoms (FILE *inputPDB)
{
	rewind (inputPDB);

	int nLines = 0;
	char lineString[2000];

	while (fgets (lineString, 2000, inputPDB) != NULL)
	{
		nLines++;
	}

	rewind (inputPDB);
	return nLines;
}

int main(int argc, char const *argv[])
{
	if (argc == 1) {
		(void)printf("\nARGUMENTS:\n~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = input datafile\n{~} argv[2] = starting atom index to replace, in datafile\n{~} argv[3] = final atom index to replace, in datafile\n{~} argv[4] = input pdb filename\n{~} argv[5] = starting atom index to substitute, from pdb file\n{~} argv[6] = final atom index to substitute, from the pdb file.\n\nNOTE: The number of atoms selected in the datafile and the pdb file must match!\n\n");
		exit (1); }

	FILE *inputData, *inputPDB, *outputData, *outputXYZ;
	inputData = fopen (argv[1], "r");
	inputPDB = fopen (argv[4], "r");

	ATOM_INDICES pdb, datafile;
	datafile.start = atoi (argv[2]);
	datafile.end = atoi (argv[3]);
	pdb.start = atoi (argv[5]);
	pdb.end = atoi (argv[6]);

	int nDataAtomsToReplace = datafile.end - datafile.start + 1, nPDBAtomsToAdd = pdb.end - pdb.start + 1, nPDBAtomsAvailable = findNPDBAtoms (inputPDB);

	if (nDataAtomsToReplace != nPDBAtomsToAdd)
	{
		printf("ERROR: Number of atoms specified are different between PDB and Datafile...\n\n");
		exit (1);
	}

	if (nPDBAtomsToAdd < nPDBAtomsAvailable)
	{
		printf("ERROR: The number of selected atoms in PDB file is greater than the available atoms.\n");
		exit (1);
	}
	
	char *outputDataString, *outputXYZstring;
	outputDataString = (char *) malloc (100 * sizeof (char));
	outputXYZstring = (char *) malloc (100 * sizeof (char));
	snprintf (outputDataString, 100, "%s.data", argv[1]);
	snprintf (outputXYZstring, 100, "%s.xyz", argv[1]);
	outputData = fopen (outputDataString, "w");
	outputXYZ = fopen (outputXYZstring, "w");

	// Read data file
	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafileInfo;

	BOUNDS datafileBoundary, pdbBoundary /*, overallBoundary*/;
	ATOMIC_MASS *mass;

	// Reading the input data file
	readData (inputData, &atoms, &bonds, &angles, &dihedrals, &impropers, &datafileInfo, &datafileBoundary, &mass);

	PDB_ATOMS *pdbCoordinates;
	pdbCoordinates = (PDB_ATOMS *) malloc (nPDBAtomsToAdd * sizeof (PDB_ATOMS));

	// Allocate memory for the char types in PDB_ATOMS variable
	pdbCoordinates = initializePDBchars (pdbCoordinates, nPDBAtomsToAdd);

	// Store all the PDB information in PDB_ATOMS variable
	pdbCoordinates = readPDB (inputPDB, pdbCoordinates, nPDBAtomsToAdd);

	pdbBoundary = findPDBboundary (pdbCoordinates, nPDBAtomsToAdd, pdbBoundary);

	// overallBoundary = getOverallBoundary (overallBoundary, datafileBoundary, pdbBoundary);

	atoms = assignMolType (atoms, datafileInfo);

	// Then replace the coordinates if the molType matches.
	// If the molType does not match, then iterate forward and check again.
	// If the end of array is reached, then start from the beginning of the array (cycle through).

	atoms = replaceAtoms (atoms, datafileInfo, pdbCoordinates, nPDBAtomsAvailable);

	printNewData (atoms, bonds, angles, dihedrals, impropers, datafileInfo, datafileBoundary, mass, pdbCoordinates, pdb, datafile, outputData, outputXYZ);

	fclose (inputData);
	fclose (inputPDB);
	fclose (outputData);
	fclose (outputXYZ);
	return 0;
}