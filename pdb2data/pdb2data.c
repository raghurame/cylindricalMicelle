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

void printNewData (DATA_ATOMS *atoms, DATA_BONDS *bonds, DATA_ANGLES *angles, DATA_DIHEDRALS *dihedrals, DATA_IMPROPERS *impropers, DATAFILE_INFO datafileInfo, BOUNDS datafileBoundary, ATOMIC_MASS *mass, PDB_ATOMS *pdbCoordinates, ATOM_INDICES pdb, ATOM_INDICES datafile)
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
		datafileBoundary.xlo, 
		datafileBoundary.xhi, 
		datafileBoundary.ylo, 
		datafileBoundary.yhi, 
		datafileBoundary.zlo, 
		datafileBoundary.zhi);

}

int main(int argc, char const *argv[])
{
	if (argc == 1) {
		(void)printf("\nARGUMENTS:\n~~~~~~~~~~\n\n{~} argv[0] = program\n{~} argv[1] = input datafile\n{~} argv[2] = starting atom index to replace, in datafile\n{~} argv[3] = final atom index to replace, in datafile\n{~} argv[4] = input pdb filename\n{~} argv[5] = starting atom index to substitute, from pdb file\n{~} argv[6] = final atom index to substitute, from the pdb file.\n\nNOTE: The number of atoms selected in the datafile and the pdb file must match!\n\n");
		exit (1); }

	FILE *inputData, *inputPDB, *outputData;
	inputData = fopen (argv[1], "r");
	inputPDB = fopen (argv[4], "r");

	ATOM_INDICES pdb, datafile;
	datafile.start = atoi (argv[2]);
	datafile.end = atoi (argv[3]);
	pdb.start = atoi (argv[5]);
	pdb.end = atoi (argv[6]);

	int nDataAtomsToReplace = datafile.end - datafile.start + 1, nPDBAtomsToAdd = pdb.end - pdb.start + 1;
	
	char *outputDataString;
	outputDataString = (char *) malloc (100 * sizeof (char));
	snprintf (outputDataString, 100, "%s.mod", argv[1]);
	outputData = fopen (outputDataString, "w");

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

	PDB_ATOMS *pdbCoordinates;
	pdbCoordinates = (PDB_ATOMS *) malloc (nPDBAtomsToAdd * sizeof (PDB_ATOMS));

	// Allocate memory for the char types in PDB_ATOMS variable
	pdbCoordinates = initializePDBchars (pdbCoordinates, nPDBAtomsToAdd);

	// Store all the PDB information in PDB_ATOMS variable
	pdbCoordinates = readPDB (inputPDB, pdbCoordinates, nPDBAtomsToAdd);

	printNewData (atoms, bonds, angles, dihedrals, impropers, datafileInfo, datafileBoundary, mass, pdbCoordinates, pdb, datafile);

	fclose (inputData);
	fclose (inputPDB);
	fclose (outputData);
	return 0;
}