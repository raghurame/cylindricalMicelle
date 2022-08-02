#ifndef STRUCT_MICELLES
#define STRUCT_MICELLES

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

typedef struct farthestPoints
{
	float x1, y1, z1, x2, y2, z2, distance;
	int sino1, sino2;
} FARTHESTPOINTS;

typedef struct com
{
	float x, y, z;
} CENTEROFMASS;

typedef struct cartesian
{
	float x, y, z;
} CARTESIAN;

#endif