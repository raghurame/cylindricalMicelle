#ifndef PACKING_MICELLE
#define PACKING_MICELLE

float calculateDistance (float x1, float y1, float z1, float x2, float y2, float z2);
FARTHESTPOINTS *calculateFarPoints (FARTHESTPOINTS *inputStructures_farPoints, COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures);
CENTEROFMASS *computeCOM (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures);
COORDINATES **translateCoordinates (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, CENTEROFMASS *comForMolecules, FARTHESTPOINTS *inputStructures_farPoints);
COORDINATES **alignMolecule (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, FARTHESTPOINTS *inputStructures_farPoints);
COORDINATES **orientSurfactants (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures, FARTHESTPOINTS *inputStructures_farPoints);
void computeLongestDimension (CARTESIAN **loDimension, CARTESIAN **hiDimension, COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures);
int countTotalAtoms (SURFACTANT *inputStructures, int nSurfactants);
int countTotalMolecules (SURFACTANT *inputStructures, int nSurfactants);
MOLECULELOG *initializeAllMoleculesLog (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants);
int findRemainingMolecules (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants);
MOLECULELOG *pickMolecule (MOLECULELOG *allMoleculeLog, SURFACTANT *inputStructures, int nSurfactants, double randomFlip, int *pickedSurfactantMolID);
CARTESIAN computeMaxSurfactantLength (CARTESIAN maxSurfactantLenght, CARTESIAN globalSurfactantlo, CARTESIAN globalSurfactanthi, float tolerance);
COORDINATES *replicateSurfactants (COORDINATES **inputCoordinates, BONDS **inputBonds, CARTESIAN globalSurfactantlo, CARTESIAN globalSurfactanthi, int nSurfactants, SURFACTANT *inputStructures);
void calculateGlobalMinMax (CARTESIAN *globalSurfactanthi, CARTESIAN *globalSurfactantlo, CARTESIAN *loDimension, CARTESIAN *hiDimension, int nSurfactants);
BONDS *addBonds (COORDINATES *outputCoordinates, COORDINATES **inputCoordinates, BONDS **inputBonds, SURFACTANT *inputStructures, int nSurfactants);


#endif