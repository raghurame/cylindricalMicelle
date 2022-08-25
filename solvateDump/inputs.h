#ifndef STRUCT_SOLVATEDUMP
#define STRUCT_SOLVATEDUMP

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers, DATAFILE_INFO *datafile, BOUNDS *datafileBoundary);
DUMP *readLastDumpFrame (const char *pipeString, int nAtoms);
int getNatoms (const char *inputFileName);

#endif