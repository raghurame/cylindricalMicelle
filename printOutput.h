#ifndef PACKING_MICELLE
#define PACKING_MICELLE

void writeCar (COORDINATES *outputCoordinates, int totalAtoms, SURFACTANT *inputStructures, int nSurfactants);
void writeMdf (BONDS *outputBonds, int totalAtoms, SURFACTANT *inputStructures, int nSurfactants);

#endif