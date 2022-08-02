#ifndef READINPUTS_MICELLE
#define READINPUTS_MICELLE

BONDS **readBonds (BONDS **inputBonds, int nSurfactants, SURFACTANT *inputStructures);
COORDINATES **readCoordinates (COORDINATES **inputCoordinates, int nSurfactants, SURFACTANT *inputStructures);
SURFACTANT *getNBonds (SURFACTANT *inputStructures, int nSurfactants);
SURFACTANT *getNAtoms (SURFACTANT *inputStructures, int nSurfactants);
SURFACTANT *storeSurfactantInformation (SURFACTANT *inputStructures, int nSurfactants, FILE *readConfig);
void checkArguments (int argc);
int checkNSurfactants (FILE *readConfig);

#endif