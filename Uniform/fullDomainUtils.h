
#ifndef __FULL_DOMAIN__ 
#define __FULL_DOMAIN__

#include "petscda.h"

void createNeumannMatrix_Full(DA da, Mat mat);

void applyDirichletMatrixCorrection_Full(DA da, Mat mat);

void setDirichletValues_Full(DA da, Vec vec);

void computeRHS_Full(DA da, Vec rhs);

double phi(int node, double psi, double eta, double gamma);

#endif



