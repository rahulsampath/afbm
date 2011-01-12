
#ifndef __FULL_DOMAIN__ 
#define __FULL_DOMAIN__

#include "petscda.h"
#include "petscdmmg.h"
#include <vector>

PetscErrorCode createMatrix_Full(DMMG dmmg, Mat J, Mat B);

void createNeumannMatrix_Full(DA da, Mat mat);

void applyDirichletMatrixCorrection_Full(DA da, Mat mat);

void setDirichletValues_Full(DA da, Vec vec);

void computeRHSterm1_Full(DA da, Vec vec);

double phi_Full(int node, double psi, double eta, double gamma);

void phi_Full(double x, double y, double z, int N,
    std::vector<unsigned int> & indices, std::vector<double> & vals);

double evalFn(double x, double y, double z);

double interpolateAtPt(Vec vec, int N, double x, double y, double z);

#endif



