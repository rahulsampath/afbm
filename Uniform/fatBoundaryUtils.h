
#ifndef __FAT_BOUNDARY__
#define __FAT_BOUNDARY__

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"

void createNeumann(unsigned int Dim);

void dirichletMatCorrection();

void dirichletVecCorrection();

void computeRHS();

#endif


