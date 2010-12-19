
#ifndef __FAT_BOUNDARY__
#define __FAT_BOUNDARY__

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"

void createNeumannMatrix_Fat(Mat mat);

void dirichletMatCorrection_Fat();

void dirichletVecCorrection_Fat();

void computeRHS_Fat();

#endif


