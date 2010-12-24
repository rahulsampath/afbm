
#ifndef __FAT_BOUNDARY__
#define __FAT_BOUNDARY__

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "linear_implicit_system.h"

void createNeumannMatrix_Fat(Mat & mat);

void dirichletMatCorrection_Fat(LinearImplicitSystem& system , Mat & stiffnesMatrix, MeshBase & mesh);

void dirichletVecCorrection_Fat();

void computeRHS_Fat();

#endif


