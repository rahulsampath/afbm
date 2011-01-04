
#ifndef __FAT_BOUNDARY__
#define __FAT_BOUNDARY__

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "linear_implicit_system.h"

void createNeumannMatrix_Fat(Mat & mat);

void dirichletMatCorrection_Fat(LinearImplicitSystem& system , Mat & stiffnesMatrix, MeshBase & mesh);

void dirichletVecCorrection_Fat();

void computeRHS_Fat(Mat  & stiffnesMatrix, Vec & dirichletVec);

void getBoundary_Fat(Vec & fatBnd);

void getDiracFunctions_Fat(LinearImplicitSystem& system, std::vector<double> solVec, std::vector<double> rhs, MeshBase & mesh);

#endif


