
#ifndef __FAT_BOUNDARY__
#define __FAT_BOUNDARY__

#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "linear_implicit_system.h"

void createNeumannMatrix_Fat(Mat & petscMat, Vec & petscVec,
    LinearImplicitSystem & system, const DofMap & dof_map, MeshBase & mesh);

void dirichletMatCorrection_Fat(Mat mat, const DofMap & dof_map, MeshBase & mesh, const int bnd_id);

void dirichletVecSetCorrection1_Fat(Vec dirichletVec, const DofMap & dof_map, MeshBase & mesh);

void dirichletVecSetCorrection2_Fat( Vec dirichletVec, Vec FullDomainSolVec, int N,
    const DofMap & dof_map, MeshBase & mesh);

void dirichletVecAddCorrection2_Fat(Mat & stiffnessMatrix, Vec correctionVec, Vec FullDomainSolutionVec, int N,
    const DofMap & dof_map, MeshBase & mesh);

void getDiracFunctions_Fat(LinearImplicitSystem& system, Vec & solVec, Vec & rhs, MeshBase & mesh);

void getBoundary_Fat(LinearImplicitSystem& system, std::vector<double> fatBnd, MeshBase & mesh);

#endif


