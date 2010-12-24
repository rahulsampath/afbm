
#include "fullDomainUtils.h"
#include <iostream>
#include <cstdlib>

#include "global.h"

#define __RG_NODE_ID__(xi, yi, zi, N) ( ( ( ((zi)*(N)) + (yi) )*(N) ) + (xi) )

extern double stencil[64];

void createNeumannMatrix_Full(DA da, Mat mat) {

  PetscInt N;
  DAGetInfo(da, PETSC_NULL, &N, PETSC_NULL, PETSC_NULL, 
      PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
  double h = 1.0/(static_cast<double>(N) - 1.0);

  MatZeroEntries(mat);

  for(int zi = 0; zi < (N - 1); zi++) {
    for(int yi = 0; yi < (N - 1); yi++) {
      for(int xi = 0; xi < (N - 1); xi++) {
        PetscInt indices[8];
        indices[0] = __RG_NODE_ID__(xi, yi, zi, N);
        indices[1] = __RG_NODE_ID__((xi + 1), yi, zi, N);
        indices[2] = __RG_NODE_ID__(xi, (yi + 1), zi, N);
        indices[3] = __RG_NODE_ID__((xi + 1), (yi + 1), zi, N);
        indices[4] = __RG_NODE_ID__(xi, yi, (zi + 1), N);
        indices[5] = __RG_NODE_ID__((xi + 1), yi, (zi + 1), N);
        indices[6] = __RG_NODE_ID__(xi, (yi + 1), (zi + 1), N);
        indices[7] = __RG_NODE_ID__((xi + 1), (yi + 1), (zi + 1), N);
        MatSetValues(mat, 8, indices, 8, indices, stencil, ADD_VALUES);
      }//end for xi
    }//end for yi
  }//end for zi

  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

  MatScale(mat, (h/2.0));

}

void applyDirichletMatrixCorrection_Full(DA da, Mat mat) {
}

void setDirichletValues_Full(DA da, Vec vec) {
}

void computeRHS_Full(DA da, Vec rhs) {
}

double phi(int node, double psi, double eta, double gamma) {

  double val = 0;

  if(node == 0) { 
    val = (1.0 - psi)*(1.0 - eta)*(1.0 - gamma)/8.0;
  } else if(node == 1) {
    val = (1.0 + psi)*(1.0 - eta)*(1.0 - gamma)/8.0;
  } else if(node == 2) {
    val = (1.0 - psi)*(1.0 + eta)*(1.0 - gamma)/8.0;
  } else if(node == 3) {
    val = (1.0 + psi)*(1.0 + eta)*(1.0 - gamma)/8.0;
  } else if(node == 4) {
    val = (1.0 - psi)*(1.0 - eta)*(1.0 + gamma)/8.0;
  } else if(node == 5) {
    val = (1.0 + psi)*(1.0 - eta)*(1.0 + gamma)/8.0;
  } else if(node == 6) {
    val = (1.0 - psi)*(1.0 + eta)*(1.0 + gamma)/8.0;
  } else if(node == 7) {
    val = (1.0 + psi)*(1.0 + eta)*(1.0 + gamma)/8.0;
  } else {
    std::cout<<"Error!!!"<<std::endl;
    exit(1);
  }

  return val;

}




