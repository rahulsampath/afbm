
#include "fullDomainUtils.h"

void createNeumannMatrix_Full(DA da, Mat mat) {
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




