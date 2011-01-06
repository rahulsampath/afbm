
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

void computeRHSterm1_Full(DA da, Vec vec) {

  PetscInt N;
  DAGetInfo(da, PETSC_NULL, &N, PETSC_NULL, PETSC_NULL, 
      PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
  double h = 1.0/(static_cast<double>(N) - 1.0);

  double gPt[] = { 0,  (sqrt(15.0)/5.0), (-sqrt(15.0)/5.0) };
  double gWt[] = { (8.0/9.0), (5.0/9.0), (5.0/9.0) };

  VecZeroEntries(vec);

  for(int zi = 0; zi < (N - 1); zi++) {
    double z0 = h*static_cast<double>(zi);
    for(int yi = 0; yi < (N - 1); yi++) {
      double y0 = h*static_cast<double>(yi);
      for(int xi = 0; xi < (N - 1); xi++) {
        double x0 = h*static_cast<double>(xi);
        PetscInt indices[8];
        indices[0] = __RG_NODE_ID__(xi, yi, zi, N);
        indices[1] = __RG_NODE_ID__((xi + 1), yi, zi, N);
        indices[2] = __RG_NODE_ID__(xi, (yi + 1), zi, N);
        indices[3] = __RG_NODE_ID__((xi + 1), (yi + 1), zi, N);
        indices[4] = __RG_NODE_ID__(xi, yi, (zi + 1), N);
        indices[5] = __RG_NODE_ID__((xi + 1), yi, (zi + 1), N);
        indices[6] = __RG_NODE_ID__(xi, (yi + 1), (zi + 1), N);
        indices[7] = __RG_NODE_ID__((xi + 1), (yi + 1), (zi + 1), N);
        for(int node = 0; node < 8; node++) {
          double val = 0.0;
          for(int p = 0; p < 3; p++) {
            double z = z0 + (0.5*h*(1.0 + gPt[p]));
            for(int n = 0; n < 3; n++) {
              double y = y0 + (0.5*h*(1.0 + gPt[n]));
              for(int m = 0; m < 3; m++) {
                double x = x0 + (0.5*h*(1.0 + gPt[m]));
                val += (gWt[m]*gWt[n]*gWt[p]*evalFn(x, y, z)*phi_Full(node, gPt[m], gPt[n], gPt[p]));
              }//end for m
            }//end for n
          }//end for p
          VecSetValue(vec, indices[node], val, ADD_VALUES);
        }//end for node
      }//end for xi
    }//end for yi
  }//end for zi

  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);

  VecScale(vec, (h*h*h/8.0));

}

void applyDirichletMatrixCorrection_Full(DA da, Mat mat) {

  PetscInt N;
  DAGetInfo(da, PETSC_NULL, &N, PETSC_NULL, PETSC_NULL, 
      PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  {
    //Left
    int xi = 0;
    for(int zi = 0; zi < N; zi++) {
      for(int yi = 0; yi < N; yi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for yi
    }//end for zi
  }

  {
    //Right
    int xi = N - 1;
    for(int zi = 0; zi < N; zi++) {
      for(int yi = 0; yi < N; yi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for yi
    }//end for zi
  }

  {
    //Front
    int yi = 0;
    for(int zi = 0; zi < N; zi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for xi
    }//end for zi
  }

  {
    //Back
    int yi = N - 1;
    for(int zi = 0; zi < N; zi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for xi
    }//end for zi
  }

  {
    //Bottom
    int zi = 0;
    for(int yi = 0; yi < N; yi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for xi
    }//end for yi
  }

  {
    //Top
    int zi = N - 1;
    for(int yi = 0; yi < N; yi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        MatSetValue(mat, curr, curr, 1.0, INSERT_VALUES);
        for(int zj = (zi - 1); zj <= (zi + 1); zj++) {
          if((zj >= 0) && (zj < N)) {
            for(int yj = (yi - 1); yj <= (yi + 1); yj++) {
              if((yj >= 0) && (yj < N)) {
                for(int xj = (xi - 1); xj <= (xi + 1); xj++) {
                  if((xj >= 0) && (xj < N)) {
                    int other = __RG_NODE_ID__(xj, yj, zj, N);
                    if(curr != other) {
                      MatSetValue(mat, curr, other, 0.0, INSERT_VALUES);
                      MatSetValue(mat, other, curr, 0.0, INSERT_VALUES);
                    }
                  }
                }//end for xj
              }
            }//end for yj
          }
        }//end for zj
      }//end for xi
    }//end for yi
  }

  MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

}

void setDirichletValues_Full(DA da, Vec vec) {

  PetscInt N;
  DAGetInfo(da, PETSC_NULL, &N, PETSC_NULL, PETSC_NULL, 
      PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  {
    //Left
    int xi = 0;
    for(int zi = 0; zi < N; zi++) {
      for(int yi = 0; yi < N; yi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for yi
    }//end for zi
  }

  {
    //Right
    int xi = N - 1; 
    for(int zi = 0; zi < N; zi++) {
      for(int yi = 0; yi < N; yi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for yi
    }//end for zi
  }

  {
    //Front
    int yi = 0;
    for(int zi = 0; zi < N; zi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for xi
    }//end for zi
  }

  {
    //Back
    int yi = N - 1; 
    for(int zi = 0; zi < N; zi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for xi
    }//end for zi
  }

  {
    //Bottom
    int zi = 0;
    for(int yi = 0; yi < N; yi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for xi
    }//end for yi
  }

  {
    //Top
    int zi = N - 1; 
    for(int yi = 0; yi < N; yi++) {
      for(int xi = 0; xi < N; xi++) {
        int curr = __RG_NODE_ID__(xi, yi, zi, N);
        VecSetValue(vec, curr, 0.0, INSERT_VALUES);
      }//end for xi
    }//end for yi
  }

  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);

}

double phi_Full(int node, double psi, double eta, double gamma) {

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

void phi_Full(double x, double y, double z, int N,
    std::vector<unsigned int> & indices, std::vector<double> & vals) {

  double h = 1.0/static_cast<double>(N - 1);

  int xi = static_cast<int>(floor(x/h));
  int yi = static_cast<int>(floor(y/h));
  int zi = static_cast<int>(floor(z/h));

  double x0 = h*static_cast<double>(xi);
  double y0 = h*static_cast<double>(yi);
  double z0 = h*static_cast<double>(zi);

  double psi = (2.0*(x - x0)/h) - 1.0;
  double eta = (2.0*(y - y0)/h) - 1.0;
  double gamma = (2.0*(z - z0)/h) - 1.0;

  indices.resize(8);

  indices[0] = __RG_NODE_ID__(xi, yi, zi, N);
  indices[1] = __RG_NODE_ID__((xi + 1), yi, zi, N);
  indices[2] = __RG_NODE_ID__(xi, (yi + 1), zi, N);
  indices[3] = __RG_NODE_ID__((xi + 1), (yi + 1), zi, N);
  indices[4] = __RG_NODE_ID__(xi, yi, (zi + 1), N);
  indices[5] = __RG_NODE_ID__((xi + 1), yi, (zi + 1), N);
  indices[6] = __RG_NODE_ID__(xi, (yi + 1), (zi + 1), N);
  indices[7] = __RG_NODE_ID__((xi + 1), (yi + 1), (zi + 1), N);

  vals.resize(8);

  for(int node = 0; node < 8; node++) {
    vals[node] = phi_Full(node, psi, eta, gamma);
  }
}

double evalFn(double x, double y, double z) {

  double val = 0.0;

  if ( ( (__SQUARE__(x - __CENTER_X__)) +
        (__SQUARE__(y - __CENTER_Y__)) +
        (__SQUARE__(z - __CENTER_Z__)) ) > (__SQUARE__(__RADIUS__)) ) {
    val = __FN__(x, y, z);
  }

  return val;
}

double interpolateAtPt(Vec vec, int N, double x, double y, double z) {

  double h = 1.0/static_cast<double>(N - 1);

  int xi = static_cast<int>(floor(x/h));
  int yi = static_cast<int>(floor(y/h));
  int zi = static_cast<int>(floor(z/h));

  double x0 = h*static_cast<double>(xi);
  double y0 = h*static_cast<double>(yi);
  double z0 = h*static_cast<double>(zi);

  double psi = (2.0*(x - x0)/h) - 1.0;
  double eta = (2.0*(y - y0)/h) - 1.0;
  double gamma = (2.0*(z - z0)/h) - 1.0;

  unsigned int indices[8];

  indices[0] = __RG_NODE_ID__(xi, yi, zi, N);
  indices[1] = __RG_NODE_ID__((xi + 1), yi, zi, N);
  indices[2] = __RG_NODE_ID__(xi, (yi + 1), zi, N);
  indices[3] = __RG_NODE_ID__((xi + 1), (yi + 1), zi, N);
  indices[4] = __RG_NODE_ID__(xi, yi, (zi + 1), N);
  indices[5] = __RG_NODE_ID__((xi + 1), yi, (zi + 1), N);
  indices[6] = __RG_NODE_ID__(xi, (yi + 1), (zi + 1), N);
  indices[7] = __RG_NODE_ID__((xi + 1), (yi + 1), (zi + 1), N);

  PetscScalar *arr;

  VecGetArray(vec, &arr);

  double val = 0.0;

  for(int node = 0; node < 8; node++) {
    val += (arr[indices[node]]*phi_Full(node, psi, eta, gamma));
  }

  VecRestoreArray(vec, &arr);

  return val;
}



