
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

double gradPhi(int node, int coord, double psi, double eta, double gamma) {

  double val = 0;

  if(node == 0) { 
    if(coord == 0) {
      val = -(1.0 - eta)*(1.0 - gamma)/8.0;
    } else if(coord == 1) {
      val = -(1.0 - psi)*(1.0 - gamma)/8.0;
    } else if(coord == 2) {
      val = -(1.0 - psi)*(1.0 - eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 1) {
    if(coord == 0) {
      val = (1.0 - eta)*(1.0 - gamma)/8.0;
    } else if(coord == 1) {
      val = -(1.0 + psi)*(1.0 - gamma)/8.0;
    } else if(coord == 2) {
      val = -(1.0 + psi)*(1.0 - eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 2) {
    if(coord == 0) {
      val = -(1.0 + eta)*(1.0 - gamma)/8.0;
    } else if(coord == 1) {
      val = (1.0 - psi)*(1.0 - gamma)/8.0;
    } else if(coord == 2) {
      val = -(1.0 - psi)*(1.0 + eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 3) {
    if(coord == 0) {
      val = (1.0 + eta)*(1.0 - gamma)/8.0;
    } else if(coord == 1) {
      val = (1.0 + psi)*(1.0 - gamma)/8.0;
    } else if(coord == 2) {
      val = -(1.0 + psi)*(1.0 + eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 4) {
    if(coord == 0) {
      val = -(1.0 - eta)*(1.0 + gamma)/8.0;
    } else if(coord == 1) {
      val = -(1.0 - psi)*(1.0 + gamma)/8.0;
    } else if(coord == 2) {
      val = (1.0 - psi)*(1.0 - eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 5) {
    if(coord == 0) {
      val = (1.0 - eta)*(1.0 + gamma)/8.0;
    } else if(coord == 1) {
      val = -(1.0 + psi)*(1.0 + gamma)/8.0;
    } else if(coord == 2) {
      val = (1.0 + psi)*(1.0 - eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 6) {
    if(coord == 0) {
      val = -(1.0 + eta)*(1.0 + gamma)/8.0;
    } else if(coord == 1) {
      val = (1.0 - psi)*(1.0 + gamma)/8.0;
    } else if(coord == 2) {
      val = (1.0 - psi)*(1.0 + eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 7) {
    if(coord == 0) {
      val = (1.0 + eta)*(1.0 + gamma)/8.0;
    } else if(coord == 1) {
      val = (1.0 + psi)*(1.0 + gamma)/8.0;
    } else if(coord == 2) {
      val = (1.0 + psi)*(1.0 + eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else {
    std::cout<<"Error!!!"<<std::endl;
    exit(1);
  }

  return val;

}

double integrand(int node1, int node2, double psi, double eta, double gamma) {
  double val = (gradPhi(node1, 0, psi, eta, gamma)*gradPhi(node2, 0, psi, eta, gamma)) +
    (gradPhi(node1, 1, psi, eta, gamma)*gradPhi(node2, 1, psi, eta, gamma)) +
    (gradPhi(node1, 2, psi, eta, gamma)*gradPhi(node2, 2, psi, eta, gamma));

  return val;
}

int main() {

  double gPt[] = { (-1.0/sqrt(3.0)), (1.0/sqrt(3.0)) };

  double mat[8][8];

  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      mat[i][j] = 0.0;
      for(int m = 0; m < 2; m++) {
        for(int n = 0; n < 2; n++) {
          for(int p = 0; p < 2; p++) {
            mat[i][j] += integrand(i, j, gPt[m], gPt[n], gPt[p]);
          }//end for p
        }//end for n
      }//end for m
    }//end for j
  }//end for i

  FILE * fp = fopen("LaplacianStencil.txt", "w");

  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      fprintf(fp, "%.13lf, ", mat[i][j]);
    }//end for j
    fprintf(fp, "\n");
  }//end for i

  fclose(fp);

}

