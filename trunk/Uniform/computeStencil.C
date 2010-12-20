
#include <cmath>
#include <cstdlib>
#include <iostream>

double dPhi(int node, int coord, double psi, double eta, double gamma) {

  double val = 0;

  if(node == 0) { 
    if(coord == 0) {
      val = -1.0*(1.0 - eta)*(1.0 - gamma)/8.0;
    } else if(coord == 1) {
      val = -1.0*(1.0 - psi)*(1.0 - gamma)/8.0;
    } else if(coord == 2) {
      val = -1.0*(1.0 - psi)*(1.0 - eta)/8.0;
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 1) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 2) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 3) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 4) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 5) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 6) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
    } else {
      std::cout<<"Error!!!"<<std::endl;
      exit(1);
    }
  } else if(node == 7) {
    if(coord == 0) {
    } else if(coord == 1) {
    } else if(coord == 2) {
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

int main() {

}

