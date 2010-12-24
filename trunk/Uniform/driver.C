
#include "petsc.h"

#include "fullDomainUtils.h"
#include "fatBoundaryUtils.h"

#define __MAIN__
#include "global.h"

int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, "options", NULL);
  
  PetscFinalize();

}

