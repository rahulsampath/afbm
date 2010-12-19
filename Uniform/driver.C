
#include "petsc.h"

int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, "options", NULL);

  PetscFinalize();

}

