
#include "petsc.h"
#include "libmesh.h"

int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, NULL, NULL);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(rank == 0) {
    LibMeshInit init(argc, argv, PETSC_COMM_SELF);
  }

  PetscFinalize();

}

