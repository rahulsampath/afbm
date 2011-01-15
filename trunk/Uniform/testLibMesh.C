
#include "petsc.h"
#include "mesh.h"
#include "equation_systems.h"
#include "linear_implicit_system.h"

int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, NULL, NULL);
  LibMeshInit init(argc, argv, PETSC_COMM_SELF);

  //Init stuff for Fat Boundary (Read mesh, build equation systems, dof-map)  
  Mesh fatBoundary (3);

  fatBoundary.read("hollowsphere.e");

  EquationSystems equation_systems (fatBoundary);

  LinearImplicitSystem & system =  equation_systems.add_system<LinearImplicitSystem> ("Poisson");
  system.add_variable ("V", FIRST);
  equation_systems.init ();

  DofMap & dof_map = system.get_dof_map();

  MeshBase & mesh = equation_systems.get_mesh();

  PetscFinalize();

}

