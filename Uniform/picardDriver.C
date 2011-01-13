
#include "petsc.h"
#include "petscdmmg.h"

#include "fullDomainUtils.h"
#include "fatBoundaryUtils.h"
#include "mesh.h"
#include "equation_systems.h"

#define __MAIN__
#include "global.h"

int main(int argc, char** argv) {

  PetscInitialize(&argc, &argv, "options", NULL);
  LibMeshInit init(argc, argv);

  //Init stuff for Fat Boundary (Read mesh, build equation systems, dof-map)  
  Mesh fatBoundary (3);

  fatBoundary.read("hollowsphere.e");

  EquationSystems equation_systems (fatBoundary);

  LinearImplicitSystem & system =  equation_systems.add_system<LinearImplicitSystem> ("Poisson");
  system.add_variable ("V", FIRST);
  equation_systems.init ();

  DofMap & dof_map = system.get_dof_map();

  MeshBase & mesh = equation_systems.get_mesh();

  //Init stuff for Full domain (Create DA/DMMG) 
  PetscInt N = 17;
  PetscOptionsGetInt(0, "-N", &N, 0);

  DA da;
  DACreate3d(MPI_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, N, N, N,
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da);

  PetscInt nlevels = 1;
  PetscOptionsGetInt(0, "-dmmg_nlevels", &nlevels, 0);

  DMMG* dmmg;
  DMMGCreate(MPI_COMM_WORLD, nlevels, NULL, &dmmg);
  DMMGSetDM(dmmg, (DM) da);

  //Create Full Domain solution vector (Initial Guess)
  Vec solFull = DMMGGetx(dmmg);
  Vec rhsFull = DMMGGetRHS(dmmg);

  VecZeroEntries(solFull);

  //Compute base RHS vec for full domain
  Vec rhsFullBase;
  VecDuplicate(rhsFull, &rhsFullBase);

  computeRHSterm1_Full(da, rhsFullBase);

  //Call DMMGSetKSP and pass CreateMat. Set RHS function = NULL 
  DMMGSetKSP(dmmg, PETSC_NULL, createMatrix_Full); 

  //Create Neumann Mat, base RHS vec for fat boundary
  Mat neumannMatFat;
  Vec rhsFat;
  createNeumannMatrix_Fat(neumannMatFat, rhsFat, system, dof_map, mesh) ;

  //Duplicate Neumann Mat and base RHS vec for fat boundary
  Mat matFat;
  Vec rhsFatBase;
  MatDuplicate( neumannMatFat, MAT_COPY_VALUES, &matFat);
  VecDuplicate( rhsVec_Fat, &duplicateVec);

  //Create Fat boundary solution vector
  Vec solutionVec_Fat;
  VecDuplicate( rhsVec_Fat, &solutionVec_Fat);

  //Create KSP for fat boundary 

  //Picard block
  {

    //Call DirichletAddCorrection2 for fat boundary and add it to base RHS
    Vec rhsVecCorrection_Fat;
    VecDuplicate( rhsVec_Fat, &rhsVecCorrection_Fat);
    dirichletVecAddCorrection2_Fat( duplicateMat, rhsVecCorrection_Fat, solutionVec_Fat, N, dof_map, mesh);
    VecAXPY(rhsVec_Fat, 1, rhsVecCorrection_Fat);

    //Call DirichletMatrixCorrection for bnd=1 and bnd=2 
    dirichletMatCorrection_Fat( neumannMat_Fat, dof_map, mesh, 1 );      
    dirichletMatCorrection_Fat( neumannMat_Fat, dof_map, mesh, 2 );

    //Call DirichletSetCorrection1 and DirichletSetCorrection2 for fat boundary
    dirichletVecSetCorrection1_Fat(rhsVec_Fat, dof_map, mesh) ;
    dirichletVecSetCorrection2_Fat(rhsVec_Fat, solFull, N, dof_map, mesh) ;

    //Solve Fat boundary problem 

    //Add dirac delta corrections to base RHS for full domain
    getDiracFunctions_Fat(system, solutionVec_Fat, rhsVec_Fat, mesh);

    //Call DirichletSetCorrection for RHS of full domain

    //Solve Full domain

  }//end for Picar block

  VecDestroy();
  VecDestroy(rhsFullBase);

  DMMGDestroy(dmmg);
  DADestroy(da);

  PetscFinalize();

}

