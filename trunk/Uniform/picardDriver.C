
#include "petsc.h"
#include "petscdmmg.h"

#include "fullDomainUtils.h"
#include "fatBoundaryUtils.h"
#include "mesh.h"
#include "equation_systems.h"

#define __MAIN__
#include "global.h"

int main(int argc, char** argv) {

  LibMeshInit init(argc, argv);

  PetscOptionsInsertFile(MPI_COMM_WORLD, "options", PETSC_TRUE);

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
  PetscInt Nc = 17;
  PetscOptionsGetInt(0, "-Nc", &Nc, 0);

  PetscInt nlevels = 1;
  PetscOptionsGetInt(0, "-dmmg_nlevels", &nlevels, 0);

  PetscInt maxIters = 10;
  PetscOptionsGetInt(0, "-max_iters", &maxIters, 0);

  PetscScalar theta = 0.9;
  PetscOptionsGetScalar(0, "-theta", &theta, 0);

  PetscScalar tol = 1.0e-10;
  PetscOptionsGetScalar(0, "-tol", &tol, 0);

  PetscInt Nf = (((Nc - 1) << (nlevels - 1)) + 1);

  DA dac;
  DACreate3d(MPI_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, Nc, Nc, Nc,
      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, PETSC_NULL, &dac);

  DMMG* dmmg;
  DMMGCreate(MPI_COMM_WORLD, nlevels, NULL, &dmmg);
  DMMGSetDM(dmmg, (DM) dac);
  DADestroy(dac);

  //Create Full Domain solution vector (Initial Guess)
  Vec solFull = DMMGGetx(dmmg);
  Vec rhsFull = DMMGGetRHS(dmmg);
  Vec solFullTmp;
  VecDuplicate(solFull, &solFullTmp);

  VecZeroEntries(solFull);

  //Compute base RHS vec for full domain
  Vec rhsFullBase;
  VecDuplicate(rhsFull, &rhsFullBase);

  computeRHSterm1_Full(DMMGGetDA(dmmg), rhsFullBase);

  //Call DMMGSetKSP and pass CreateMat. Set RHS function = NULL 
  DMMGSetKSP(dmmg, PETSC_NULL, createMatrix_Full); 

  //Create Neumann Mat, base RHS vec for fat boundary
  Mat neumannMatFat;
  Vec rhsFatBase;
  createNeumannMatrix_Fat(neumannMatFat, rhsFatBase, system, dof_map, mesh) ;

  //Duplicate Neumann Mat and base RHS vec for fat boundary
  Mat matFat;
  Vec rhsFat;
  MatDuplicate(neumannMatFat, MAT_COPY_VALUES, &matFat);
  VecDuplicate(rhsFatBase, &rhsFat);

  //Call DirichletMatrixCorrection for bnd=1 and bnd=2 
  dirichletMatCorrection_Fat( matFat, dof_map, mesh, 1 );      
  dirichletMatCorrection_Fat( matFat, dof_map, mesh, 2 );

  //Create Fat boundary solution vector
  Vec solFat;
  VecDuplicate(rhsFatBase, &solFat);

  //Create KSP for fat boundary 
  KSP kspFat;
  KSPCreate(MPI_COMM_WORLD, &kspFat);
  KSPSetOperators(kspFat, matFat, matFat, SAME_PRECONDITIONER);
  KSPSetOptionsPrefix(kspFat, "fat_");
  KSPSetFromOptions(kspFat);
  KSPSetUp(kspFat);

  //Picard block
  int iter = 0;
  for(; iter < maxIters; iter++) {
    std::cout<<"iter = "<<iter<<std::endl;
    //Store for later
    VecCopy(solFull, solFullTmp);

    //Call DirichletAddCorrection2 for fat boundary and add it to base RHS
    dirichletVecAddCorrection2_Fat( neumannMatFat, rhsFat, solFull, Nf, dof_map, mesh);
    VecAXPY(rhsFat, 1, rhsFatBase);

    //Call DirichletSetCorrection1 and DirichletSetCorrection2 for fat boundary
    dirichletVecSetCorrection1_Fat(rhsFat, dof_map, mesh) ;
    dirichletVecSetCorrection2_Fat(rhsFat, solFull, Nf, dof_map, mesh);

    //Solve Fat boundary problem 
    KSPSolve(kspFat, rhsFat, solFat);

    //Add dirac delta corrections to base RHS for full domain
    getDiracFunctions_Fat(system, solFat, rhsFull, mesh, Nf);
    VecAXPY(rhsFull, 1, rhsFullBase);

    //Call DirichletSetCorrection for RHS of full domain
    setDirichletValues_Full(DMMGGetDA(dmmg), rhsFull);

    //Solve Full domain
    DMMGSolve(dmmg);

    //Update solFull
    VecAXPBY(solFull, theta, (1.0 - theta), solFullTmp);

    //Check convergence
    double diffNorm;
    VecAXPY(solFullTmp, -1, solFull);
    VecNorm(solFullTmp, NORM_2, &diffNorm);
    if(diffNorm < tol) {
      break;
    }
  }//end for Picard block


  VecDestroy(solFat);
  VecDestroy(rhsFatBase);
  VecDestroy(rhsFat);

  MatDestroy(matFat);
  MatDestroy(neumannMatFat);

  VecDestroy(rhsFullBase);
  VecDestroy(solFullTmp);

  KSPDestroy(kspFat);

  DMMGDestroy(dmmg);

  mesh.clear();

  system.clear();

  equation_systems.clear();

  fatBoundary.clear();

}

