//----------------------------------------------- 
// AFBM Unstructured Poisson - Rahul & Srikanth  
//----------------------------------------------- 

#include "fatBoundaryUtils.h"

#include <iostream>
#include <cstdlib>

#include "petscmat.h"
#include "petscvec.h"
//------------------------------------------------------
// Basic include file needed for the mesh functionality.
//------------------------------------------------------
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "petsc_matrix.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "elem.h"

//create the Poisson Volume Operator
void createNeumannMatrix_Fat(std::string meshName, Mat petscMat)
{

      Mesh fatBoundary (3);

      fatBoundary.read(meshName);

      EquationSystems equation_systems (fatBoundary);

      LinearImplicitSystem & system = 
        equation_systems.add_system<LinearImplicitSystem> ("Poisson");

      system.add_variable ("V", FIRST);
      equation_systems.init ();

      const MeshBase & mesh = equation_systems.get_mesh();

      const unsigned int V_var = system.variable_number ("V");

      FEType fe_temp_type = system.variable_type(V_var);

      AutoPtr<FEBase> fe_temp  (FEBase::build(3, fe_temp_type));

      QGauss qrule (3, fe_temp_type.default_quadrature_order());

      fe_temp->attach_quadrature_rule (&qrule);

      const std::vector<Real>& JxW = fe_temp->get_JxW();

      const std::vector<Point>& q_point = fe_temp->get_xyz();

      const std::vector<std::vector<RealGradient> >& dphi = fe_temp->get_dphi();

      const std::vector<std::vector<Real> >& phi = fe_temp->get_phi();

      const DofMap & dof_map = system.get_dof_map();

      DenseMatrix<Number> Ke;
      DenseVector<Number> Fe;

      std::vector<unsigned int> dof_indices;

      MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

      for ( ; el != end_el; ++el)
      {    
        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices);

        const unsigned int n_dofs   = dof_indices.size();

        fe_temp->reinit  (elem);

        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          const Real x = q_point[qp](0);
          const Real y = q_point[qp](1);
          const Real z = q_point[qp](2);

          // Assemble the Poisson Operator equation 
          for (unsigned int i=0; i<n_dofs; i++)
            for (unsigned int j=0; j<n_dofs; j++)
              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

          // Build the element right hand side
          for (unsigned int i=0; i<phi.size(); i++)
            Fe(i) += JxW[qp]*phi[i][qp];


        } // end of the quadrature point qp-loop

        system.matrix->add_matrix (Ke, dof_indices);  //                 
        system.rhs->add_vector    (Fe, dof_indices);  // 
      } // end of element loop

      system.matrix->close();

      petscMat = (dynamic_cast <PetscMatrix<Number>* > (system.matrix))->mat();
//      petscVec = (system.rhs)->vec();

}

//DirichletMatrixCorrection
void dirichletMatCorrection_Fat()
{

}

//DirichletVectorCorrection
void dirichletVecCorrection_Fat()
{
}

//ComputeRHSCorrection
void computeRHS_Fat()
{
}

