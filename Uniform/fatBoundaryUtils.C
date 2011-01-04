//----------------------------------------------- 
// AFBM Unstructured Poisson - Rahul & Srikanth  
//----------------------------------------------- 

#include "fatBoundaryUtils.h"
#include "fullDomainUtils.h"

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
#include "boundary_info.h"

#include "global.h"

//create the Poisson Volume Operator
void createNeumannMatrix_Fat(Mat & petscMat)
{

  Mesh fatBoundary (3);

  fatBoundary.read("hollowsphere.e");

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

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      const Real x = __CENTER_X__ + q_point[qp](0);
      const Real y = __CENTER_Y__ + q_point[qp](1);
      const Real z = __CENTER_Z__ + q_point[qp](2);

      // Assemble the Poisson Operator equation 
      for (unsigned int i = 0; i < n_dofs; i++) {
        for (unsigned int j = 0; j < n_dofs; j++) {
          Ke(i,j) += (JxW[qp]*dphi[i][qp]*dphi[j][qp]);
        }
      }

      // Build the element right hand side
      for (unsigned int i = 0; i < phi.size(); i++) {
        Fe(i) += (__FN__(x, y, z)*JxW[qp]*phi[i][qp]);
      }

    } // end of the quadrature point qp-loop

    system.matrix->add_matrix (Ke, dof_indices);  //                 
    system.rhs->add_vector    (Fe, dof_indices);  // 
  } // end of element loop

  system.matrix->close();

  Mat OriginalMat = (dynamic_cast <PetscMatrix<Number>* > (system.matrix))->mat();
  PetscErrorCode  ierr;
  ierr = MatDuplicate( OriginalMat , MAT_COPY_VALUES, &petscMat);
  //      petscVec = (system.rhs)->vec();

}

//DirichletMatrixCorrection
void dirichletMatCorrection_Fat(LinearImplicitSystem& system , Mat & stiffnesMatrix, MeshBase & mesh)
{

  MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 

  const DofMap & dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;

  for( ; el != end_el; ++el) {
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);

    const unsigned int n_dofs   = dof_indices.size();

    for(unsigned int s = 0; s < elem->n_sides(); s++) {
      if ( elem->neighbor(s) == NULL ) {
        AutoPtr<Elem> side (elem->build_side(s));
        for(unsigned int ns = 0; ns < side->n_nodes(); ns++) {
          for(unsigned int n = 0; n < elem->n_nodes(); n++) {

            if( elem->node(n) == side->node(ns) ) { 
              for(unsigned int m = 0; m < elem->n_nodes(); m++) {
                if(m!=n){
                  MatSetValue ( stiffnesMatrix, m , n, 0.0, INSERT_VALUES);
                  MatSetValue ( stiffnesMatrix, n , m, 0.0, INSERT_VALUES);
                }
              }//end for m
              MatSetValue ( stiffnesMatrix, n , n, 1.0, INSERT_VALUES);
            }

          }//end for n
        }//end for ns
      }

    }//end for s
  }//end for el
  system.matrix->close();

}

//DirichletVectorCorrection
void dirichletVecCorrection_Fat(LinearImplicitSystem& system , Vec dirichletVec, MeshBase & mesh)
{

  MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 

  const DofMap & dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;

  PetscScalar dirichletValue;
  int idx;

  for( ; el != end_el; ++el) {
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);

    const unsigned int n_dofs   = dof_indices.size();

    for(unsigned int s = 0; s < elem->n_sides(); s++) {
      if ( elem->neighbor(s) == NULL ) {
        AutoPtr<Elem> side (elem->build_side(s));
        for(unsigned int ns = 0; ns < side->n_nodes(); ns++) {
          for(unsigned int n = 0; n < elem->n_nodes(); n++) {

            idx = dof_indices[n];
            VecGetValues( dirichletVec, 1 , &idx , &dirichletValue );
            if( elem->node(n) == side->node(ns) ) { 
              system.rhs->set(dof_indices[n], dirichletValue );
              //       MatSetValue ( & stiffnesMatrix, n , n, 1.0);
            }

          }//end for n
        }//end for ns
      }

    }//end for s
  }//end for el

}

//ComputeRHSCorrection
void computeRHS_Fat(Mat  & stiffnessMatrix, Vec & dirichletVec)
{
  Vec* correctionVec; 
  VecDuplicate( dirichletVec, correctionVec ); 
  MatMult( stiffnessMatrix, dirichletVec, * correctionVec );
}

//get coordinates of the FatBoundary
void getBoundary_Fat(LinearImplicitSystem& system, std::vector<double> fatBnd, MeshBase & mesh)
{

  MeshBase::const_node_iterator       node     = mesh.local_nodes_begin();
  const MeshBase::const_node_iterator end_node = mesh.local_nodes_end(); 

  const DofMap & dof_map = system.get_dof_map();
  std::vector<unsigned int> dof_indices;

  PetscScalar dirichletValue;
  int idx;

}

void getDiracFunctions_Fat(LinearImplicitSystem& system, Vec & solVec, Vec & rhs, MeshBase & mesh)
{

  MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 

  const unsigned int V_var = system.variable_number ("V");

  FEType fe_face_type = system.variable_type(V_var);

  AutoPtr<FEBase> fe_face  (FEBase::build(3, fe_face_type));

  QGauss qrule (2, fe_face_type.default_quadrature_order());

  fe_face->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe_face->get_JxW();

  const std::vector<Point>& q_point = fe_face->get_xyz();

  const std::vector<std::vector<RealGradient> >& dphi = fe_face->get_dphi();

  const DofMap & dof_map = system.get_dof_map();

  const std::vector<Point>& face_normal = (fe_face->get_normals());

  std::vector<unsigned int> dof_indices;
  int idx, N;
  std::vector<unsigned int> indices;
  std::vector<double> vals;
  PetscScalar solValue;

  VecZeroEntries(rhs);

  for( ; el != end_el; ++el) {

    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs   = dof_indices.size();

    for(unsigned int s = 0; s < elem->n_sides(); s++) {

      fe_face->reinit(elem, s);

      const short int bnd_id = (mesh.boundary_info)->boundary_id (elem, s);
      if(  bnd_id == 2 ){

        for(unsigned int qp = 0; qp < qrule.n_points(); qp++) 
        {
          double px = __CENTER_X__ + q_point[qp](0);
          double py = __CENTER_Y__ + q_point[qp](1);
          double pz = __CENTER_Z__ + q_point[qp](2);

          phi_Full(px, py, pz, N, indices, vals);

          for (unsigned int j=0; j<indices.size(); j++)
          {
            double dVdx = 0.;
            double dVdy = 0.;
            double dVdz = 0;
            double dVdn = 0;

            for (unsigned int l=0; l<dphi.size(); l++)
            {
              idx = dof_indices[l];
              VecGetValues( solVec, 1 , &idx , &solValue );
              dVdx +=  dphi[l][qp](0)*solValue;
              dVdy +=  dphi[l][qp](1)*solValue;
              dVdz +=  dphi[l][qp](2)*solValue;
            }//nodes

            dVdn  = (dVdx*face_normal[qp](0) + dVdy*face_normal[qp](1) + dVdz*face_normal[qp](2) )*vals[j]*JxW[qp];
            VecSetValue(rhs, indices[j], dVdn , ADD_VALUES);
          }//full domain indices
        }//quadrature points

      }// interior fat boundary

    }//elem side iterator
  }// elem iterator

  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);

}


