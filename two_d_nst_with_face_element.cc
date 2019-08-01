//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//Driver for Navier Stokes in backward step domain -- meshed with triangle

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Dev header
#include "new_nst_sing_face_element.h"
#include "fluid_traction_elements.h"

using namespace std;

using namespace oomph;

// Hierher we introduce u



//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 // Dimensionless domain values
 // ---------------------------

 /// Dimless width of inflow channel
 double H_up = 1.0; 

 /// Dimless length of channel upstream of step
 double L_up = 2.0; 

 /// \short Dimless length of channel downstream of step
 double L_down = 2.0; 

 /// Dimless width of outflow channel
 double H_down = 2.0;

 /// Quadratic nromalized inflow profile: x-velocity as a fct of y coordinate
 double quadratic_flow(const double& y)
  {
  // Bottom of inflow is at y=0; top is at H_up
  return 4.0*y*(H_up - y);
  }

 /// Radius of internal boundary (surrounding the singularity)
 double Radius_of_internal_boundary = 0.5;

 // Bit of a hack but it facilitates reuse...
 #include "unstructured_backward_step_mesh.h"

 double u_additional(const Vector<double>& x);

 /// \short "Singular" function (in cartesian) 
 /// and gradient (in polar coordinates, bc it's easier for the stress tensor
 /// computation)
 void singular_fct_and_gradient(const Vector<double>& x,
                                Vector<double>& u,
                                Vector<Vector<double> >& du_dpolar)
 {

  // Definition of the constants
  const double K = 1.0;
  const long double firstLambda = 1.5444837367824622;
  const double alpha = 3.0 * MathematicalConstants::Pi / 4.0;

  // The dynamic viscosity (water's one i.e. 1e-3)
  const double mu = 1;

  // Radius & polar angle (hireher be careful, the L_up term is designed to
  // recentered the singularity on the corner)
  double y = x[1];
  double r = sqrt((x[0]-L_up)*(x[0]-L_up) + y*y);
  double phi = atan2(y, x[0] - L_up);

  // Define another phi, to rotate the corner and get it the right position
  double phi2 = phi - MathematicalConstants::Pi / 4.0;

  // // Stream function computation
  // // Hierher enable this line if you want to get access to it
  // double radius_part = K*pow(r, firstLambda);

  double first_cos = std::cos(firstLambda*phi2)*
                      std::cos((firstLambda - 2.0)*alpha);
  double second_cos = std::cos((firstLambda - 2.0)*phi2)*
                        std::cos(firstLambda*alpha);

  // // Hierher renable if you want to check the validity of \psi
  // double stream_function = radius_part*(first_cos - second_cos);


  // Velocities w.r.t radius and angle
  double u_r = K*pow(r, firstLambda - 1.0)
                *((firstLambda - 2.0)*std::cos(alpha*firstLambda)
                *std::sin((phi2*(firstLambda - 2.0)))
                - firstLambda*std::cos(alpha*(-2.0 + firstLambda))
                *std::sin(phi2*firstLambda));

  double u_phi = - K*firstLambda*pow(r, firstLambda - 1.0)
                  *((first_cos - second_cos));

  double pressure = 0.0, du_dr = 0.0, du_dphi = 0.0, 
          dv_dr = 0.0, dv_dphi = 0.0;

  // Make sure we don't get a NaN!
  if (r != 0.0)
  {

    // Possible solution of the Stokes flow
    // Since the problem is linear, 
    // could every linear combination fit with the finite element solution?
    double pressure1 = 4.0*K*(firstLambda - 1.0)*pow(r, firstLambda - 2.0)
                        *std::cos(alpha*firstLambda)
                        *std::sin(phi2*(firstLambda - 2.0));

    pressure = mu*pressure1;


    // Compute the derivatives of the velocities
    // du_r/dr, du_r/dphi, du_phi/dr, du_phi/dphi
    du_dr = K*(firstLambda - 1.0)*pow(r, firstLambda - 2.0)
                  *((firstLambda - 2.0)*std::cos(alpha*firstLambda)
                  *std::sin(phi2*(firstLambda - 2.0))
                  - firstLambda*std::cos(alpha*(-2.0 + firstLambda))
                  *std::sin(phi2*firstLambda));

    du_dphi = K*pow(r, firstLambda - 1.0)*( - firstLambda*firstLambda
                        *std::cos(firstLambda*phi2)
                        *std::cos(alpha*(firstLambda - 2.0)) 
                      + std::cos(alpha*firstLambda)
                        *std::cos(phi2*(firstLambda - 2.0))*(firstLambda - 2.0)
                        *(firstLambda - 2.0));

    dv_dr = - K*firstLambda*(firstLambda - 1.0)*pow(r, firstLambda - 2.0)
                    *(first_cos - second_cos);

    dv_dphi = - K*firstLambda*pow(r, firstLambda - 1.0)
                    *(std::cos(alpha*firstLambda)
                      *std::sin(phi2*(firstLambda - 2.0))*(firstLambda - 2.0) 
                    - firstLambda*std::sin(firstLambda*phi2)
                      *std::cos(alpha*(firstLambda - 2.0)));
  } 

  double u_x = u_r*std::cos(phi) - u_phi*std::sin(phi);
  double u_y = u_r*std::sin(phi) + u_phi*std::cos(phi);

  // A singular fonction that solves Stokes flow
  u[0] = u_x;
  u[1] = u_y;
  u[2] = pressure;

  // The associated gradient in polar coordinates
  du_dpolar[0][0] = du_dr;
  du_dpolar[0][1] = du_dphi;
  du_dpolar[1][0] = dv_dr;
  du_dpolar[1][1] = dv_dphi;

  // Hierher discard this .B.
  // double dudr = -2.0/3.0*pow(r,-1.0/3.0)*
  //  sin(2.0/3.0 * (phi + MathematicalConstants::Pi/2.0));
  // double dudphi = pow(r,2.0/3.0)*2.0/3.0*
  //  cos(2.0/3.0 * (phi + MathematicalConstants::Pi/2.0));
  
  // du_dx[0][0] = dudr*cos(phi) - 1.0/r*dudphi*sin(phi);
  // du_dx[1][0]  = dudr*sin(phi) + 1.0/r*dudphi*cos(phi);

  // std::fill(u.begin(), u.end(), 2019);
  // std::fill(du_dx[0].begin(), du_dx[0].end(), 0.0);


 }


 /// \short "Singular" function
  Vector<double> singular_fct(const Vector<double>& x)
 {
  // hierher just taylored to dim 2!!
  Vector<double> u(3, 0);
  Vector<Vector<double> > du_dpolar(2, Vector<double>(2, 0.0));
  singular_fct_and_gradient(x, u, du_dpolar);
  return u;
 }

 /// \short Gradient of "singular" function (in polar coordinates)
 Vector<Vector<double> > gradient_of_singular_fct(const Vector<double>& x)
 {
  Vector<double> u(3, 0); // hierher just taylored to dim 2 !!
  Vector<Vector<double> > du_dpolar(2, u);
  singular_fct_and_gradient(x, u, du_dpolar);
  return du_dpolar;
 }

 /// Exact solution
 Vector<double> u_exact(const Vector<double>& x)
 {
  Vector<double> u(3, 0);
  if (!CommandLineArgs::command_line_flag_has_been_set
       ("--suppress_sing_in_exact_soln"))
   {
    Vector<double> temp = singular_fct(x);
    for (unsigned i = 0; i < u.size(); ++i)
    {
      u[i] += temp[i];
    }
   }
  if (CommandLineArgs::command_line_flag_has_been_set
       ("--add_sin_cosh_to_exact_soln"))
   {
    double temp2 = u_additional(x);
    for (unsigned j = 0; j < u.size(); ++j)
    {
      u[j] += temp2;
    }
   }
  return u;
 }

 double u_additional(const Vector<double>& x)
 {
  return (sin(MathematicalConstants::Pi * x[0]/2)
    *sinh(MathematicalConstants::Pi * x[1]/2));
 }


 // Compute stress tensor coeff \T_ij in cartesian at global coordinates x
 // based on the Moffat singular solution
 Vector<Vector<double> > singular_traction_stress_tensor_nst(
                                              const Vector<double> &x)
 {
   // Polar coordinates w.r.t. the origin
   double r = sqrt((x[0]-L_up)*(x[0]-L_up) + x[1]*x[1]);
   double phi = atan2(x[1], x[0] - 2.0);

   // Get the values of the singular function at our current location
   Vector<double> sing_fct_value = singular_fct(x);
   Vector<Vector<double> > grad_u_sing = gradient_of_singular_fct(x);


   double u_x = sing_fct_value[0];
   double u_y = sing_fct_value[1];
   double pressure = sing_fct_value[2];

   // Compute polar velocities
   double u_r = u_x*std::cos(phi) + u_y*std::sin(phi);
   double u_theta = - u_x*std::sin(phi) + u_y*std::cos(phi);




   Vector<Vector<double> > local_stress(2, Vector<double>(2, 0.0));
   double D_rr, D_rt, D_tr, D_tt;
   double D_xx_hat, D_xy_hat, D_yx_hat, D_yy_hat;
   double D_xx, D_xy, D_yx, D_yy;

   // Deformation tensor in polar coordinates
   // r stands for radius and t for theta
   D_rr = 2.0 * grad_u_sing[0][0];
   D_rt = (grad_u_sing[0][1]/r + grad_u_sing[1][0] 
             - u_theta/r);
   D_tr = D_rt;
   D_tt = 2.0 * (grad_u_sing[1][1] + u_r)/r;


   // 'Half' mapping of deformation tensor in cartesian coord
   // We will also need to multiply this matrix by the transpose of
   // the rotation matrix, to get the real cartesian D_ij
   D_xx_hat = std::cos(phi)*D_rr - std::sin(phi)*D_tr;
   D_xy_hat = std::cos(phi)*D_rt - std::sin(phi)*D_tt;
   // D_yx_hat = D_xy_hat;
   D_yx_hat = std::sin(phi)*D_rr + std::cos(phi)*D_tr;
   D_yy_hat = std::sin(phi)*D_rt + std::cos(phi)*D_tt;


   // We finally get our correct strain tensor D_ij in cartesian
   D_xx = std::cos(phi)*D_xx_hat - std::sin(phi)*D_xy_hat;
   D_xy = std::sin(phi)*D_xx_hat + std::cos(phi)*D_xy_hat;
   D_yx = std::cos(phi)*D_yx_hat - std::sin(phi)*D_yy_hat;
   D_yy = std::sin(phi)*D_yx_hat + std::cos(phi)*D_yy_hat;


   // Compute T_ij = -p delta_ij * D_ij
   // local_stress[0][0] = - sing_fct_value[2] + D_xx;
   local_stress[0][0] = - pressure + D_xx;
   local_stress[0][1] = D_xy;
   local_stress[1][0] = D_yx;
   // local_stress[1][1] = - sing_fct_value[2] + D_yy;
   local_stress[1][1] = - pressure + D_yy;

   // Avoid Nan value at the corner
   if (r == 0.0)
   { 
     local_stress[0][0] = 0.0;
     local_stress[0][1] = 0.0;
     local_stress[1][0] = 0.0;
     local_stress[1][1] = 0.0;
   }

   return local_stress;
 }

 /// Traction required at boundary 2
 void prescribed_traction(const double &t,
                          const Vector<double>& x,
                          const Vector<double> &n,
                          Vector<double>& traction)
  {
    Vector<Vector<double> > traction_stress_tensor = 
                              singular_traction_stress_tensor_nst(x);

    traction[0] = traction_stress_tensor[0][0];
    traction[1] = traction_stress_tensor[0][1];

    // Renable if you want to see the effet of the traction amplitude
    // traction[0] = 0.0;
    // traction[1] = 0.0;
  } 


 // Take a x coordinates vector, return the cartesian singular velocities
 Vector<double> moffat_singular_function(const Vector<double>& x)
  {  

    // Definition of the constant
    const double K = 1.0;
    const long double firstLambda = 1.544483736;
    const double alpha = 3.0 * MathematicalConstants::Pi / 4.0;

    // The dynamic viscosity (here water's one)
    const double mu = 1;

    // Radius & polar angle (hierher be careful, the L_up term is designed to
    // recentered the singularity on the corner)
    double y = x[1];
    double r = sqrt((x[0]-L_up)*(x[0]-L_up) + y*y);
    double phi = atan2(y, x[0] - L_up);

    // Define another phi, to rotate the corner and get it the right position
    double phi2 = phi - MathematicalConstants::Pi / 4.0;

    // Stream function
    double first_cos = std::cos(firstLambda*phi2)*
                        std::cos((firstLambda - 2.0)*alpha);
    double second_cos = std::cos((firstLambda - 2.0)*phi2)*
                          std::cos(firstLambda*alpha);

    // double radius_part = K*pow(r, firstLambda);
    // double stream_function = radius_part*(first_cos - second_cos);
    // double stream_function = K*pow(r, firstLambda)*
    // (std::cos((firstLambda - 2.0)*alpha)*std::cos(firstLambda*phi/alpha) 
    //   - std::cos(firstLambda*alpha)*std::cos((firstLambda - 2.0)*phi/alpha));

    // Velocities w.r.t radius and angle
    double u_r = K*pow(r, firstLambda - 1.0)
                  *((firstLambda - 2.0)*std::cos(alpha*firstLambda)
                  *std::sin((phi2*(firstLambda - 2.0)))
                  - firstLambda*std::cos(alpha*(-2.0 + firstLambda))
                  *std::sin(phi2*firstLambda));

    double u_phi = - K*firstLambda*pow(r, firstLambda - 1.0)
                    *((first_cos - second_cos));


    double pressure;

    // Make sure we don't get a NaN!
    if (r != 0.0)
    {
      double pressure1 = 4.0*K*pow(r, firstLambda - 2.0)
                          *std::cos(alpha*firstLambda)
                          *std::sin(phi2*(firstLambda - 2.0))
                          *(firstLambda - 1) 
                        - 2*K*pow(r,(firstLambda - 2.0))
                          *std::sin(firstLambda*phi2)
                          *std::cos(alpha*(firstLambda - 2.0))
                          *(2.0*firstLambda*firstLambda 
                              - 5.0*firstLambda + 2.0);

      // double pressure2 = -(2.0*K*pow(r,firstLambda - 2.0)
      //                     *(2.0*std::cos(alpha*firstLambda)
      //                         *std::sin(phi2*(firstLambda - 2.0)) 
      //                       + 2.0*firstLambda*firstLambda
      //                         *std::sin(firstLambda*phi2)
      //                         *std::cos(alpha*(firstLambda - 2.0)) 
      //                       - firstLambda*std::cos(alpha*firstLambda)
      //                         *std::sin(phi2*(firstLambda - 2.0)) 
      //                       - firstLambda*std::sin(firstLambda*phi2)
      //                         *std::cos(alpha*(firstLambda - 2.0))))
      //                     /(firstLambda - 2.0);

      pressure = mu*pressure1;
    } 
    else 
    {
      pressure = 0.0;
    }


    Vector<double> u_s(4, 0.0);

    // We get back to the cartesian coordinates
    u_s[0] = u_r*std::cos(phi) - u_phi*std::sin(phi);
    u_s[1] = u_r*std::sin(phi) + u_phi*std::cos(phi);
    u_s[2] = pressure;

    // u_s[0] = stream_function;
    // u_s[1] = u_r*std::cos(phi) - u_phi*std::sin(phi);
    // u_s[2] = u_r*std::sin(phi) + u_phi*std::cos(phi);
    // u_s[3] = pressure;

    return u_s;
  }

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// NavierStokes in backward-facing step. Dirichlet boundary conditions along
/// all boundaries.
//====================================================================
template<class ELEMENT>
class StepProblem : public Problem
{

public:


  /// Constructor
  StepProblem();

  /// Destructor 
 ~StepProblem()
  {
   // hierher: at some point delete things properly and do memory leak
   // check .B.
   delete Bulk_mesh_pt->spatial_error_estimator_pt();
   delete Bulk_mesh_pt;
  }
 
 /// Update the after solve (empty)
 void actions_after_newton_solve() {
  delete_face_elements(); // Hierher be carefull!
 }
 
 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
 
 // Perform actions after mesh adaptation
 void actions_after_adapt()
  {
   // Recreate face elements
   // hierher renable create_face_elements();
   
   // Complete problem setup
   complete_problem_setup();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 
 /// Perform actions after mesh adaptation (empty)
 void actions_before_adapt()
  {
   // Kill face elements
   // hierher renable delete_face_elements();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 
 /// Access function for the specific mesh
 RefineableTriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  

  /// \short Create traction elements on boundary b of the Mesh pointed
  /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
  /// surface_mesh_pt
 void create_traction_elements(const unsigned &b, 
                                Mesh* const &bulk_mesh_pt,
                                Mesh* const &surface_mesh_pt);

 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Do what it says
 void complete_problem_setup();
 
 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// sDelete face elements and flush meshes
 void delete_face_elements()
  {
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }
   
   // // hierher 
   // // Wipe the mesh
   // Face_mesh_for_singularity_integral_pt->flush_element_and_node_storage();
   // n_element = Singular_fct_element_mesh_pt->nelement();
   // for(unsigned e=0;e<n_element;e++)
   // {
   // delete Singular_fct_element_mesh_pt->element_pt(e);
   // }

   // Singular_fct_element_mesh_pt->flush_element_and_node_storage();

  }



 // hierher renable!
 /// Create face elements
 void create_face_elements()
  {
   
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }

    // Hierher renable to set the amplitude .B.
    // // We set the amplitude of the sing func with that
    // dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
    //       Singular_fct_element_mesh_pt->element_pt(0))
    //       ->set_amplitude_of_singular_fct(3.0);

   /// Map keeps a running count of duplicate nodes already created;
   /// existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
   std::map<Node*,Node*> existing_duplicate_node_pt;


   // Create the face elements needed to compute the amplitude of
   // the singular function,

   // hierher change this to include the outer boundaries only!
  for(unsigned i_bound = 0; i_bound < 6; i_bound++)
  {
      unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
      for(unsigned e = 0; e < n_element; e++)
      {
       // Create Pointer to bulk element adjacent to the boundary
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
         (Bulk_mesh_pt->boundary_element_pt(i_bound, e));
       
       // Get Face index of boundary in the bulk element
       int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound, e);
       
       // Create corresponding face element
       FluxElementForSingularityEquation<ELEMENT>* flux_element_pt =
        new FluxElementForSingularityEquation<ELEMENT>(bulk_elem_pt,face_index);

        // We pass the pointer of singular function to the face element
       flux_element_pt->navier_stokes_sing_el_pt()
         ->unscaled_singular_fct_pt()
        = &Global_Physical_Variables::singular_fct;

        // We pass the pointer of the gradient of the singular function to 
        // the face element 
       flux_element_pt->navier_stokes_sing_el_pt()
         ->gradient_of_unscaled_singular_fct_pt() =
           &Global_Physical_Variables::gradient_of_singular_fct;

       // Attach it to the mesh
       Face_mesh_for_singularity_integral_pt->add_element_pt(flux_element_pt);
    }
  }

  // Update the pointer to the face elements
  dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
   Singular_fct_element_mesh_pt->element_pt(0))->
  set_mesh_of_face_elements(Face_mesh_for_singularity_integral_pt);

  
  // Be cautious with this code, can be removed .B.
  
   // Now loop over bulk elements 
   //----------------------------------------------------
   // and swap over any of their nodes have been replaced
   //----------------------------------------------------
   unsigned n_el = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < n_el; e++)
    {
     ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->element_pt(e));
   
     // Loop over all nodes and check if they're amongst the replaced
     // ones
     unsigned nnod=bulk_el_pt->nnode();
     for (unsigned j = 0; j < nnod; j++)
      {
       Node* nod_pt = bulk_el_pt->node_pt(j);
     
       // Find this original node in the map; if we find it
       // it's already been duplicated
       std::map<Node*,Node*>::iterator it =
        existing_duplicate_node_pt.find(nod_pt);
       if (it != existing_duplicate_node_pt.end())
        {
         // Use the existing duplicate node
         bulk_el_pt->node_pt(j) = (*it).second;
        }
      }   
    }
}


// Hierher mesh definition

 /// Pointer to the bulk mesh
 RefineableTriangleMesh<ELEMENT> *Bulk_mesh_pt;

 /// \short Face elements used to compute the amplitude of the singular
 /// function
 Mesh* Face_mesh_for_singularity_integral_pt;
 
 /// Mesh for (single) element containing singular fct
 Mesh* Singular_fct_element_mesh_pt;

 /// Pointer to the "surface" mesh (for traction elements)
 Mesh* Surface_mesh_pt;

 /// \short Enumeration for IDs of FaceElements (used to figure out
 /// who's added what additional nodal data...)
 enum{Flux_jump_el_id, BC_el_id};
 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for StepProblem problem
//========================================================================
template<class ELEMENT>
StepProblem<ELEMENT>::StepProblem()
{

  // Build the mesh
  double uniform_element_area=0.01;
  Bulk_mesh_pt = Global_Physical_Variables::build_the_mesh<ELEMENT>
   (uniform_element_area);

  // // Let's have a look at the boundary enumeration (don't need this)
  // Bulk_mesh_pt->output_boundaries("RESLT/boundaries.dat");

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
  
  // Set element size limits
  Bulk_mesh_pt->max_element_size() = 0.01;
  Bulk_mesh_pt->min_element_size() = 1e-30;
  Bulk_mesh_pt->max_permitted_error() = 0.005;
  Bulk_mesh_pt->min_permitted_error() = 0.0;

  // Hierher traction elements implementation (dev)
  // Create "surface mesh" that will contain only the prescribed-traction 
  // elements. The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_pt = new Mesh;

  // Create prescribed-traction elements from all elements that are 
  // adjacent to boundary 2, but add them to a separate mesh.
  create_traction_elements(2, Bulk_mesh_pt, Surface_mesh_pt);



  // Add sub-mesh
  add_sub_mesh(Bulk_mesh_pt);
  // Add traction elements sub mesh 
  add_sub_mesh(Surface_mesh_pt);

  // hierher renable
  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    
    // Create element that stores the singular fct and its amplitude
    // ---------------------------------------------------------------
    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt=
     new ScalableSingularityForNavierStokesElement<ELEMENT>;

     // cout << "Pointer of singularity is ready: \033[7;41m" 
     // << el_pt << "\033[0m\n" << std::endl;
    
    // Pass fct pointers:
    el_pt->unscaled_singular_fct_pt()
     = &Global_Physical_Variables::singular_fct;
    el_pt->gradient_of_unscaled_singular_fct_pt() =
     &Global_Physical_Variables::gradient_of_singular_fct;

    // Add to mesh
    Singular_fct_element_mesh_pt = new Mesh;
    Singular_fct_element_mesh_pt->add_element_pt(el_pt);
    add_sub_mesh(Singular_fct_element_mesh_pt); // reach oomph error

    Face_mesh_for_singularity_integral_pt = new Mesh; 
    create_face_elements();
    
    add_sub_mesh(Face_mesh_for_singularity_integral_pt);
    // // I think it is not needed

    }


  // Build global mesh
  build_global_mesh();
  
  // Complete problem setup
  complete_problem_setup();
    
  // Setup equation numbering scheme
  oomph_info << "Number of equations: " 
             << this->assign_eqn_numbers() 
             << std::endl;
  
} // end_of_constructor


//==start_of_complete======================================================
/// Set boundary condition, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::complete_problem_setup()
{
 
 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   
   // Loop over the elements to set up element-specific
   // things that cannot be handled by constructor
  
   // hierher renable this one first!
   // Bulk elements
   unsigned n_el = Bulk_mesh_pt->nelement();
   for (unsigned e = 0; e < n_el; e++)
    {
     ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->element_pt(e));
     
     // Tell the bulk element about the singular fct
     bulk_el_pt->navier_stokes_sing_el_pt() =
      dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
       Singular_fct_element_mesh_pt->element_pt(0));
    }
  }
 
 // Apply bcs
 apply_boundary_conditions();

 /// Assign traction elements values
 // Loop over the flux elements to pass pointer to prescribed traction function
 // and pointer to global time object
 unsigned n_traction_el = Surface_mesh_pt->nelement();

 for(unsigned e = 0; e < n_traction_el; ++e)
  {

    // Upcast from GeneralisedElement to traction element
    NavierStokesTractionElement<ELEMENT> *el_pt = 
     dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
      Surface_mesh_pt->element_pt(e));
    // Set the pointer to the prescribed traction function
    el_pt->traction_fct_pt() = 
     &Global_Physical_Variables::prescribed_traction;

  }
 // How many traction elements?
 cout << "Number of traction elements: " << n_traction_el << std::endl;
 
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::apply_boundary_conditions()
{
 
 //Set the boundary conditions for this problem: All nodes are
 //free by default -- just pin the ones that have Dirichlet conditions
 //here.
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound = 0; ibound < num_bound; ibound++)
  {   
   // Leave internal boundary alone
   if (ibound != 6 && ibound != 2)
    { 
     unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod = 0; inod < num_nod; inod++)
      {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);
       Vector<double> x_position(2);
       x_position[0] = nod_pt->x(0); 
       x_position[1] = nod_pt->x(1); 
       
       // if (!( (x_position[0] == 4.0)
       //        && (x_position[1] == 1.0 || x_position[1] == -1.0) ))
       //  {
         Vector<double> singular_velocities;
         singular_velocities = Global_Physical_Variables
                                  ::singular_fct(x_position);

         nod_pt->pin(0);
         nod_pt->pin(1);
         nod_pt->set_value(0, singular_velocities[0]);
         nod_pt->set_value(1, singular_velocities[1]);

       // }


       // // x velocity is pinned at and on no slip boundaries
       // // i.e. everywhere apart from outflow boundary 2
       // if (ibound != 2)
       //  {
       //   Bulk_mesh_pt->boundary_node_pt(ibound, inod)->pin(0);
       //  }
       
       // // y velocity is zero either because of no slip or
       // // horizontal flow
       // Bulk_mesh_pt->boundary_node_pt(ibound, inod)->pin(1);

      }
    }
  } // end loop over boundaries



//  Hierher renable
 // // Now set boundary values on inflow
 // unsigned ibound = 0;
 // unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
 // for (unsigned inod = 0; inod < num_nod; inod++) 
 // {
 //   Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);
 //   double y = nod_pt->x(1);
 //   double u = Global_Physical_Variables::quadratic_flow(y);
 //   nod_pt->set_value(0, u);
 //  }



  // // Hierher quick test: impose u = u_s on all boundaries != 6
  // unsigned num_bound = Bulk_mesh_pt->nboundary();
  // for(unsigned ibound = 0; ibound < num_bound; ibound++)
  //  {   
  //   // Leave internal boundary alone
  //   if (ibound != 6)
  //    { 
  //     unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
  //     for (unsigned inod = 0; inod < num_nod; inod++)
  //      {
  //       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);

  //       Vector<double> x_position(2);
  //       x_position[0] = nod_pt->x(0); 
  //       x_position[1] = nod_pt->x(1); 

  //       // if (!( (x_position[0] == Global_Physical_Variables::L_up)
  //       //        && (x_position[1] == 0.0) ))
  //       {
  //         Vector<double> singular_velocities;
  //         singular_velocities = Global_Physical_Variables::singular_fct(
  //                                   x_position);

  //         double phi = atan2(x_position[1], 
  //                             x_position[0] - Global_Physical_Variables::L_up);

  //         double u_x = singular_velocities[0]*std::cos(phi) 
  //                     - singular_velocities[1]*std::sin(phi);
  //         double u_y = singular_velocities[0]*std::sin(phi) 
  //                     + singular_velocities[1]*std::cos(phi);


          
  //         if (!( (x_position[0]  == Global_Physical_Variables::L_up)
  //                && (x_position[1] == 0.0) ))
  //          {
  //           nod_pt->pin(0);
  //          }
  //         else
  //          {
  //           oomph_info << "at the corner...\n";
  //          }
  //         nod_pt->pin(1);

  //         nod_pt->set_value(0, u_x);
  //         nod_pt->set_value(1, u_y);
  //       }
  //       // else {
  //       //   // nod_pt->pin(0);
  //       //   // nod_pt->pin(1);

  //       //   // nod_pt->set_value(0, 0.0);
  //       //   // nod_pt->set_value(1, 0.0);
  //       // }


  //      }
  //    }
  //  } // end loop over boundaries


  // hierher
  // Assign solution everywhere
  { 
   unsigned num_nod = Bulk_mesh_pt->nnode();
   for (unsigned inod = 0; inod < num_nod; inod++)
    {
     Node* nod_pt = Bulk_mesh_pt->node_pt(inod);
     
     Vector<double> x_position(2);
     x_position[0] = nod_pt->x(0); 
     x_position[1] = nod_pt->x(1); 
     
     Vector<double> singular_velocities;
     singular_velocities = Global_Physical_Variables::singular_fct(x_position);

     nod_pt->set_value(0, singular_velocities[0]);
     nod_pt->set_value(1, singular_velocities[1]);

     
     if (nod_pt->nvalue() == 3)
      {
       double p = 0.0; // dummy infinity 
       if (!( (x_position[0] == Global_Physical_Variables::L_up)
              && (x_position[1] == 0.0) ))
        {
         p = singular_velocities[2];
        }
       nod_pt->set_value(2, p);
      }
     
    }
  }

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 10;

  // Output solution
  sprintf(filename, "%s/soln%i.dat", doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();


  // Output solution just using vertices so we can see the mesh
  sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  npts = 2;
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();
  

  // // hierher renable (?)
  // // Plot "extended solution" showing contributions; also work out
  // // average element size
  // double av_el_size=0.0;
  // sprintf(filename,"%s/extended_soln%i.dat",doc_info.directory().c_str(),
  //         doc_info.number());
  // some_file.open(filename);
  // unsigned nel=Bulk_mesh_pt->nelement();
  // for (unsigned e=0;e<nel;e++)
  //  {
  //   npts=20;
  //   ELEMENT* el_pt=
  //    dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
  //   el_pt->output_with_various_contributions(some_file,npts);
  //   av_el_size+=el_pt->size();
  // }
  // some_file.close();
  // av_el_size/=double(nel);



  // Doc error norm:
  oomph_info << "\nav el size, av h, Ndof, # bulk els, Norm of error    : "   
// hierher renable             
              // << av_el_size << " " 
             // << sqrt(av_el_size) << " " 
             << ndof() << " " 
             << Bulk_mesh_pt->nelement() << " " 
            // << sqrt(error) <<" "
             << std::endl;
  
  

  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    oomph_info 
     << "Amplitude of singular function: "
     << dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
      Singular_fct_element_mesh_pt->element_pt(0))->
     amplitude_of_singular_fct() << std::endl;

     oomph_info << "Compared with the perimeter i.e.: " << 12.0 << std::endl;


    // Output face elements used to compute amplitude of singularity
    sprintf(filename,"%s/hierher_face_elements%i.dat",
            doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    unsigned nel = Face_mesh_for_singularity_integral_pt->nelement();
    for (unsigned e = 0; e < nel; e++)
     {
      npts = 5;
      FluxElementForSingularityEquation<ELEMENT>* el_pt =
       dynamic_cast<FluxElementForSingularityEquation<ELEMENT>*>
       (Face_mesh_for_singularity_integral_pt->element_pt(e));
      //  Vector<double> s(2);
      //  oomph_info << el_pt->stress_tensor_nst(s, 0, 0) << std::endl;
      el_pt->output(some_file, npts);
     }
    some_file.close();

   }



} // end_of_doc_solution


//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the 
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::create_traction_elements(const unsigned &b,
                                                  Mesh* const &bulk_mesh_pt,
                                                  Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e = 0; e < n_element; e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b, e));
   
   // What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

   // Build the corresponding prescribed-flux element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
    NavierStokesTractionElement<ELEMENT>(bulk_elem_pt, face_index);
   
   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);
  } //end of loop over bulk elements adjacent to boundary b
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for backward step with impedance outflow bc
//=====================================================================
int main(int argc, char **argv)
{
  // Little test hierher (should be discarded)
  double N_test(200.0);
  ofstream some_file;
  char filename[100];

  // Output the analytical solution on a grid 1 x 1 centered on the corner
  sprintf(filename,"RESLT/moffat_singular_function.dat");
  some_file.open(filename);

  some_file << "ZONE I=" << N_test << ", J=" << N_test << std::endl;

  for (double kline = 0.0; kline < N_test; ++kline)
  {
    for (double kcolumn = 0.0; kcolumn < N_test; ++kcolumn)
    {
      Vector<double> x(2);

      // Column index (begin at 1, end at 3)
      x[0] = 2*(kcolumn/N_test - 0.5) + 2.0;
      // Line index
      x[1] = 2*(- kline/N_test + 0.5);

      some_file << x[0] << " " << x[1] << " ";

      Vector<double> u_s(4, 0.0);
      if (x[0] <= 2.0 && x[1] < 0)
      {
        some_file << 0 << " " << 0 << " " << std::endl;
      }
      else {
        u_s = Global_Physical_Variables::moffat_singular_function(x);
        some_file << u_s[0] << " " << u_s[1] << " " << std::endl;
      }
      
    }
  }

  some_file.close();

 // Store command line arguments
 CommandLineArgs::setup(argc, argv);
  
 // Don't subtract off singularity
 CommandLineArgs::specify_command_line_flag(
  "--dont_use_singularity");

 // Use sin cosh in exact solution
 CommandLineArgs::specify_command_line_flag(
  "--add_sin_cosh_to_exact_soln");

 // Suppress singular term in exact solution
 CommandLineArgs::specify_command_line_flag(
  "--suppress_sing_in_exact_soln");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
 // Set up doc info
 // ---------------
 
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 

  // Build the problem with 
 // StepProblem<ProjectableTaylorHoodElement<TTaylorHoodElement<2> > > 
  StepProblem<ProjectableTaylorHoodElement<MyTTaylorHoodElement<2> > >
  problem;

  
 // Doc initial guess
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Solve, refine uniformly and keep going
  unsigned max_adapt = 1;
  for (unsigned i = 0; i < max_adapt; i++)
   {
    // Solve the bloody thing
    problem.newton_solve();
    problem.doc_solution(doc_info);
    doc_info.number()++;
    if (i != (max_adapt - 1)) problem.refine_uniformly();
   }

}

