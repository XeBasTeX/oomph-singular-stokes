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
// Header file for elements that are used to ...hierher
#ifndef OOMPH_NAVIER_STOKES_SING_FACE_ELEMENTS_HEADER
#define OOMPH_NAVIER_STOKES_SING_FACE_ELEMENTS_HEADER

namespace oomph
{

 //============================================================================
 /// TemplateFreeScalableSingularityForNavierStokesElement defines the
 /// elements managing the singular function : it is essentially a pointer
 /// to the singular function, its gradient and its amplitude
 //============================================================================
 class TemplateFreeScalableSingularityForNavierStokesElement :
  public virtual GeneralisedElement
 {
   public:

  /// Typedef for function that returns the singular solutinon u[i], where
  /// the entries are u,v,[w],p
  typedef Vector<double>(*UnscaledSingSolnFctPt) (const Vector<double>& x);
  

  /// Typedef for gradient of singular function:
  /// returned_vector[i][j]=d u_i/d x_j
  typedef Vector<Vector<double> >(*GradientOfUnscaledSingSolnFctPt) 
   (const Vector<double>& x);
  
  /// Constructor
  TemplateFreeScalableSingularityForNavierStokesElement()
  {
  add_internal_data(new Data(1)); // data to store amplitude
  }

  /// Destructor
  ~TemplateFreeScalableSingularityForNavierStokesElement()
  {
    /* Nothing here, should investigate on memory leak risk*/ 
  }

  /// Function to get pointer to unscaled version of singular function
  UnscaledSingSolnFctPt& unscaled_singular_fct_pt()
   {return Unscaled_singular_fct_pt;}
  
  /// \short Function to get pointer to unscaled version of gradient of
  /// singular function
  GradientOfUnscaledSingSolnFctPt& gradient_of_unscaled_singular_fct_pt() 
   {return Gradient_of_unscaled_singular_fct_pt;}
  
  /// \short Function to compute unscaled version of singular fct u[i], where
  /// the entries are u,v,[w],p
  Vector<double> unscaled_singular_fct(const Vector<double>& x) const
   {
    if(Unscaled_singular_fct_pt == 0)
     {
      throw OomphLibError(
       "The pointer hasn't be set",
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
      abort();
     }
    return Unscaled_singular_fct_pt(x);
   }
  
  /// \short Compute unscaled version of gradient of singular function
  /// returned_vector[i][j]=d u_i/d x_j
  Vector<Vector<double> > gradient_of_unscaled_singular_fct(
   const Vector<double>& x) const
   {
    Vector<double> grad;
    if(Gradient_of_unscaled_singular_fct_pt == 0)
     {
      throw OomphLibError(
       "The pointer hasn't be set",
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
      abort();
     }
    return Gradient_of_unscaled_singular_fct_pt(x);
   }
  
  /// \short Compute scaled version of singular function u[i], where
  /// the entries are u,v,[w],p
  Vector<double> singular_fct(const Vector<double>& x) const
  {
    Vector<double> z;
    for (unsigned i = 0; i < unscaled_singular_fct(x).size(); ++i)
    {
      z.push_back(amplitude_of_singular_fct()*unscaled_singular_fct(x)[i]);
    }
    return z;
  }
  
  /// \short Compute scaled version of gradient of singular function
  /// returned_vector[i][j]=d u_i/d x_j
  Vector<Vector<double> > gradient_of_singular_fct(
   const Vector<double>& x) const
   {
    // create a 2-dim vector named grad
    Vector<Vector<double> > grad(gradient_of_unscaled_singular_fct(x));
    unsigned n = grad.size();
    for(unsigned i = 0; i < n; ++i)
     {
      unsigned m = grad[i].size();
      for (unsigned j = 0; i < m; ++j)
      {
        grad[i][j]*= amplitude_of_singular_fct();
      }
     }
    return grad;
    }

    ///Access the amplitude of the singular function
    double amplitude_of_singular_fct() const
    {
      return data_that_stores_amplitude_of_singular_fct()
        ->value(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Set the amplitude of thz singular function
    void set_amplitude_of_singular_fct(const double& value)
    {
      data_that_stores_amplitude_of_singular_fct()
        ->set_value(index_of_value_that_stores_amplitude_of_singular_fct(),
                    value);
    }

    /// Pin amplitude of singular function
    void pin_amplitude_of_singular_fct()
    {
      data_that_stores_amplitude_of_singular_fct()
        ->pin(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Pointer to data that stores the amplitude of singular function
    Data* data_that_stores_amplitude_of_singular_fct() const
    {
      return internal_data_pt(0);
    }

    ///Gives the index of the amplitude value : default is 0
    unsigned index_of_value_that_stores_amplitude_of_singular_fct() const 
    {
      return 0;
    }

  private:

    /// Pointer to singular function
    UnscaledSingSolnFctPt Unscaled_singular_fct_pt;

    /// Pointer to gradient of singular function
    GradientOfUnscaledSingSolnFctPt Gradient_of_unscaled_singular_fct_pt;
  };

 
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// FluxElementForSingularityEquation is a class of face elements 
  ///used to compute the contribution to the residuals from the the singular function
  // hierher shouldn't be called flux and why are they called Equation?
  //===========================================================================
  template <class ELEMENT>
  class FluxElementForSingularityEquation :
   public virtual FaceGeometry<ELEMENT>, 
  public virtual FaceElement
  {
   
  public:
    //Pointer to compute singular function related stuff
    TemplateFreeScalableSingularityForNavierStokesElement*& navier_stokes_sing_el_pt() 
    { 
      return Navier_stokes_sing_el_pt; 
    } 

   /// \short Function pointer to the prescribed-flux function fct(x,f(x)) -- 
   /// x is a Vector! 
    typedef void (*NavierStokesPrescribedFluxFctPt) // hierher shouldn't be called flux!
    (const Vector<double>& x, double& flux);

   /// \short Constructor, takes the pointer to the "bulk" element and the 
   /// index of the face to which the element is attached.
   FluxElementForSingularityEquation(FiniteElement* const &bulk_el_pt, 
                      const int& face_index);
   
   ///\short  Broken empty constructor
   FluxElementForSingularityEquation()
    {
     throw OomphLibError(
      "Don't call empty constructor for FluxElementForSingularityEquation",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }

   /// Broken copy constructor
   FluxElementForSingularityEquation(const FluxElementForSingularityEquation&
                                     dummy) 
    { 
     BrokenCopy::broken_copy("FluxElementForSingularityEquation");
    } 
   
   /// Broken assignment operator
   void operator=(const FluxElementForSingularityEquation&) 
    {
     BrokenCopy::broken_assign("FluxElementForSingularityEquation");
    }

   /// \short Specify the value of nodal zeta from the face geometry
   /// The "global" intrinsic coordinate of the element when
   /// viewed as part of a geometric object should be given by
   /// the FaceElement representation, by default (needed to break
   /// indeterminacy if bulk element is SolidElement)
   double zeta_nodal(const unsigned &n, const unsigned &k,           
                            const unsigned &i) const 
    {return FaceElement::zeta_nodal(n, k, i);}

   /// Add the element's contribution to its residual vector
   inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
     // Call the generic residuals function with flag set to 0
     // using a dummy matrix argument
     fill_in_generic_residual_contribution_navier_stokes_flux(
      residuals, GeneralisedElement::Dummy_matrix, 0);
    }


   /// \short Add the element's contribution to its residual vector and its 
   /// Jacobian matrix
   inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian)
    {
     // Call the generic routine with the flag set to 1
     fill_in_generic_residual_contribution_navier_stokes_flux(residuals,
                                                              jacobian, 1);
    }

    // // Former Stress tensor hierher, can discard
    // // Hierher usefull function for stress tensor's computation
    // // We should improve the computation with a case!
    // unsigned kronecker_delta(const unsigned i, const unsigned j) {
    //   if(i == j) {
    //     return 1;
    //   }
    //   else {
    //     return 0;
    //   }
    // }

    // Transform polar velocities in cartesian coordinates
    Vector<double> velocities_from_polar_to_cartesian(const double &phi,
                                                        const double &u_r,
                                                        const double &u_phi) 
                                                          const
      {

        double u_x = u_r*std::cos(phi) - u_phi*std::sin(phi);
        double u_y = u_r*std::sin(phi) + u_phi*std::cos(phi);

        Vector<double> cartesian_velocities(Dim);
        cartesian_velocities[0] = u_x;
        cartesian_velocities[1] = u_y;

        return cartesian_velocities;
      }


    // Transform cartesian velocities in polar coordinates
    Vector<double> velocities_from_cartesian_to_polar(const double &phi,
                                                        const double &u_x,
                                                        const double &u_y) 
                                                          const
      {

        double u_r = u_x*std::cos(phi) + u_y*std::sin(phi);
        double u_theta = - u_x*std::sin(phi) + u_y*std::cos(phi);

        Vector<double> polar_velocities(Dim);
        polar_velocities[0] = u_r;
        polar_velocities[1] = u_theta;

        return polar_velocities;
      }


    // Compute stress tensor coeff \T_ij at local coordinates s
    // based on the finite element solution
    Vector<Vector<double> > fe_stress_tensor_nst(const Vector<double> &s)
                               const
      {

        // // Mu = 1 (dynamic viscosity)
        // double mu = 1.0;

        // Get gradient of FE solution from bulk element 
         ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt()); 
         Vector<double> flux(Dim); 
         Vector<double> s_bulk(Dim); 
         s_bulk = local_coordinate_in_bulk(s);

         // Vector<double> x(Dim);
         // for(unsigned i = 0; i < Dim; i++) 
         //  {
         //   x[i] = this->interpolated_x(s, i);
         //  }

          Vector<Vector<double> > local_stress(2, Vector<double>(2, 0.0));

          // compute T_ij with finite element data

          // local_stress[0][0] = - bulk_el_pt->TTaylorHoodElement<2>
          //                         ::interpolated_p_nst(s_bulk)
          //                       + 2*bulk_el_pt->TTaylorHoodElement<2>::
          //                         interpolated_dudx_nst(s_bulk, 0, 0);

          local_stress[0][0] = 2.0*bulk_el_pt->TTaylorHoodElement<2>::
                                  interpolated_dudx_nst(s_bulk, 0, 0);



          local_stress[0][1] = bulk_el_pt->TTaylorHoodElement<2>::
                                  interpolated_dudx_nst(s_bulk, 0, 1) 
                                + bulk_el_pt->TTaylorHoodElement<2>::
                                  interpolated_dudx_nst(s_bulk, 1, 0);

          local_stress[1][0] = local_stress[0][1];

          // local_stress[1][1] = - bulk_el_pt->TTaylorHoodElement<2>
          //                         ::interpolated_p_nst(s_bulk) 
          //                       + 2*bulk_el_pt->TTaylorHoodElement<2>
          //                         ::interpolated_dudx_nst(s_bulk, 1, 1);

          local_stress[1][1] = 2.0*bulk_el_pt->TTaylorHoodElement<2>
                                  ::interpolated_dudx_nst(s_bulk, 1, 1) ;

         return local_stress;
      }

  // Compute stress tensor coeff \T_ij at local coordinates s
  // based on the Moffat singular solution
  Vector<Vector<double> > singular_stress_tensor_nst(const Vector<double> &s)
                            const
    {

      // // Get gradient of FE solution from bulk element 
      //  ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt()); 
      //  Vector<double> flux(Dim); 
      //  Vector<double> s_bulk(Dim); 
      //  s_bulk = local_coordinate_in_bulk(s);


       Vector<double> x(Dim);
       for(unsigned i = 0; i < Dim; i++) 
        {
         x[i] = this->interpolated_x(s, i);
        }

        // Polar coordinates w.r.t. the origin
        double r = sqrt((x[0]-2.0) * (x[0]-2.0) + x[1]*x[1]);
        double phi = atan2(x[1], x[0] - 2.0);

        // Get the values of the singular function at our current location
        Vector<double> sing_fct_value =
         Navier_stokes_sing_el_pt->unscaled_singular_fct(x);


        double u_x = sing_fct_value[0];
        double u_y = sing_fct_value[1];
        double pressure = sing_fct_value[2];

        // Compute polar velocities
        double u_r = u_x*std::cos(phi) + u_y*std::sin(phi);
        double u_theta = - u_x*std::sin(phi) + u_y*std::cos(phi);

        // Get gradient of singular function
        // Hierher we use a little trick to initialize grad_u_sing 
        // with sing_fct_value
        Vector<Vector<double> > grad_u_sing(Dim, sing_fct_value);
        grad_u_sing = Navier_stokes_sing_el_pt
                      ->gradient_of_unscaled_singular_fct(x);


       // Should make this multi-dim (3D) compatible
       // stackoverflow.com/questions/24580714/why-is-this-not-a-constant-expression



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

   /// Output function -- forward to broken version in FiniteElement
   /// until somebody decides what exactly they want to plot here...
   void output(std::ostream &outfile) {FiniteElement::output(outfile);}

   /// Output with various contributions
   void  output(std::ostream &outfile, 
                const unsigned &nplot)
    {
     // Vector of local coordinates
     Vector<double> s(Dim - 1);
     
     // Tecplot header info
     outfile << this->tecplot_zone_string(nplot);
     
     // Loop over plot points
     unsigned num_plot_points = this->nplot_points(nplot);
     for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
       // Get local coordinates of plot point
       this->get_s_plot(iplot, nplot, s);
       
       
       Vector<double> x(Dim);
       for(unsigned i = 0; i < Dim; i++) 
        {
         x[i] = this->interpolated_x(s, i);
         outfile << x[i] << " ";
        }

       // Compute outer unit normal at the specified local coordinate
       // Vector<double> unit_normal(Dim);
       // outer_unit_normal(s, unit_normal);
       // for(unsigned i = 0; i < Dim; i++) 
       //  {
       //   outfile << unit_normal[i] << " ";
       //  }
       

       // Vector<double> sing_fct_value =
       //  Navier_stokes_sing_el_pt->unscaled_singular_fct(x);

        // we compute phi w.r.t. the corner (hence the minus 2)
        // double phi = atan2(x[1], x[0] - 2.0);

        // Vector<double> sing_fct_cartesian = 
        //                   velocities_from_polar_to_cartesian(phi, sing_fct_value[0], sing_fct_value[1]);

        // Vector<Vector<double> > sing_grad_value =
        //  Navier_stokes_sing_el_pt->gradient_of_unscaled_singular_fct(x);

         Vector<Vector<double> > T_s, T_fe;
         T_s = singular_stress_tensor_nst(s);
         T_fe = fe_stress_tensor_nst(s);

         // Error in norm 1 i.e. \sum abs(T_fe - T_s)

         outfile << "" << T_fe[1][1] << " " << T_s[1][1];

         // outfile << "" << std::abs(T_fe[0][0] - T_s[0][0]) 
         //          + std::abs(T_fe[0][1] - T_s[0][1])
         //          + std::abs(T_fe[1][0] - T_s[1][0])
         //          + std::abs(T_fe[1][1] - T_s[1][1]);

         // outfile << " && " << T_fe[0][0] << " " << T_fe[0][1] 
         //          << " " << T_fe[1][0] << " " << T_fe[1][1];
         // outfile << " && " << T_s[0][0] << " " << T_s[0][1]
         //          << " " << T_s[1][0] << " " << T_s[1][1];
          // outfile << " && " << T_fe[0][0] - T_s[0][0] 
          //          << " " << T_fe[0][1] - T_s[0][1]
          //          << " " << T_fe[1][0] - T_s[1][0]
          //          << " " << T_fe[1][1] - T_s[1][1];

        // End the outstream by skipping a line
        outfile << std::endl;
      }
     
     // Write tecplot footer (e.g. FE connectivity lists)
     this->write_tecplot_zone_footer(outfile, nplot);
     
    }


   /// C-style output function -- forward to broken version in FiniteElement
   /// until somebody decides what exactly they want to plot here...
   void output(FILE* file_pt) {FiniteElement::output(file_pt);}

   /// \short C-style output function -- forward to broken version in 
   /// FiniteElement until somebody decides what exactly they want to plot 
   /// here...
   void output(FILE* file_pt, const unsigned &n_plot)
    {FiniteElement::output(file_pt, n_plot);}

   double get_contribution_integral();




  protected:

   /// \short Function to compute the shape and test functions and to return 
   /// the Jacobian of mapping between local and global (Eulerian)
   /// coordinates
   inline double shape_and_test(const Vector<double> &s, Shape &psi,
                                Shape &test)
    const
    {
     //Find number of nodes
     unsigned n_node = nnode();

     //Get the shape functions
     shape(s, psi);

     //Set the test functions to be the same as the shape functions
     for(unsigned i = 0; i < n_node; i++) {test[i] = psi[i];}

     //Return the value of the jacobian
     return J_eulerian(s);
    }


   /// \short Function to compute the shape and test functions and to return 
   /// the Jacobian of mapping between local and global (Eulerian)
   /// coordinates
   inline double shape_and_test_at_knot(const unsigned &ipt,
                                        Shape &psi, Shape &test)
    const
    {
     //Find number of nodes
     unsigned n_node = nnode();

     //Get the shape functions
     shape_at_knot(ipt, psi);

     //Set the test functions to be the same as the shape functions
     for(unsigned i = 0; i < n_node; i++) {test[i] = psi[i];}

     //Return the value of the jacobian
     return J_eulerian_at_knot(ipt);
    }

  private:

   /// \short Add the element's contribution to its residual vector.
   /// flag=1(or 0): do (or don't) compute the contribution to the
   /// Jacobian as well. 
   void fill_in_generic_residual_contribution_navier_stokes_flux(
    Vector<double> &residuals, DenseMatrix<double> &jacobian, 
    const unsigned& flag);

   /// The spatial dimension of the problem
   unsigned Dim;

   /// The index at which the unknown is stored at the nodes
   /// Hierher we changed U_index_nst to U_index_nst!!! Be cautious
   unsigned U_index_nst;

   /// hierher
   TemplateFreeScalableSingularityForNavierStokesElement*
    Navier_stokes_sing_el_pt;
  }; 

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  FluxElementForSingularityEquation<ELEMENT>::
  FluxElementForSingularityEquation(FiniteElement* const &bulk_el_pt, 
                     const int &face_index) : 
    FaceGeometry<ELEMENT>(), FaceElement()
    { 
     // Let the bulk element build the FaceElement, i.e. setup the pointers 
     // to its nodes (by referring to the appropriate nodes in the bulk
     // element), etc.
     bulk_el_pt->build_face_element(face_index, this);

     // Initialising the pointer to the singularity function
     this->Navier_stokes_sing_el_pt = 
      new TemplateFreeScalableSingularityForNavierStokesElement;
   
  #ifdef PARANOID
     {
      //Check that the element is not a refineable 3d element
      ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
      //If it's three-d
      if(elem_pt->dim() == 3)
       {
        //Is it refineable
        RefineableElement* ref_el_pt =
         dynamic_cast<RefineableElement*>(elem_pt);
        if(ref_el_pt != 0)
         {
          if (this->has_hanging_nodes())
           {
            throw OomphLibError(
             "This flux element will not work correctly if nodes are hanging\n",
             OOMPH_CURRENT_FUNCTION,
             OOMPH_EXCEPTION_LOCATION);
           }
         }
       }
     }
  #endif

     // Extract the dimension of the problem from the dimension of 
     // the first node
     Dim = this->node_pt(0)->ndim();

     //Set up U_index_nst. Initialise to zero, which probably won't change
     //in most cases, oh well, the price we pay for generality
     U_index_nst = 0;

     //Cast to the appropriate NavierStokesEquation so that we can
     //find the index at which the variable is stored
     //We assume that the dimension of the full problem is the same
     //as the dimension of the node, if this is not the case you will have
     //to write custom elements, sorry
     switch(Dim)
      {
       //One dimensional problem
      case 1:
      {
       NavierStokesEquations<1>* eqn_pt = 
        dynamic_cast<NavierStokesEquations<1>*>(bulk_el_pt);
       //If the cast has failed die
       if(eqn_pt == 0)
        {
         std::string error_string =
          "Bulk element must inherit from NavierStokesEquations.";
         error_string += 
          "Nodes are one dimensional, but cannot cast the bulk element to\n";
         error_string += "NavierStokesEquations<1>\n.";
         error_string += 
          "If you desire this functionality, you must implement it yourself\n";
         
         throw OomphLibError(error_string,
                             OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
        }
       //Otherwise read out the value
       else
        {
         //Read the index from the (cast) bulk element
         U_index_nst = eqn_pt->u_index_nst(face_index);
        }
      }
      break;
      
      //Two dimensional problem
      case 2:
      {
       NavierStokesEquations<2>* eqn_pt = 
        dynamic_cast<NavierStokesEquations<2>*>(bulk_el_pt);
       //If the cast has failed die
       if(eqn_pt == 0)
        {
         std::string error_string =
          "Bulk element must inherit from NavierStokesEquations.";
         error_string += 
          "Nodes are two dimensional, but cannot cast the bulk element to\n";
         error_string += "NavierStokesEquations<2>\n.";
         error_string += 
          "If you desire this functionality, you must implement it yourself\n";
         
         throw OomphLibError(error_string,
                             OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
       else
        {
         //Read the index from the (cast) bulk element.
         U_index_nst = eqn_pt->u_index_nst(face_index);
        }
      }
      break;
      
      //Three dimensional problem
      case 3:
      {
       NavierStokesEquations<3>* eqn_pt = 
        dynamic_cast<NavierStokesEquations<3>*>(bulk_el_pt);
       //If the cast has failed die
       if(eqn_pt == 0)
        {
         std::string error_string =
          "Bulk element must inherit from NavierStokesEquations.";
         error_string += 
          "Nodes are three dimensional, but cannot cast the bulk element to\n";
         error_string += "NavierStokesEquations<3>\n.";
         error_string += 
          "If you desire this functionality, you must implement it yourself\n";
         
         throw OomphLibError(error_string,
                             OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
         
        }
       else
        {
         //Read the index from the (cast) bulk element.
         U_index_nst = eqn_pt->u_index_nst(face_index);
        }
      }
      break;

      //Any other case is an error
      default:
       std::ostringstream error_stream; 
       error_stream <<  "Dimension of node is " << Dim 
                    << ". It should be 1,2, or 3!" << std::endl;
       
       throw OomphLibError(error_stream.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
       break;
      }
    }



  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
  void FluxElementForSingularityEquation<ELEMENT>::
  fill_in_generic_residual_contribution_navier_stokes_flux(
   Vector<double> &residuals, DenseMatrix<double> &jacobian, 
   const unsigned& flag)
  {
  }

  //===========================================================================
  /// Hierher calculate the contribution of the face element to the integral
  //===========================================================================
  template<class ELEMENT>
  double FluxElementForSingularityEquation<ELEMENT>::
  get_contribution_integral()
  {
   //Find out how many nodes there are
   const unsigned n_node = nnode();
    
   //Set up memory for the shape and test functions
   Shape psif(n_node), testf(n_node);
   
   //Set the value of Nintpt
   const unsigned n_intpt = integral_pt()->nweight();
   
   //Set the Vector to hold local coordinates
   Vector<double> s(Dim-1);
   
   // Saves result of integration
   double integral_result = 0.0;

   // Loop over the integration points
   //--------------------------------
   for(unsigned ipt = 0; ipt < n_intpt; ipt++)
    {

     //Assign values of s
     for(unsigned i = 0; i < (Dim-1); i++)
      {s[i] = integral_pt()->knot(ipt, i);}
     
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Find the shape and test functions and return the Jacobian
     //of the mapping
     double J = shape_and_test(s, psif, testf);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;

     // compute outer normal unit vector
     Vector<double> unit_normal(Dim);
     outer_unit_normal(s, unit_normal);

     // Get the gradient of u_fe and global coordinates
     Vector<double> flux(Dim);
     Vector<double> s_bulk(Dim);
     Vector<double> x(Dim); 

     for(unsigned i = 0; i < Dim; i++)  
        { 
         x[i] = this->interpolated_x(s, i); 
        } 

     // Get the local bulk coordinates    
     s_bulk = local_coordinate_in_bulk(s);

     // dynamic_cast<ELEMENT*>(this->bulk_element_pt())->get_flux(s_bulk,flux);
     // cout << "My Flux is " << flux[0]<<","<< flux[1] <<" I live there : " 
     // <<  x[0] << "," << x[1] << std::endl;

     // URGENT: find the problem with get_flux .B.
     // // Get the Gradient of the FE part of the solution
     // dynamic_cast<ELEMENT*>(this->bulk_element_pt())->TTaylorHoodElement<2>::get_flux(s_bulk, flux);

     // Get the values of the singular function at our current location
     Vector<double> sing_fct_value =
      Navier_stokes_sing_el_pt->unscaled_singular_fct(x);

     // Get gradient of singular function
     // Hierher we use a little trick to initialize grad_u_sing with u_sing
     Vector<Vector<double> > grad_u_sing(Dim, sing_fct_value);
     grad_u_sing = Navier_stokes_sing_el_pt->
                                    gradient_of_unscaled_singular_fct(x);

      // Hierher should discard this (former computation of perimeter)
      integral_result += W*1.0;


      // Hierher machinery of the residual
     // Now we compute the contribution of the integral
     for(unsigned i = 0; i < Dim; i++)
     {
      for (unsigned j = 0; j < sing_fct_value.size(); ++j)
      {
    //     // Hierher renable .B.
    //     integral_result += W*(unit_normal[i]*flux[i]*sing_fct_value[j]
    //       - unit_normal[i]*grad_u_sing[i][j]*
    //         (dynamic_cast<ELEMENT*>(this->bulk_element_pt())
    //           ->raw_interpolated_u_nst(s_bulk))[i]);
    //           // hierher former raw_interpolated_u_poisson
      }
     } 
    }
   return integral_result;
  }



  //======================================================================
  /// \short Class for elements that handle singularities
  /// in NavierStokes equations. Templated by bulk element within
  /// which we impose regularity on the FE solution by insisting that
  /// the slope of the solution at a specified local coordinate, and in
  /// in a specified direction is zero. Nodal values of that element
  /// become external data for the current element whose equation
  /// (zero slope of the FE solution, as discussed) determines the 
  /// amplitude of the singular function.
  //======================================================================
  template<class BULK_ELEMENT> 
  class ScalableSingularityForNavierStokesElement : 
    public virtual TemplateFreeScalableSingularityForNavierStokesElement
   {
   
     public:

    /// Constructor
     ScalableSingularityForNavierStokesElement() :
    Bulk_element_pt(0), Face_element_mesh_pt(0)
     {
     }


    /// Set pointer to mesh containing the FaceElements (and flush
    /// the previuos ones first!)
    void set_mesh_of_face_elements(Mesh* const& face_mesh_pt)
     {
      Face_element_mesh_pt = face_mesh_pt;
      flush_external_data();
      unsigned nel = face_mesh_pt->nelement();
      oomph_info << "Number of face elements: " << nel << std::endl;
      for (unsigned e = 0; e < nel; e++)
       {
        FiniteElement* el_pt =
         dynamic_cast<FluxElementForSingularityEquation<BULK_ELEMENT>*>(
          face_mesh_pt->element_pt(e))->bulk_element_pt();
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
         {
          add_external_data(el_pt->node_pt(j));
         }
       }
     }
    
    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
     //Call the generic residuals function with flag set to 0
     //using a dummy matrix argument
     fill_in_generic_residual_contribution_nst_sing_fct(
      residuals, GeneralisedElement::Dummy_matrix, 0);
    }
    
     private:

    /// Add the element's contribution to its residual vector
    inline void fill_in_generic_residual_contribution_nst_sing_fct(
     Vector<double> &residuals,
     DenseMatrix<double> &jacobian,
     const unsigned& flag)
    {

     if (ndof() == 0) return;

  #ifdef PARANOID
     
  #endif
     
     // hierher paranoid check null pointers and zero sized vectors 
     int c_local_eqn = internal_local_eqn(0, 0);
     residuals[c_local_eqn] = 0.0;
     if (c_local_eqn >= 0)
      {
       unsigned n_element = Face_element_mesh_pt->nelement();

       residuals[c_local_eqn] = amplitude_of_singular_fct();

       for(unsigned e = 0; e < n_element; e++)
       {
        //hierher calculate contribution of each face element 
        //and add it to the residuals
        residuals[c_local_eqn] +=
           - dynamic_cast<FluxElementForSingularityEquation<BULK_ELEMENT>*>
            (Face_element_mesh_pt->finite_element_pt(e))
              ->get_contribution_integral();
       }
      }
    }

     private:
    
    /// Pointer to bulk element where FE solution is regularised
    BULK_ELEMENT* Bulk_element_pt;

    /// Pointer to mesh of face elements that contribute to the surface
    /// integral that determines the amplitude of the unkown function
    Mesh* Face_element_mesh_pt;
    
   };

   ////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////

   // hierher really need to tidy this up! Should only need one class 
   // for T and Q
   //
   //====================================================================
   /// New class. Mainly overloads output-related functions to add
   /// "singular function" (which is assumed to satisfy the Stokes
   /// equation; therefore no change to the governing (bulk) equations) 
   /// to the FE solution. 
   //====================================================================
   template<unsigned DIM>
   class MyTTaylorHoodElement : public virtual TTaylorHoodElement<DIM>
   {

   public:

    /// Constructor
    MyTTaylorHoodElement():Navier_stokes_sing_el_pt(0)
     {
     }

    /// \short Return FE representation of function value u_nst(s) 
    /// plus scaled singular fct (if provided) at local coordinate s
    inline Vector<double> interpolated_u_nst(const Vector<double> &s) const
     {
      Vector<double> u_fe;
      for (unsigned j = 0; j < DIM + 1; ++j)
      {
        u_fe[j] = TTaylorHoodElement<DIM>::interpolated_u_nst(s, j);
        // hierher renable
         if (Navier_stokes_sing_el_pt != 0) 
          {      
           Vector<double> x(DIM); 
           for(unsigned i = 0; i < DIM; i++)  
            { 
             x[i] = this->interpolated_x(s, i); 
            } 
            double sing_value_at_j_component =
             Navier_stokes_sing_el_pt->singular_fct(x)[j];
            u_fe[j] += sing_value_at_j_component; 
          } 
        // hierher end renable
      }
      
      return u_fe;
     } 

    /// \short Return FE representation of function value u_fe
    inline Vector<double> raw_interpolated_u_nst(const Vector<double> &s) const
     {
       Vector<double> u_fe;
       for (unsigned j = 0; j < DIM + 1; ++j)
       {
        u_fe[j] = TTaylorHoodElement<DIM>::interpolated_u_nst(s, j);
        } 
       return u_fe;
     } 

     // This one gives a tecplot output. For the moment: nevermind

    /// Output with various contributions
    void  output_with_various_contributions(std::ostream &outfile, 
                                            const unsigned &nplot)
     {
      //Vector of local coordinates
      Vector<double> s(DIM);
      
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
      
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
       {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);
        
        Vector<double> x(DIM);
        for(unsigned i = 0; i < DIM; i++) 
         {
          x[i] = this->interpolated_x(s,i);
          outfile << x[i] << " ";
         }
        Vector<Vector<double> > u_sing(DIM + 1, DIM);
        // hierher renable
        if (Navier_stokes_sing_el_pt != 0) 
          {
           u_sing = Navier_stokes_sing_el_pt->gradient_of_singular_fct(x); 
          }
        // hierher end renable

        // Currently we only return u-velocity component. Should be update
        outfile << this->interpolated_u_nst(s) << " "
                << TTaylorHoodElement<DIM>::interpolated_u_nst(s) << " "
                << u_sing[0][0] << " "
                << std::endl;   
       }
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
      
     }


    // hierher renable
     /// Pointer to element that stores singular fct 
     TemplateFreeScalableSingularityForNavierStokesElement*& 
     navier_stokes_sing_el_pt() 
     { 
       return Navier_stokes_sing_el_pt; 
     } 
    // hierher end renable

    

   private:

     /// Pointer to element that stores singular fct 
     TemplateFreeScalableSingularityForNavierStokesElement* Navier_stokes_sing_el_pt; 
    
   };


   ////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////



    //=======================================================================
    /// Face geometry for the MyTTaylorHoodElement elements: The spatial
    /// dimension of the face elements is one lower than that of the
    /// bulk element but they have the same number of points
    /// along their 1D edges.
    //=======================================================================
    template<unsigned DIM>
    class FaceGeometry<MyTTaylorHoodElement<DIM> >:
     public virtual TElement<DIM-1,3>
    {

      public:
     
     /// \short Constructor: Call the constructor for the
     /// appropriate lower-dimensional TElement
     FaceGeometry() : TElement<DIM-1,3>() {}

    };

    // ////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////////////////////////////////////


    // //=======================================================================
    // /// Face geometry for the 1D MyTTaylorHoodElement elements: Point elements
    // //=======================================================================
    // template<unsigned NNODE_1D>
    // class FaceGeometry<MyTTaylorHoodElement<1> >: 
    //  public virtual PointElement
    // {

    //   public:
     
    //  /// \short Constructor: Call the constructor for the
    //  /// appropriate lower-dimensional TElement
    //  FaceGeometry() : PointElement() {} 

    // };




 }
#endif
