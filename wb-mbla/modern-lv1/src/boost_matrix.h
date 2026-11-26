
#ifndef BOOST_MATRIX_H
#define BOOST_MATRIX_H

/*
  This is collection of MATRIX which will be used the whole calculation.
  Without this collection, we have to allocate matrix all the time during the calculation, which will be deleted after the lifetime of the matrix.

 */

#include "parameters.h"
#include "matrix.h"

class AUX_MATRIX
{
public:
  MATRIX phi_yt_div_phib_arr;
  MATRIX Ieta_yt_arr;
  MATRIX ID_yt_arr;
  MATRIX shear_rate_nodim; // for shear-induced migration
  MATRIX shear_rate_nodim_pre; // for shear-induced migration (using the previous iteration step)
  MATRIX shear_stress_nodim;
  MATRIX shear_stress_nodim_pre;

  MATRIX A_mat, y_mat, d_mat, R_arr, Q_arr, f_arr;
  MATRIX P_div_DLP_arr, uz_div_u0_arr, vw_div_vw0_arr, Ln_arr;

  MATRIX report_outlet_profile; // it will record the last iteration at the outlet profile
  // MATRIX tmp_diffusivity_yt_arr;

  MATRIX F_info_arr;
  
  AUX_MATRIX()
  {
    cout << "ERROR: AUX_MATRIX class constructor needs PARAMTER_SET as input\n";
  }

  AUX_MATRIX(PARAMETER_SET& PARS)
  {
    // related with Ny
    // related with Ny. But now the information should be stored for all Nz.
    // therefore, the below two matrix now has Ny x Nz dimension
    phi_yt_div_phib_arr.initial(PARS.Nz, PARS.Ny, 0.);
    Ieta_yt_arr.initial(PARS.Nz, PARS.Ny, 0.);
    ID_yt_arr.initial(PARS.Nz, PARS.Ny, 0.);

    /*
      Note that the characteristic shear rate is its expected (or actual) value at the entrance-wall. 
      If u(0,0) is specified as the boundary condition -- its not a currernt version --,
      the characteristic shear rate is 2 u(0,0)/R.
      If P(0) is specified as the boundary condition -- current version --,
      the characteristic shear rate is 2 u^HP / R (see PARAMETER_SET initializer)
     */
    shear_rate_nodim.initial(PARS.Nz, PARS.Ny, 0.); // new member variable
    
    shear_rate_nodim_pre.initial(PARS.Nz, PARS.Ny, 0.); // new member variable

    shear_stress_nodim.initial(PARS.Nz, PARS.Ny, 0.); // new member variable
    
    shear_stress_nodim_pre.initial(PARS.Nz, PARS.Ny, 0.); // new member variable

    report_outlet_profile.initial(PARS.Ny, 6, 0.); // yt, phi_yt, dg, tau, diffusivity, viscosity
    
    // related with Nz
    P_div_DLP_arr.initial(PARS.Nz, 1, 0.);
    uz_div_u0_arr.initial(PARS.Nz, 1, 0.);
    vw_div_vw0_arr.initial(PARS.Nz, 1, 0.);
    Ln_arr.initial(PARS.Nz, 1, 0.);
    
    F_info_arr.initial(PARS.Nz, 4, 0.); // z, F1, F2, F2_modi

    // related with N_eff
    // long N_eff = PARS.Nz - 2;
    long N_eff = PARS.Nz;
    
    A_mat.initial(N_eff, N_eff, 0.);
    y_mat.initial(N_eff, 1, 0.);
    d_mat.initial(N_eff, 1, 0.);

    R_arr.initial(PARS.Nz, 1, 0.);
    Q_arr.initial(PARS.Nz, 1, 0.);
    f_arr.initial(PARS.Nz, 1, 0.);

    
  }
    
  ~AUX_MATRIX()
  {}
};

#endif
