
#include "ported_lib.h"


int gen_z_div_L_arr(MATRIX& z_div_L_arr, PARAMETER_SET& PARS)
{
  double dz_div_L = PARS.dz / PARS.L_channel;
  z_div_L_arr(0) = 0.;
  for(long i=1; i<PARS.Nz; i++)
    z_div_L_arr(i) = z_div_L_arr(i-1) + dz_div_L;
  return 0;
}

int gen_yt_arr(MATRIX& yt_arr, PARAMETER_SET& PARS)
{
  /*
    NOTE:
    The transversal coordinate here uses as \tilde{y} = y/R. This is of importance for the all the cross-sectional integration where the bulk dimensionless transversal coordinate is selected INSTEAD OF the stretched coordinate. 
   */
  double dy = PARS.dr;
  double dyt = dy/PARS.R_channel;

  double log10_dyp_st = -6;
  double d_log10_dyp = 6./double(PARS.Nr_bulk-1);
  
  PARS.Ny = PARS.Nr_bulk;
  yt_arr.initial(PARS.Ny, 1, 0.);
  yt_arr(0) = pow(10., log10_dyp_st);
  
  double tmp_log10_dyp = log10_dyp_st;
  for (long i=1; i<PARS.Ny; i++)
    {
      tmp_log10_dyp += d_log10_dyp;
      yt_arr(i) = pow(10., tmp_log10_dyp);
    }

  double tmp_yt_0 = yt_arr(0);
  for (long i=0; i<PARS.Ny; i++)
    {
      // make zero at the wall
      yt_arr(i) -= tmp_yt_0;
    }
  
  return 0;

}



int
update_Pi_div_DLP_arr(MATRIX& Pi_div_DLP_arr, MATRIX& phiw_set, MATRIX& shear_rate_nodim_NOT_USED, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  /*
    Note that we do not account the shear-rate dependence on Pi (osmotic pressure)
    since this will be only the matter for the larger particles where the osmotic pressure contribution is negligible.
    
    For future, the appropriate pressure contribution will be addressed here for some specific applications.
   */
  double sf = PARS.Pi_ast/PARS.DLP;
  for (long i=0; i<PARS.Nz; i++)
    {
      Pi_div_DLP_arr(i) = DISP.fcn_Pi_div_Piast(phiw_set(i)*PARS.phi_bulk, PARS, shear_rate_nodim_NOT_USED(i, 0))*sf;
    }
  return 0;
}

double REL_ERR_data(MATRIX& data_set_2, MATRIX& data_set_1, long skip_rows)
{
  double re = 0.;
  for (long i=skip_rows; i<data_set_2.size; i++)
    {
      re += pow(1. - data_set_2(i)/data_set_1(i), 2.0);
    }
  return sqrt(re/double(data_set_2.size)); // divided by number of elements?
}

double ABS_ERR_data(MATRIX& data_set_2, MATRIX& data_set_1, long skip_rows)
{
  double re = 0.;
  for (long i=skip_rows; i<data_set_2.size; i++)
    {
      re += pow(data_set_2(i) - data_set_1(i), 2.0);
    }
  return sqrt(re/double(data_set_2.size - skip_rows)); // divided by number of elements?
}



int get_df_from_f(MATRIX& df_arr, MATRIX& z_div_L_arr, MATRIX& f_arr)
{
  long Nz = df_arr.size;
  double h = z_div_L_arr(1) - z_div_L_arr(0);
  df_arr(0) = (f_arr(1) - f_arr(0))/h;			// forward difference
  df_arr(Nz-1) = (f_arr(Nz-1) - f_arr(Nz-2))/h;		// backward difference
  for (long i=1; i<Nz-1; i++)
    {
      df_arr(i) = (f_arr(i+1) - f_arr(i-1))/(2.*h);	// central difference
    }
  return 0;
}

double get_nodim_Lpf(double f) // f = delta_f_div_R/epsilon_f
{
  return 1./(1. + f);
}

// Hydraulic resistances add in series: Lp_effective/Lp_clean=1/(1+Rcake/Rmem).
// The clean inlet is imposed separately from downstream cake.
int get_nodim_Lpf_arr(MATRIX& Ln_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS)
{
  Ln_arr(0) = 1.; // at the inlet, this value must be the same as Lp0, which emans Lp(z=0)/Lp0 = 1.
  for (long i=1; i<PARS.Nz; i++) 
    {
      Ln_arr(i) = get_nodim_Lpf(Rcake_div_Rmem_arr(i));
    }
  return 0;
}

int init_get_P(MATRIX& A_mat, MATRIX& y_mat, MATRIX &d_mat, MATRIX& R_arr, MATRIX& Q_arr, MATRIX& f_arr, PARAMETER_SET& PARS)
{
  // long N_eff = PARS.Nz - 2;
  // long N_eff = PARS.Nz - 1; // the inlet condition now account also for Neumann boundary conditions.
  long N_eff = PARS.Nz;
  A_mat.initial(N_eff, N_eff, 0.);
  y_mat.initial(N_eff, 1, 0.);
  d_mat.initial(N_eff, 1, 0.);

  R_arr.initial(N_eff, 1, 0.);
  Q_arr.initial(N_eff, 1, 0.);
  f_arr.initial(N_eff, 1, 0.);

  
  return 0;
}



// Dense LU of a tridiagonal axial finite-difference assembly. The final-row
// osmotic RHS term is inherited; see the documented endpoint-pressure limit.
int get_P_div_DLP_arr_with_fouling(MATRIX& P_div_DLP_arr, PARAMETER_SET& PARS, MATRIX& Ln_arr, MATRIX& Pi_div_DLP_arr, MATRIX& A_mat, MATRIX& y_mat, MATRIX &d_mat, MATRIX& R_arr, MATRIX& Q_arr, MATRIX& f_arr, MATRIX& Rcake_div_Rmem_arr)
{
  double h = PARS.dz/PARS.L_channel;
  double tP_in = PARS.Pin/PARS.DLP;
  double tP_out = PARS.Pout/PARS.DLP;
  double tP_perm = PARS.Pper/PARS.DLP;
  

  long np = 0;
  long n;
  for(n=0; n<PARS.Nz; n++)
    {
      // Note that R_arr keeps zero values.
      // double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem_arr(n);
      // double pre_df = pow(1. - cake_thickness_nodim, -3.0);
      double pre_df = 1.;
      Q_arr(n) = pow(PARS.k, 2.0)*Ln_arr(n)*pre_df;
      f_arr(n) = pow(PARS.k, 2.0)*Ln_arr(n)*Pi_div_DLP_arr(n)*pre_df;
    }
  
  double Rn, Qn, fn, an, bn, cn, dn;
  if(PARS.inlet_BC=="velocity")
    {
      // Revised Neumann boundary condition using fictitious node approach
      A_mat(0, 0) = 2. + pow(h, 2.0)*Q_arr(0);
      A_mat(0, 1) = -2.;
      d_mat(0) = 2.*h*1.; // 1 because the dimensionless quantity u(0,0)/u^0 = 1 for inlet_BC=="velocity"

    }
  else
    {
      A_mat(0, 0) = 1.;
      A_mat(0, 1) = 0.;
      d_mat(0) = tP_in - tP_perm;
    }
  
  for(n=1; n<PARS.Nz; n++)
    {
      // np = n - 1;
      np = n; // no reduced space for the inlet 
      Rn = R_arr(np);
      Qn = Q_arr(np);
      fn = f_arr(np);

      an = -1. + (h/2.)*Rn;
      bn =  2. + pow(h, 2.0)*Qn;
      cn = - 1. - (h/2.)*Rn;
      dn = fn * pow(h, 2.0);
      // if (n==1)
      // 	{
      // 	  dn += tP_in - tP_perm; // inlet BC
      // 	  A_mat(np, np) = bn;
      // 	  A_mat(np, np+1) = cn;
      // 	}
      if(n==PARS.Nz - 1)
	{
	  /*
	    Using the reduced space representation for Dirichlet BC
	    Therefore, the matrix will not used the last node n==PARS.Nz - 1
	   */
	  dn += tP_out - tP_perm; // outlet BC
	  // A_mat(np, np-1) = an;
	  A_mat(np, np) = 1.;
	}
      else
	{
	  /*
	    For the interior nodes: 0< np < PARS.Nz - 1
	   */
	  A_mat(np, np-1) = an;
	  A_mat(np, np) = bn;
	  A_mat(np, np+1) = cn;
	}
      d_mat(np) = dn;
    }
  // P_div_DLP_arr(0) = tP_in;
  /*
    Original Python code: P_arr[1:-1] = lin.solve(A_mat, d_mat) + tP_perm
  */
  // gsl_sol_Ax_b(P_div_DLP_arr, A_mat, d_mat, 1);
  /*
    starting index is now 0 instead of 1.
    This is because of the code compatibility for the both of Dirichlet and Neumann boundary conditions. 
    The former BC could be reduced as did for the outlet node, but the latter BC required to be solved.
  */
  gsl_sol_Ax_b(P_div_DLP_arr, A_mat, d_mat, 0);  // this will use the LU decomposition.
  
  for (long i=0; i<PARS.Nz; i++)
    {
      P_div_DLP_arr(i) += tP_perm;
    }
  // P_div_DLP_arr(PARS.Nz - 1) = tP_out;
  
  return 0;
}

int get_uz_div_u0_arr_with_fouling(MATRIX& uz_div_u0_arr, MATRIX& z_div_L_arr, MATRIX& P_div_DLP_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS)
{
  double dP = 0.;
  double h = PARS.dz/PARS.L_channel;
  // for i=0: forward difference
  if (PARS.inlet_BC == "velocity")
    {
      uz_div_u0_arr(0) = 1.;


      // dP = (P_div_DLP_arr(2) - P_div_DLP_arr(0))/(2.*h);       // centeral difference
      // uz_div_u0_arr(1) = -dP;
    }
  else
    {
      dP = (P_div_DLP_arr(1) - P_div_DLP_arr(0))/h;
      uz_div_u0_arr(0) = -dP;

      // dP = (P_div_DLP_arr(2) - P_div_DLP_arr(1))/h;
      // uz_div_u0_arr(1) = -dP;
    }

  // 1<=i<=PARS.Nz-2: centeral difference
  for(long i=1; i<PARS.Nz - 1; i++)
    {
      double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem_arr(i);
      double pre_df = pow(1. - cake_thickness_nodim, -2.0);
      // double pre_df = 1.;
      dP = (P_div_DLP_arr(i+1) - P_div_DLP_arr(i-1))/(2.*h);       // centeral difference
      uz_div_u0_arr(i) = -dP*pre_df;  // dimensionless
    }

  /*
  // old version
  // The previous version used the backward difference scheme
  // for i=Nz-1: backward difference
  dP = (P_div_DLP_arr(PARS.Nz-1) - P_div_DLP_arr(PARS.Nz - 2))/h;
  uz_div_u0_arr(PARS.Nz - 1) = -dP;
  */

  // New version for the outlet used the backward difference scheme with the second order
  // Note that the outlet-pressure boundary condition is used


  double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem_arr(PARS.Nz - 1);
  double pre_df = pow(1. - cake_thickness_nodim, -2.0);
  // double pre_df = 1.;
  uz_div_u0_arr(PARS.Nz - 1) = pre_df*(-1./(2.*h))*(3.*P_div_DLP_arr(PARS.Nz - 1) - 4. * P_div_DLP_arr(PARS.Nz - 2) + P_div_DLP_arr(PARS.Nz - 3));
      
  return 0;
}

int get_vw_div_vw0_arr_with_fouling(MATRIX& vw_div_vw0_arr, MATRIX& P_div_DLP_arr, MATRIX& Pi_div_DLP_arr, MATRIX& Ln_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS)
{
  double tP_perm = PARS.Pper/PARS.DLP;
  for (long n=0; n<PARS.Nz; n++)
    {
      // double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem_arr(n);
      // double pre_df = pow(1. - cake_thickness_nodim, +1.0);
      double pre_df = 1.;
      vw_div_vw0_arr(n) = pre_df*(1./PARS.alpha_ast)*Ln_arr(n)*(P_div_DLP_arr(n) - tP_perm - Pi_div_DLP_arr(n));
    }
  return 0;
}

double get_f_div_fb_modi_RK(double yt, double dyt, double f, double df, double int_INV_D_pre, double vw_div_vw0, double dg_1, double dg_2, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{

  /*
    Due to the precision problem, we made phi_yt array into phi_yt_div_phib array.
    Therefore, the calculation of collective diffusion coefficient must take re-mapping with phi_bulk.

    Note that dg_2 does not mean the full step. It must be adjusted properly (half-step for k2 and k3 in RK 45 method).
  */
  double y_new = yt + dyt;
  double f_new = f + df;
  double int_INV_D = int_INV_D_pre;
  // if (f_new < 0)
  //   cout << "Error: f_new is negative! This is usually happened when Ny is too small. " << f << '\t' << df << '\t' << f_new << endl;
  if (df != 0 ) // it is related with half-step in RK45 method
    int_INV_D += (dyt/2.)*(1./DISP.fcn_Dc_div_D0(f*PARS.phi_bulk, PARS, dg_1) + 1./DISP.fcn_Dc_div_D0(f_new*PARS.phi_bulk, PARS, dg_2));
      
  return (-1./PARS.epsilon_d)*(vw_div_vw0/DISP.fcn_Dc_div_D0(f_new*PARS.phi_bulk, PARS, dg_1))*(f_new - (1. - exp(-(vw_div_vw0/PARS.epsilon_d)*int_INV_D)));
}

int get_phi_yt_div_phib_with_given_phiw(long index_z, MATRIX& phi_yt_div_phib_arr, MATRIX& yt_arr, double given_phiw_div_phib, double vw_div_vw0, MATRIX& shear_rate_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  /*
    Although the matrix dimension has been revised from Ny x 1 to Nz x Ny, the processing of each function is still based on Ny only (that means for the given z).
    Therefore, the treatment here will not be updated unless the index_z for array index.
   */
  // std::cout << given_phiw_div_phib << endl;
  phi_yt_div_phib_arr(index_z, 0) = given_phiw_div_phib; // given value is normalized by phi_bulk
  // printf("%4.3e\n", phi_yt_div_phib_arr(index_z, 0)*PARS.phi_bulk);      

  double int_INV_D_pre = 0.;
  double pre_df = 1.;
  /*
    The new coordinate is from 0 to 1 even with the cake layer.
    Therefore, pre_df is not necessary to account the reduced cross-sectional area.
    Instead, a different place (cross-ectional average) will take account it.
  */

  double y_1 = 0, y_2 = 0., dy = 0., y_h = 0., phi_1 = 0., phi_2 = 0.;
  double k1=0., k2=0., k3=0., k4=0.;

  for(long i=1; i<PARS.Ny; i++)
    {
      y_1 = yt_arr(i-1);      
      y_2 = yt_arr(i);
      dy = (y_2 - y_1)*pre_df;  // adaptive dy
      y_h = y_1 + dy/2.;
      phi_1 = phi_yt_div_phib_arr(index_z, i-1);
      double shear_rate_nodim_1 = shear_rate_nodim_pre_ZY(index_z, i-1);
      double shear_rate_nodim_2 = shear_rate_nodim_pre_ZY(index_z, i);
      double shear_rate_nodim_1_half = (shear_rate_nodim_1 + shear_rate_nodim_2)/2.; // simple half expectation
      
      k1 = dy * get_f_div_fb_modi_RK(y_1, 0., phi_1, 0.,	int_INV_D_pre, vw_div_vw0, shear_rate_nodim_1, shear_rate_nodim_1, DISP, PARS);
      k2 = dy * get_f_div_fb_modi_RK(y_1, dy/2., phi_1, k1/2., int_INV_D_pre, vw_div_vw0, shear_rate_nodim_1, shear_rate_nodim_1_half, DISP, PARS);
      k3 = dy * get_f_div_fb_modi_RK(y_1, dy/2., phi_1, k2/2., int_INV_D_pre, vw_div_vw0, shear_rate_nodim_1, shear_rate_nodim_1_half, DISP, PARS);
      k4 = dy * get_f_div_fb_modi_RK(y_1, dy, phi_1, k3,	int_INV_D_pre, vw_div_vw0, shear_rate_nodim_1, shear_rate_nodim_2, DISP, PARS);
      phi_2 = phi_1 + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4);
      phi_yt_div_phib_arr(index_z, i) = phi_2;

	
      int_INV_D_pre += (dy/2.) * (1./DISP.fcn_Dc_div_D0(phi_2*PARS.phi_bulk, PARS, shear_rate_nodim_pre_ZY(index_z, i)) + 1./DISP.fcn_Dc_div_D0(phi_1*PARS.phi_bulk, PARS, shear_rate_nodim_pre_ZY(index_z, i-1)));

      /*
	Due to the precision problem, we made phi_yt array into phi_yt_div_phib array.
	Therefore, the calculation of collective diffusion coefficient must take re-mapping with phi_bulk.
       */
    }
  return 0;
}

int get_Ieta_yt(long index_z, MATRIX& Ieta_yt_arr, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& shear_stress_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  /*
    NOTE: shear_rate and viscosity are coupled.
    Therefore, as the equilibrium transport properties, we calculate viscosity for the given concentration (previous step) and shear rate (previous step).
    Note that shear rate calculation used for the given viscosity in that context, independently.
   */
  double re = 0.;
  Ieta_yt_arr(index_z, 0) = re;
  double tmp_f1 = 1./DISP.fcn_eta_div_eta0(phi_yt_div_phib_arr(index_z, 0)*PARS.phi_bulk, PARS, shear_stress_nodim_pre_ZY(index_z, 0));
  double tmp_f2 = tmp_f1;
  double pre_df = 1.;
  // pre_df = 1. - delta_f_div_R
  double dy = 0.;
  for (long i = 1; i<PARS.Ny; i++)
    {
      dy = (yt_arr(i) - yt_arr(i-1))*pre_df;
      tmp_f1 = tmp_f2;
      tmp_f2 = 1./DISP.fcn_eta_div_eta0(phi_yt_div_phib_arr(index_z, i)*PARS.phi_bulk, PARS, shear_stress_nodim_pre_ZY(index_z, i));
      
      Ieta_yt_arr(index_z, i) = Ieta_yt_arr(index_z, i-1) + 0.5*dy*(tmp_f1 + tmp_f2);
    }


  for (long i = 1; i<PARS.Ny; i++)
    {
      Ieta_yt_arr(index_z, i) /= Ieta_yt_arr(index_z, PARS.Ny - 1); // normalize to 1 at centerline
    }
  
  // if (index_z == 100)
  //   std::printf("Ieta(%4.3e) = %4.3e\n",  yt_arr(PARS.Ny-1), Ieta_yt_arr(index_z, PARS.Ny-1));
  return 0;
}

int get_ID_yt(long index_z, MATRIX& ID_yt_arr, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& shear_rate_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  double re = 0.;
  ID_yt_arr(index_z, 0) = re;
  double tmp_f1 = 1./DISP.fcn_Dc_div_D0(phi_yt_div_phib_arr(index_z, 0)*PARS.phi_bulk, PARS, shear_rate_nodim_pre_ZY(index_z, 0));
  double tmp_f2 = tmp_f1;
  double pre_df = 1.;
  // pre_df = 1. - delta_f_div_R
  double dy = 0.;
  for (long i = 1; i<PARS.Ny; i++)
    {
      dy = (yt_arr(i) - yt_arr(i-1))*pre_df;
      tmp_f1 = tmp_f2;
      tmp_f2 = 1./DISP.fcn_Dc_div_D0(phi_yt_div_phib_arr(index_z, i)*PARS.phi_bulk, PARS, shear_rate_nodim_pre_ZY(index_z, i));
      ID_yt_arr(index_z, i) = ID_yt_arr(index_z, i-1) + 0.5*dy*(tmp_f1 + tmp_f2);
    }
  return 0;
}

double get_F1(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS)
{
  double re = 0.;
  double tmp_f1 = (1. - yt_arr(0))*(2. - yt_arr(0))*Ieta_yt_arr(index_z, 0)*exp(-(vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, 0));
  double tmp_f2 = tmp_f1;
  double pre_df = 1.;
  double dy = 0.;
  for (long i=1; i<PARS.Ny; i++)
    {
      dy = (yt_arr(i) - yt_arr(i-1))*pre_df;
      tmp_f1 = tmp_f2;
      tmp_f2 = (1. - yt_arr(i))*(2. - yt_arr(i))*Ieta_yt_arr(index_z, i)*exp(-(vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, i));
      
      re += 0.5 * dy * (tmp_f1 + tmp_f2);
      
    }
  return re*uz_div_u0;
}



double get_F2(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS)
{
  // vw_div_vw0 is not necessary
  double re = 0.;
  double tmp_f1 = (1. - yt_arr(0))*(2. - yt_arr(0))*Ieta_yt_arr(index_z, 0);
  double tmp_f2 = tmp_f1;
  double pre_df = 1.;
  // double pre_df = 1. - delta_f_div_R;

  double dy = 0.;
  for (long i=1; i<PARS.Ny; i++)
    {
      dy = (yt_arr(i) - yt_arr(i-1))*pre_df;
      tmp_f1 = tmp_f2;
      tmp_f2 = (1. - yt_arr(i))*(2. - yt_arr(i))*Ieta_yt_arr(index_z, i);
      re += 0.5 * dy * (tmp_f1 + tmp_f2);
    }
  return re*uz_div_u0;
}

double get_F_modi(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS)
{
  double re = 0.;
  double tmp_f1 = (1. - yt_arr(0))*(2. - yt_arr(0))*Ieta_yt_arr(index_z, 0)*exp(-(vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, 0))*((vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, 0));
  double tmp_f2 = tmp_f1;
  double pre_df = 1.;

  double dy = 0.;
  for (long i=1; i<PARS.Ny; i++)
    {
      dy = (yt_arr(i) - yt_arr(i-1))*pre_df;
      tmp_f1 = tmp_f2;
      tmp_f2 = (1. - yt_arr(i))*(2. - yt_arr(i))*Ieta_yt_arr(index_z, i)*exp(-(vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, i))*((vw_div_vw0/PARS.epsilon_d)*ID_yt_arr(index_z, i));
      re += 0.5 * dy * (tmp_f1 + tmp_f2);
    }
  return re*uz_div_u0;
}

// At phi_freeze, search resistance against the inherited flux functional.
// This thresholded cake closure has no counterpart in the public Python code.
double cal_Rcake_div_Rmem(long index_z, long n_iter, MATRIX& yt_arr, double F2_0, double F2_0_modi, double uz_div_u0, double vw_div_vw0, double P_div_DLP, double Rcake_div_Rmem, double phiw_div_phib, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_cake)
{
  /*
    About shear-induced migration: the information of shear rate and stress are provided by individual arguments (originated from AUX_h1).
    The information from AUX_h2 will not be used for those shear rate and stress since it is only important for the aspect of concentration, Ieta_yt_arr, and ID_yt_arr (using the individually provided information of shear rate and stress)
   */

  double re = PARS.phi_freeze/PARS.phi_bulk;

  get_phi_yt_div_phib_with_given_phiw(index_z, AUX_cake.phi_yt_div_phib_arr, yt_arr, phiw_div_phib, vw_div_vw0, shear_rate_nodim_pre_ZY, DISP, PARS);
  get_Ieta_yt(index_z, AUX_cake.Ieta_yt_arr, yt_arr, AUX_cake.phi_yt_div_phib_arr, shear_stress_nodim_pre_ZY, DISP, PARS);
  get_ID_yt(index_z, AUX_cake.ID_yt_arr, yt_arr, AUX_cake.phi_yt_div_phib_arr, shear_rate_nodim_pre_ZY, DISP, PARS);


  double tmp_nodim_R_0 = 0.;
  double tmp_nodim_R_1 = 0.;
  double tmp_nodim_R_2 = tmp_nodim_R_1;

  // double DEL_nodim_R_tmp = 0.01*PARS.CAKE_convert_unit; // WIP
  double DEL_nodim_R_tmp = 0.01;
  // if(PARS.CAKE_convert_unit > 0.)
  //   DEL_nodim_R_tmp = 0.01*PARS.CAKE_convert_unit; // WIP
  

  double f_0 = PARS.phi_bulk*(F2_0 - F2_0_modi);
  double Pi_div_DLP = DISP.fcn_Pi_div_Piast(PARS.phi_freeze, PARS, shear_rate_nodim_pre_ZY(0, 0))*PARS.Pi_ast/PARS.DLP;
  long tmp_cnt = 0;
  double f_i = -1.;
  double f_i_pre = f_i;
  /*
    comments on f_i = -1
    this is for the first identification with large step size gives positive value
    if the initial step size is h0 and the calculated f_i(h0) > 0,
    then the sign identifier f_i*f_i_pre must be negative
    in order to reduce step size from h0 to h0/10 (let say).
  */
  // double delta_f_div_R_new = 0.;

  double tmp_F1 = 0., tmp_F2_Z = 0., tmp_F2_Z_modi = 0.;
  double effective_nodim_R_range = 1e2;

  double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem;
  double pre_df = pow(1. - cake_thickness_nodim, 2.0);
  // double pre_df = 1.; 
  /*
    This located outside since each time of Rcake update will change pre_df value.
    This is non-logical since we update pre_df value from the previous iteration.
   */
  // This is a numerical search guard, not a different cake closure. A valid
  // search follows the inherited steps unchanged; stalled/NaN cases fail clearly.
  long search_steps = 0;
  while(1)
    {
      if (++search_steps > 200000) throw std::runtime_error("cake resistance search did not terminate");

      tmp_nodim_R_1 = tmp_nodim_R_2;

      // double pre_df = 1.0;
      f_i_pre = f_i;
      tmp_nodim_R_2 += DEL_nodim_R_tmp;
      double tmp_vw_div_vw0 = PARS.Lp*get_nodim_Lpf(tmp_nodim_R_2)*(P_div_DLP*PARS.DLP - PARS.Pper - Pi_div_DLP*PARS.DLP)/PARS.vw0;
      
      tmp_F1 = pre_df*get_F1(index_z, uz_div_u0, tmp_vw_div_vw0, yt_arr, AUX_cake.phi_yt_div_phib_arr, AUX_cake.Ieta_yt_arr, AUX_cake.ID_yt_arr, PARS);
      tmp_F2_Z = pre_df*get_F2(index_z, uz_div_u0, tmp_vw_div_vw0, yt_arr, AUX_cake.phi_yt_div_phib_arr, AUX_cake.Ieta_yt_arr, AUX_cake.ID_yt_arr, PARS);
      tmp_F2_Z_modi = pre_df*get_F_modi(index_z, uz_div_u0, tmp_vw_div_vw0, yt_arr, AUX_cake.phi_yt_div_phib_arr, AUX_cake.Ieta_yt_arr, AUX_cake.ID_yt_arr, PARS);

      f_i = (PARS.phi_freeze - PARS.phi_bulk)*tmp_F1 + PARS.phi_bulk*(tmp_F2_Z - tmp_F2_Z_modi) - f_0;
      // printf("no_predf: tmp_F1=%4.3e, tmp_F2_Z = %4.3e, tmp_F2_Z_modi = %4.3e, f_0=%4.3e\n", tmp_F1/pre_df, tmp_F2_Z/pre_df, tmp_F2_Z_modi/pre_df, f_0);

      /*
	At this point, I guess instability is due to the increased cake thickness does not contribute directly on the axial velocity profile. What if the velocity field is simultaneously updated?
       */
      
      f_i /= f_0;
      if (!std::isfinite(f_i)) throw std::runtime_error("nonfinite cake residual");


      if (tmp_nodim_R_2> effective_nodim_R_range || (f_i != -1 && f_i - f_i_pre < 0))
	{
	  if (DEL_nodim_R_tmp <=PARS.TOL_delta_f)
	    {
	      printf("Warning: f(%04ld) - f_pre(%04ld) = %4.3e < 0\n",index_z, index_z, f_i - f_i_pre);
	      printf("DEL_nodim_R_tmp=%4.3e, PARS.TOL_delta_f = %4.3e, tmp_nodim_R_2 = %4.3e\n", DEL_nodim_R_tmp, PARS.TOL_delta_f, tmp_nodim_R_2);
      printf("pre_df = %4.3e, tmp_F1=%4.3e, tmp_F2_Z = %4.3e, tmp_F2_Z_modi = %4.3e, f_0=%4.3e\n", pre_df, tmp_F1, tmp_F2_Z, tmp_F2_Z_modi, f_0);
      printf("uz_div_u0 = %4.3e\n", uz_div_u0);
      
	      tmp_nodim_R_2 = Rcake_div_Rmem;//delta_f_div_R/PARS.epsilon_f; // CHECK
	      break;
	    }
	  
	  // this is sometimes the initial guess is too large
	  DEL_nodim_R_tmp /= 10.;
	  tmp_nodim_R_2 = tmp_nodim_R_0;
	  f_i = f_i_pre;
	  tmp_cnt = 0; // it reset tmp_cnt which is the identifier for exceeding numbers
	}
      
      if (f_i * f_i_pre < 0 )
	{
	  if (DEL_nodim_R_tmp <= PARS.TOL_delta_f)
	    {
	      break;
	    }
	  else
	    {
	      DEL_nodim_R_tmp /= 10.;
	      tmp_nodim_R_2 = tmp_nodim_R_1;
	      f_i = f_i_pre;
	      tmp_cnt = 0; // it reset tmp_cnt which is the identifier for exceeding numbers
	    }
	}
      tmp_cnt += 1;
    }
  Rcake_div_Rmem = tmp_nodim_R_1;
  return Rcake_div_Rmem;
}


double process_at_z_modi(long index_z, MATRIX& yt_arr, double phiw_div_phib, double phiw_div_phib_pre, double F2_0, double F2_0_modi, double uz_div_u0, double vw_div_vw0, double Rcake_div_Rmem, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_h2)
{
  /*
    About shear-induced migration: the information of shear rate and stress are provided by individual arguments (originated from AUX_h1).
    The information from AUX_h2 will not be used for those shear rate and stress since it is only important for the aspect of concentration, Ieta_yt_arr, and ID_yt_arr (using the individually provided information of shear rate and stress)
   */
  get_phi_yt_div_phib_with_given_phiw(index_z, AUX_h2.phi_yt_div_phib_arr, yt_arr, phiw_div_phib, vw_div_vw0, shear_rate_nodim_pre_ZY, DISP, PARS);

  get_Ieta_yt(index_z, AUX_h2.Ieta_yt_arr, yt_arr, AUX_h2.phi_yt_div_phib_arr, shear_stress_nodim_pre_ZY, DISP, PARS);
  get_ID_yt(index_z, AUX_h2.ID_yt_arr, yt_arr, AUX_h2.phi_yt_div_phib_arr, shear_rate_nodim_pre_ZY, DISP, PARS);

  // double pre_df = 1.; // 1. - delta_f_div_R
  double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem;
  double pre_df = pow(1. - cake_thickness_nodim, 2.0);
  // double pre_df = 1.;
  // this account the reduced cross-sectional area
  // therefore, the assumption of thickness is of importance
  // the application is only for the cross-sectional averaged value with z>0
    
  double F1 = pre_df*get_F1(index_z, uz_div_u0, vw_div_vw0, yt_arr, AUX_h2.phi_yt_div_phib_arr, AUX_h2.Ieta_yt_arr, AUX_h2.ID_yt_arr, PARS);
  
  double F2_Z = pre_df*get_F2(index_z, uz_div_u0, vw_div_vw0, yt_arr, AUX_h2.phi_yt_div_phib_arr, AUX_h2.Ieta_yt_arr, AUX_h2.ID_yt_arr, PARS);
  double F2_Z_modi = pre_df*get_F_modi(index_z, uz_div_u0, vw_div_vw0, yt_arr, AUX_h2.phi_yt_div_phib_arr, AUX_h2.Ieta_yt_arr, AUX_h2.ID_yt_arr, PARS);

  // CHECK
  AUX_h2.F_info_arr(index_z, 1) = F1;
  AUX_h2.F_info_arr(index_z, 2) = F2_Z;
  AUX_h2.F_info_arr(index_z, 3) = F2_Z_modi;  
  
  double re = (1. + (F2_0 - F2_0_modi - (F2_Z - F2_Z_modi))/F1);
  return re;
}


long phiw_update(MATRIX& phiw_set_2, MATRIX& Rcake_div_Rmem_new_arr, MATRIX& z_div_L_arr, MATRIX& phiw_set_1, long ind_onset_cake_pre, MATRIX& Rcake_div_Rmem_arr, MATRIX& Pi_div_DLP_arr, MATRIX& yt_arr, long n_iter, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_h1, AUX_MATRIX& AUX_h2)
{

  /*
    The usage of AUX_h1, AUX_h2, and AUX_h3:
    AUX_h1 is the provided condition (previous iteration step)
    AUX_h2 is the typical CP-layer computation and cake layer (present iteration step)

    For the future, this function must be revised accordingly.
   */
  
  z_div_L_arr(0) = 0.;
  for (long i=1; i<PARS.Nz; i++)
    {
      z_div_L_arr(i) = z_div_L_arr(i-1) + PARS.dz/PARS.L_channel;
    }


  MATRIX& shear_rate_nodim_pre_ZY = AUX_h1.shear_rate_nodim_pre; // avoiding confusion. We use AUX_h1 for the shear rate.
  MATRIX& shear_stress_nodim_pre_ZY = AUX_h1.shear_stress_nodim_pre; // avoiding confusion. We use AUX_h1 for the shear stress.

  get_nodim_Lpf_arr(AUX_h1.Ln_arr, Rcake_div_Rmem_arr, PARS);
  
  get_P_div_DLP_arr_with_fouling(AUX_h1.P_div_DLP_arr, PARS, AUX_h1.Ln_arr, Pi_div_DLP_arr, AUX_h1.A_mat, AUX_h1.y_mat, AUX_h1.d_mat, AUX_h1.R_arr, AUX_h1.Q_arr, AUX_h1.f_arr, Rcake_div_Rmem_arr);
  get_uz_div_u0_arr_with_fouling(AUX_h1.uz_div_u0_arr, z_div_L_arr, AUX_h1.P_div_DLP_arr, Rcake_div_Rmem_arr, PARS);
  get_vw_div_vw0_arr_with_fouling(AUX_h1.vw_div_vw0_arr, AUX_h1.P_div_DLP_arr, Pi_div_DLP_arr, AUX_h1.Ln_arr, Rcake_div_Rmem_arr, PARS);

  double uz_div_u0_z0 = AUX_h1.uz_div_u0_arr(0);
  double vw_div_vw0_z0 = AUX_h1.vw_div_vw0_arr(0);

  get_phi_yt_div_phib_with_given_phiw(0, AUX_h1.phi_yt_div_phib_arr, yt_arr, phiw_set_1(0), vw_div_vw0_z0, shear_rate_nodim_pre_ZY, DISP, PARS);

  get_Ieta_yt(0, AUX_h1.Ieta_yt_arr, yt_arr, AUX_h1.phi_yt_div_phib_arr, shear_stress_nodim_pre_ZY, DISP, PARS);
  get_ID_yt(0, AUX_h1.ID_yt_arr, yt_arr, AUX_h1.phi_yt_div_phib_arr, shear_rate_nodim_pre_ZY, DISP, PARS);


  double F1_0 = get_F1(0, uz_div_u0_z0, vw_div_vw0_z0, yt_arr, AUX_h1.phi_yt_div_phib_arr, AUX_h1.Ieta_yt_arr, AUX_h1.ID_yt_arr, PARS);
  
  double F2_0 = get_F2(0, uz_div_u0_z0, vw_div_vw0_z0, yt_arr, AUX_h1.phi_yt_div_phib_arr, AUX_h1.Ieta_yt_arr, AUX_h1.ID_yt_arr, PARS);
  double F2_0_modi = get_F_modi(0, uz_div_u0_z0, vw_div_vw0_z0, yt_arr, AUX_h1.phi_yt_div_phib_arr, AUX_h1.Ieta_yt_arr, AUX_h1.ID_yt_arr, PARS);

  // this is just a record purpose
  AUX_h2.F_info_arr(0, 1) = F1_0;
  AUX_h2.F_info_arr(0, 2) = F2_0;
  AUX_h2.F_info_arr(0, 3) = F2_0_modi;  
  
  double tmp_Rcake_div_Rmem = 0.;

  long ind_onset_cake = PARS.Nz + 1;

  double tmp_phiw = 0.;
  double phi_f_div_phi_b = PARS.phi_freeze/PARS.phi_bulk;

  bool safety_rule = FALSE;
  // bool discontinuity = FALSE;
  double dz_nodim = PARS.dz/PARS.L_channel; // dimensionless step size  
  double TOLERANCE_DISCONTINUOUS = 10.;

    
  for (long i=1; i<PARS.Nz; i++)
    {

      phiw_set_2(i) = process_at_z_modi(i, yt_arr, phiw_set_1(i), phiw_set_2(i), F2_0, F2_0_modi, AUX_h1.uz_div_u0_arr(i), AUX_h1.vw_div_vw0_arr(i), Rcake_div_Rmem_arr(i), shear_rate_nodim_pre_ZY, shear_stress_nodim_pre_ZY, DISP, PARS, AUX_h2);
      
      AUX_h1.F_info_arr.copy_obj(AUX_h2.F_info_arr); // CHECK

      
      safety_rule = FALSE;
      // discontinuity = FALSE;
      // Preserve the inherited history-dependent cake-persistence guard:
      // previous cake at i and current upstream cake prevent a raw candidate
      // from dropping below freeze here. This precedes relaxation/clamping;
      // it is not a general monotonicity validator. See FIXED_POINT_METHOD.org.
      if (phiw_set_2(i) < phi_f_div_phi_b && phiw_set_1(i) >= phi_f_div_phi_b && phiw_set_2(i-1) >= phi_f_div_phi_b)
	{
	  // this will prevent discontinuity of cake formation
	  // possibly due to numerical error induces just 'slightly' small values of wall concentration
	  // which identified as NO CAKE layer formed (although the previous z_i formed the cake layer)
	  safety_rule = TRUE;
	}
	
      if (safety_rule)
	
	{
	  phiw_set_2(i) = phiw_set_2(i-1);
	  
	  PARS.count_negative_phiw ++; // this will not be important for the iteration in the beginning, but it is of important to when it influence strongly on the cake layer formation, and its unstable prediction.
	}
      tmp_phiw = (1. - PARS.weight)*phiw_set_1(i) + PARS.weight*phiw_set_2(i);
      if(tmp_phiw < 0.)
	{
	  // cout << "WARNING: negative phiw = " << phiw_set_2(i)*PARS.phi_bulk << endl;
	  tmp_phiw = phiw_set_1(i);
	}
      
      if (tmp_phiw > phi_f_div_phi_b && ind_onset_cake > PARS.Nz)
	{
	  // when it is onset of cake formation
	  ind_onset_cake = i;

	}

      // if (ind_onset_cake <= PARS.Nz) // MONOTONIC PHIW
      if(tmp_phiw >= phi_f_div_phi_b)
	{

	  phiw_set_2(i) = phi_f_div_phi_b;

	  tmp_Rcake_div_Rmem = cal_Rcake_div_Rmem(i, n_iter, yt_arr, F2_0, F2_0_modi, AUX_h1.uz_div_u0_arr(i), AUX_h1.vw_div_vw0_arr(i), AUX_h1.P_div_DLP_arr(i), Rcake_div_Rmem_arr(i), phiw_set_2(i), shear_rate_nodim_pre_ZY, shear_stress_nodim_pre_ZY, DISP, PARS, AUX_h2);
	  
	  Rcake_div_Rmem_new_arr(i) = (1. - PARS.weight_f)*Rcake_div_Rmem_arr(i) + PARS.weight_f*tmp_Rcake_div_Rmem;

	}
      else
	{
	  // when there is no cake 
	  phiw_set_2(i) = tmp_phiw;
	}
    }
  long ind_onset_cake_min = check_min(ind_onset_cake, ind_onset_cake_pre);
  // // MONOTONIC PHIW


  update_shear_rate_nodim_ZY(AUX_h1.shear_rate_nodim, shear_rate_nodim_pre_ZY, shear_stress_nodim_pre_ZY, AUX_h2.Ieta_yt_arr, yt_arr, AUX_h2.phi_yt_div_phib_arr, AUX_h1.uz_div_u0_arr, DISP, PARS);
  update_shear_stress_nodim_ZY(AUX_h1.shear_stress_nodim, shear_stress_nodim_pre_ZY, AUX_h1.shear_rate_nodim, AUX_h2.phi_yt_div_phib_arr, DISP, PARS);      
  
  
  return ind_onset_cake;
}


int update_shear_rate_nodim_ZY(MATRIX& shear_rate_nodim_ZY, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, MATRIX& Ieta_yt_arr_ZY, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& uz_div_u0_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  for(long index_z=1; index_z < PARS.Nz; index_z++) // start with 1 because inlet is without CP layer
    {
      for(long i=0; i<PARS.Ny-1; i++)
	{
	  double phi = phi_yt_div_phib_arr(index_z, i)*PARS.phi_bulk;
	  double eta_phi = DISP.fcn_eta_div_eta0(phi, PARS, shear_stress_nodim_pre_ZY(index_z, i));
	  shear_rate_nodim_ZY(index_z, i) = (-Ieta_yt_arr_ZY(index_z, i) + (2. - yt_arr(i))/eta_phi)*uz_div_u0_arr(index_z);
	  if(shear_rate_nodim_ZY(index_z, i) < 0)
	    {
	      // cout << "Warning shear_rate is negative. Taking absolute value of it for safety reason: " << shear_rate_nodim_ZY(index_z, i)<< endl;
	      shear_rate_nodim_ZY(index_z, i) *= -1.;
	    }
	  shear_rate_nodim_ZY(index_z, i) *= (PARS.Pe_a/2.);
	}


      shear_rate_nodim_ZY(index_z, PARS.Ny-1) = 0.; // this will be the center
    }

  return 0;
}


int update_shear_stress_nodim_ZY(MATRIX& shear_stress_nodim_ZY, MATRIX& shear_stress_nodim_pre_ZY, MATRIX& shear_rate_nodim_ZY, MATRIX& phi_yt_div_phib_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS)
{
  for(long index_z=1; index_z < PARS.Nz; index_z++) // start with 1 because inlet is without CP layer
    {
      for(long i=0; i<PARS.Ny-1; i++)
	{
	  double phi = phi_yt_div_phib_arr(index_z, i)*PARS.phi_bulk;
	  double eta_phi = DISP.fcn_eta_div_eta0(phi, PARS, shear_stress_nodim_pre_ZY(index_z, i));
	  shear_stress_nodim_ZY(index_z, i) = shear_rate_nodim_ZY(index_z, i) * eta_phi/(6.*M_PI);
	}


      shear_rate_nodim_ZY(index_z, PARS.Ny-1) = 0.; // this will be the center
    }

  return 0;
}


int
shear_rate_initialization
(MATRIX& shear_rate_nodim_ZY, MATRIX& yt_arr, PARAMETER_SET& PARS)
{
  /* 
     The initializer uses the pure Hagen-Poiseiulle flow profile.
     This means that the axial velocity along center-lines is constant, and the shear-rate at the wall does not change along the z-direction.
     The procedure, of course, the over-estimation of shear rate since the viscosity effect is not accounted, as the CP-layer is not accounted in the begining of FPI. 
     Note that it is a benefit to have over-estimation of shear-rate in the beginning since the it usually in the way to supress the CP-layer, and, therefore, the possible diverging tendency in FPI is suppressed.

     Note that the definition of characteristic shear rate is 2*u^HP/R for inlet-pressure boundary condition.
   */
  for (long j=0; j<PARS.Ny; j++)
    {
      for (long i=0; i<PARS.Nz; i++)
	{
	  // WIP: it must be dimensionless quantity. 
	  shear_rate_nodim_ZY(i,j) = (1. - yt_arr(j));

	}
    }
  return 0;
}

int
shear_stress_initialization
(MATRIX& shear_stress_nodim_ZY, MATRIX& shear_rate_nodim_ZY, MATRIX& yt_arr, PARAMETER_SET& PARS, DISP_TYPE& DISP)
{
  /*
    Initial condition for shear stress: tau^0 = eta(phi=phi_b, tau=0)*dg^0

   */
  double pre_initial_shear_stress_nodim_ZY = 0.;
  for (long j=0; j<PARS.Ny; j++)
    {
      for (long i=0; i<PARS.Nz; i++)
	{
	  // the prefactor 1/(6*pi) is due to dimensionless tau:
	  // dimensionless shear stress = relative viscosity * dimensionless shear rate / (6 * pi)
	  // since base unit for shear rate is D0/a^2
	  // and base unit for shear stress is kT/a^3
	  // shear_stress_nodim_ZY(i,j) = (1./(6.*M_PI))*shear_rate_nodim_ZY(i,j)*DISP.fcn_eta_div_eta0(PARS.phi_bulk, PARS, pre_initial_shear_stress_nodim_ZY);
	  shear_stress_nodim_ZY(i,j) = (1./(6.*M_PI))*shear_rate_nodim_ZY(i,j)*DISP.fcn_eta_div_eta0(PARS.phi_bulk, PARS, pre_initial_shear_stress_nodim_ZY);
	  
	}
    }
  return 0;
}




int fprint_analysis_with_fouling(MATRIX& z_div_L_arr, MATRIX& phiw_div_phib_arr, MATRIX& Rcake_div_Rmem_arr, MATRIX& yt_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS, string fn_out, AUX_MATRIX& AUX)
{

    
  MATRIX Pi_div_DLP_arr(PARS.Nz, 1, 0.);
  update_Pi_div_DLP_arr(Pi_div_DLP_arr, phiw_div_phib_arr, AUX.shear_rate_nodim_pre, DISP, PARS);

  get_nodim_Lpf_arr(AUX.Ln_arr, Rcake_div_Rmem_arr, PARS);
  
  get_P_div_DLP_arr_with_fouling(AUX.P_div_DLP_arr, PARS, AUX.Ln_arr, Pi_div_DLP_arr, AUX.A_mat, AUX.y_mat, AUX.d_mat, AUX.R_arr, AUX.Q_arr, AUX.f_arr, Rcake_div_Rmem_arr);
  get_uz_div_u0_arr_with_fouling(AUX.uz_div_u0_arr, z_div_L_arr, AUX.P_div_DLP_arr, Rcake_div_Rmem_arr, PARS);
  get_vw_div_vw0_arr_with_fouling(AUX.vw_div_vw0_arr, AUX.P_div_DLP_arr, Pi_div_DLP_arr, AUX.Ln_arr, Rcake_div_Rmem_arr, PARS);

  MATRIX re(PARS.Nz, 10, 0.);
  for(long i=0; i<PARS.Nz; i++)
    {
      

      // calculation cross-sectional particle-flux conservation
      get_phi_yt_div_phib_with_given_phiw(i, AUX.phi_yt_div_phib_arr, yt_arr, phiw_div_phib_arr(i), AUX.vw_div_vw0_arr(i), AUX.shear_rate_nodim, DISP, PARS);
      get_Ieta_yt(i, AUX.Ieta_yt_arr, yt_arr, AUX.phi_yt_div_phib_arr, AUX.shear_stress_nodim, DISP, PARS);
      get_ID_yt(i, AUX.ID_yt_arr, yt_arr, AUX.phi_yt_div_phib_arr, AUX.shear_rate_nodim, DISP, PARS);

      // double pre_df = 1.; // 1. - delta_f_div_R
      double cake_thickness_nodim = PARS.CAKE_convert_unit * Rcake_div_Rmem_arr(i);
      double pre_df = pow(1. - cake_thickness_nodim, 2.0);
      // double pre_df = 1.;
      double F1 = pre_df*get_F1(i, AUX.uz_div_u0_arr(i), AUX.vw_div_vw0_arr(i), yt_arr, AUX.phi_yt_div_phib_arr, AUX.Ieta_yt_arr, AUX.ID_yt_arr, PARS);
      double F2_Z = pre_df*get_F2(i, AUX.uz_div_u0_arr(i), AUX.vw_div_vw0_arr(i), yt_arr, AUX.phi_yt_div_phib_arr, AUX.Ieta_yt_arr, AUX.ID_yt_arr, PARS);
      double F2_Z_modi = pre_df*get_F_modi(i, AUX.uz_div_u0_arr(i), AUX.vw_div_vw0_arr(i), yt_arr, AUX.phi_yt_div_phib_arr, AUX.Ieta_yt_arr, AUX.ID_yt_arr, PARS);
      double Phi_z_div_Phi_ast = ((phiw_div_phib_arr(i) - 1.)*F1 + (F2_Z - F2_Z_modi));

      re(i, 0) = z_div_L_arr(i)		*PARS.L_channel;
      re(i, 1) = phiw_div_phib_arr(i)	*PARS.phi_bulk;
      re(i, 2) = AUX.P_div_DLP_arr(i)	*PARS.DLP;
      re(i, 3) = AUX.vw_div_vw0_arr(i)	*PARS.vw0;
      re(i, 4) = AUX.uz_div_u0_arr(i)	*PARS.u0;
      re(i, 5) = Pi_div_DLP_arr(i)	*PARS.DLP;
      re(i, 6) = Phi_z_div_Phi_ast	*2. * M_PI * pow(PARS.R_channel, 2.0) * PARS.phi_bulk * PARS.u0;
      re(i, 7) = Rcake_div_Rmem_arr(i)	*PARS.R_channel;
      re(i, 8) = AUX.P_div_DLP_arr(i)	*PARS.DLP - PARS.Pper;
      re(i, 9) = cake_thickness_nodim;

    }

  ofstream fout(fn_out, ofstream::out);
  if (!fout) throw std::runtime_error("cannot create output table: " + fn_out);
  for (long j=0; j<re.size; ++j)
    if (!std::isfinite(re(j))) throw std::runtime_error("nonfinite final output; refine/review the numerical case");  
  re.fprint_basic(fout);
  finish_output(fout, fn_out);
  return 0;
}




