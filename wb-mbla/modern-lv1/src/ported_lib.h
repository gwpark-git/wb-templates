#ifndef PORTED_LIB_H
#define PORTED_LIB_H


#include "parameters.h"
#include "dispersion_properties.h"
#include "matrix.h"
#include "boost_matrix.h"
extern "C"{
#include <gsl/gsl_const_mksa.h> // physical constants
#include <gsl/gsl_math.h> // mathematical constants
}

#define TRUE 1
#define FALSE 0

int fprint_analysis_with_fouling(MATRIX& z_div_L_arr, MATRIX& phiw_div_phib_arr, MATRIX& Rcake_div_Rmem_arr, MATRIX& yt_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS, string fn_out, AUX_MATRIX& AUX);

int gen_z_div_L_arr(MATRIX& z_div_L_arr, PARAMETER_SET& PARS);

int gen_yt_arr(MATRIX& yt_arr, PARAMETER_SET& PARS);

int update_Pi_div_DLP_arr(MATRIX& Pi_div_DLP_arr, MATRIX& phiw_set, MATRIX& shear_rate_nodim_NOT_USED, DISP_TYPE& DISP, PARAMETER_SET& PARS);

double REL_ERR_data(MATRIX& data_set_2, MATRIX& data_set_1, long skip_rows=0);

double ABS_ERR_data(MATRIX& data_set_2, MATRIX& data_set_1, long skip_rows=0);

int get_df_from_f(MATRIX& df_arr, MATRIX& z_div_L_arr, MATRIX& f_arr);

double get_nodim_Lpf(double f);

int get_nodim_Lpf_arr(MATRIX& Ln_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS);

int init_get_P(MATRIX& A_mat, MATRIX& y_mat, MATRIX &d_mat, MATRIX& R_arr, MATRIX& Q_arr, MATRIX& f_arr, PARAMETER_SET& PARS);

int get_P_div_DLP_arr_with_fouling(MATRIX& P_div_DLP_arr, PARAMETER_SET& PARS, MATRIX& Ln_arr, MATRIX& Pi_div_DLP_arr, MATRIX& A_mat, MATRIX& y_mat, MATRIX &d_mat, MATRIX& R_arr, MATRIX& Q_arr, MATRIX& f_arr, MATRIX& Rcake_div_Rmem_arr);

int get_uz_div_u0_arr_with_fouling(MATRIX& uz_div_u0_arr, MATRIX& z_div_L_arr, MATRIX& P_div_DLP_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS);

int get_vw_div_vw0_arr_with_fouling(MATRIX& vw_div_vw0_arr, MATRIX& P_div_DLP_arr, MATRIX& Pi_div_DLP_arrr, MATRIX& Ln_arr, MATRIX& Rcake_div_Rmem_arr, PARAMETER_SET& PARS);

double get_f_div_fb_modi_RK(double yt, double dyt, double f, double df, double int_INV_D_pre, double vw_div_vw0, double dg_1, double dg_2, DISP_TYPE& DISP, PARAMETER_SET& PARS);

int get_phi_yt_div_phib_with_given_phiw(long index_z, MATRIX& phi_yt_div_phib_arr, MATRIX& yt_arr, double given_phiw_div_phib, double vw_div_vw0, MATRIX& shear_rate_nodim_pre_ZY, DISP_TYPE& DISP,  PARAMETER_SET& PARS);

int get_Ieta_yt(long index_z, MATRIX& Ieta_yt_arr, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& shear_rate_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS);

int get_ID_yt(long index_z, MATRIX& ID_yt_arr, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& shear_rate_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS);

double get_F1(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS);

double get_F2(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS);

double get_F_modi(long index_z, double uz_div_u0, double vw_div_vw0, MATRIX& yt_arr, MATRIX& phi_yt_arr, MATRIX& Ieta_yt_arr, MATRIX& ID_yt_arr, PARAMETER_SET& PARS);

double cal_Rcake_div_Rmem(long index_z, long n_iter, MATRIX& yt_arr, double F2_0, double F2_0_modi, double uz_div_u0, double vw_div_vw0, double P_div_DLP, double Rcake_div_Rmem, double phiw_div_phib, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_cake);

double process_at_z_modi(long index_z, MATRIX& yt_arr, double phiw_div_phib, double phiw_div_phib_pre, double F2_0, double F2_0_modi, double uz_div_u0, double vw_div_vw0, double Rcake_div_Rmem, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_h2);

long phiw_update(MATRIX& phiw_set_2, MATRIX& Rcake_div_Rmem_new_arr, MATRIX& z_div_L_arr, MATRIX& phiw_set_1, long ind_onset_cake_pre, MATRIX& Rcake_div_Rmem_arr, MATRIX& Pi_div_DLP_arr, MATRIX& yt_arr, long n_iter, DISP_TYPE& DISP, PARAMETER_SET& PARS, AUX_MATRIX& AUX_h1, AUX_MATRIX& AUX_h2);

int shear_rate_initialization(MATRIX& shear_rate_nodim_ZY, MATRIX& yt_arr, PARAMETER_SET& PARS);

int shear_stress_initialization (MATRIX& shear_stress_nodim_ZY, MATRIX& shear_rate_nodim_ZY, MATRIX& yt_arr, PARAMETER_SET& PARS, DISP_TYPE& DISP);

int update_shear_rate_nodim_ZY(MATRIX& shear_rate_nodim_ZY, MATRIX& shear_rate_nodim_pre_ZY, MATRIX& shear_stress_nodim_pre_ZY, MATRIX& Ieta_yt_arr_ZY, MATRIX& yt_arr, MATRIX& phi_yt_div_phib_arr, MATRIX& uz_div_u0_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS);

int update_shear_stress_nodim_ZY(MATRIX& shear_stress_nodim_ZY, MATRIX& shear_stress_nodim_pre_ZY, MATRIX& shear_rate_nodim_ZY, MATRIX& phi_yt_div_phib_arr, DISP_TYPE& DISP, PARAMETER_SET& PARS);

inline long check_min(long a, long b)
{
  if (a < b)
    return a;
  return b;
	   
}

#endif
