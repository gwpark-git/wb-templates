
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <math.h>
// #include <mkl.h>
#include <string>

#include "file_IO.h"

extern "C"{
#include <gsl/gsl_const_mksa.h> // physical constants
#include <gsl/gsl_math.h> // mathematical constants
}

#define TRUE 1
#define FALSE 0

struct PARAMETER_SET
{
  // bool linear_P_profile;
  bool DETAILED_REPORT = FALSE;
  long N_iter;
  double weight;
  double weight_f;

  double R_channel, L_channel;
  double Lp, eta0, T;
  
  double k, preU;
  double DLP, DTP_HP, DTP;
  double Pin, Pout, Pper;

  double a, a_H;
  double Va, kT;
  double D0;

  double u0, vw0;
  double alpha_ast, beta_ast;
  double dg_ast, tau_ast;

  long Nz, Nr_bulk, Ny;
  long count_negative_phiw = 0;
  double dz, dr;
  double gamma;
  double epsilon_d;
  double phi_freeze, phi_bulk;
  double phiw_init; // this is the initial condition except at the inlet

  // double phi_cake, epsilon_f, delta_f_ast, porosity, zeta_f_KC; // no more necessary for the current modeling
  bool CAKE_THICKNESS=FALSE;
  double CAKE_convert_unit=0.; // (h'/R)*(kappa_cake / kappa_mem) * xi_tilde
  double phi_cake = 0.64;
  double h_mem=0.;
  
  double R_mem;

  double Pi_ast;

  double TOL_phiw;
  double TOL_delta_f;
  
  string str_Dc_div_D0;
  string str_eta_div_eta0;
  string str_Pi;

  string inlet_BC; // this will record either the velocity or pressure boundary conditions
  bool legacy_inlet_velocity = false; // First-order pressure-gradient replay for frozen outputs.

  double shear_rate_ast = 0.;
  double Pe_a; // single-particle shear Peclet number of inlet-center
  
  int init(COND& given_cond);
  int print();
};

double get_Pin(double DLP, double Pout);
double get_Pper(double DLP, double DTP, double k, double Pout);


#endif 

