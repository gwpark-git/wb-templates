
#ifndef TRANSPORT_PROPERTIES_SID_H
#define TRANSPORT_PROPERTIES_SID_H

#include "transport_properties_SPHS.h" // some function will be used internally
#include "parameters.h"

namespace SID
{

  // double fcn_unity(double phi, PARAMETER_SET& PARS, double shear_rate_nodim, double relative_viscosity = 0.);
  double fcn_unity(double phi, PARAMETER_SET& PARS, double shear_rate_nodim);
  
  double eta_div_eta0_KD(double phi, double phi_max, double q);
  double eta_div_eta0_test(double phi, PARAMETER_SET& PARS, double shear_stress_nodim, double temp = 0.);  
  double eta_div_eta0_test_EQ(double phi, PARAMETER_SET& PARS, double shear_rate_nodim, double relative_viscosity = 0.);
  
  double delta_g_iso_contact(double phi, double shear_rate_nodim);
  double Dc_short_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_nodim, double relative_viscosity = 0.);


  
}
#endif
