
#ifndef TRANSPORT_PROPERTIES_SPHS_H
#define TRANSPORT_PROPERTIES_SPHS_H

#include "parameters.h"

namespace SPHS
{
  double fcn_unity(double phi, PARAMETER_SET& PARS);
  double eta_div_eta0(double phi, PARAMETER_SET& PARS);
  double Dc_short_div_D0(double phi, PARAMETER_SET& PARS);
  double Ds_div_S0(double phi, double gamma);
  double lambda_t(double gamma);
  double eta_inf_div_eta0(double phi, double gamma);
  double Saito_fcn(double phi, double gamma);
  double lambda_V(double gamma);
  double Del_eta_noHI_div_eta0(double phi);
  double Gamma_S(double phi, double gamma);
  double K(double phi, double gamma);
  double lambda_K(double gamma);
  double S0_CS(double phi);

}

  namespace HS
  {
    double Saito_fcn(double phi, double gamma=1.0);
    double eta_inf_div_eta0(double phi, double gamma=1.0);
    double eta_div_eta0(double phi, PARAMETER_SET& PARS);      
  }


#endif
