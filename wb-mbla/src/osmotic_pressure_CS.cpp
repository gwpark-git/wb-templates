
#include "osmotic_pressure_CS.h"

double CS::fcn_zero(double phi, PARAMETER_SET& PARS)
{
  return 0.;
}


double CS::Pi_CS(double phi, PARAMETER_SET& PARS)
{
  return CS::get_Pi(phi, PARS);
  
}

double CS::CS_Z(double phi)
{
  return (1 + phi + pow(phi, 2.0) - pow(phi, 3.0))/pow(1.-phi, 3.0);
}

double CS::get_Pi(double phi, PARAMETER_SET& PARS)
{
  double rho = phi/PARS.Va;
  return rho*PARS.kT*CS::CS_Z(phi);
}

double CS::get_Pi_div_Piast(double phi, PARAMETER_SET& PARS)
{
  return phi * CS::CS_Z(phi);
}

double CS::CS_S0(double phi)
{
  return pow(1. - phi, 4.0) / (pow(1. + 2.*phi, 2.0) + pow(phi, 3.0) *(phi - 4.));
}


