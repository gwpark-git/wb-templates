
#ifndef OSMOTIC_PRESSURE_CS
#define OSMOTIC_PRESSURE_CS

#include "parameters.h"

namespace CS
{
  double fcn_zero(double phi, PARAMETER_SET& PARS);
  double Pi_CS(double phi, PARAMETER_SET& PARS);
  double CS_Z(double phi);
  double get_Pi(double phi, PARAMETER_SET& PARS);
  double get_Pi_div_Piast(double phi, PARAMETER_SET& PARS);
  double CS_S0(double phi);
}


#endif
