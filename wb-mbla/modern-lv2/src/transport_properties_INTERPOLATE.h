
#ifndef TRANSPORT_PROPERTIES_INTERPOLATE_H
#define TRANSPORT_PROPERTIES_INTERPOLATE_H

#include "parameters.h"
#include "matrix.h" 
extern "C"
{
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
}

namespace INTERP
{
  int interpolate_spline(gsl_spline* spline_cubic, MATRIX& x_vec, MATRIX& y_vec);
  int interpolate_spline(gsl_spline* spline_cubic, MATRIX& xy_mat);
  
  double Dc_div_D0(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline);
  double eta_div_eta0(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline);
  double Pi_div_Piast(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline);

  
}


#endif
