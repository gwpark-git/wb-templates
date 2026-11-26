
#include "transport_properties_INTERPOLATE.h"

namespace {
// Convert GSL domain or nonfinite results into a numerical input error.
double eval(gsl_spline* spline, double phi, gsl_interp_accel* accelerator) {
  double value = 0.;
  if (gsl_spline_eval_e(spline, phi, accelerator, &value) != GSL_SUCCESS || !std::isfinite(value))
    throw std::runtime_error("property spline evaluation outside its domain or nonfinite");
  return value;
}
}


int INTERP::
interpolate_spline(gsl_spline* spline_cubic, MATRIX& x_vec, MATRIX& y_vec)
{
  const size_t N = x_vec.size;

  if (gsl_spline_init(spline_cubic, x_vec.data, y_vec.data, N) != GSL_SUCCESS)
    throw std::runtime_error("invalid cubic property table");

  return 0;
}

int INTERP::
interpolate_spline(gsl_spline* spline_cubic, MATRIX& xy_mat)
{
  const size_t N = xy_mat.rows;
  MATRIX x_vec(N, 1, 0.);
  MATRIX y_vec(N, 1, 0.);

  for(long i=0; i<N; i++)
    {
      x_vec(i) = xy_mat(i, 0);
      y_vec(i) = xy_mat(i, 1);
    }

  if (gsl_spline_init(spline_cubic, x_vec.data, y_vec.data, N) != GSL_SUCCESS)
    throw std::runtime_error("invalid cubic property table");
  return 0;
}


// with gsl_spline
// Interpolate the supplied relative diffusivity on the capped concentration domain.
double INTERP::
Dc_div_D0(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline)
{
  if (phi < PARS.phi_freeze)
    {
      return eval(spline_cubic, phi, acc_spline);
    }
  return eval(spline_cubic, PARS.phi_freeze, acc_spline);
}

// Interpolate relative viscosity and cap the query at the configured freeze fraction.
double INTERP::
eta_div_eta0(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline)
{
  
  if (phi < PARS.phi_freeze)
    {
      return eval(spline_cubic, phi, acc_spline);
    }
  return eval(spline_cubic, PARS.phi_freeze, acc_spline);
}

// Interpolate wall osmotic pressure in Pa and normalize it by the case pressure scale.
double INTERP::
Pi_div_Piast(double phi, PARAMETER_SET& PARS, gsl_spline* spline_cubic, gsl_interp_accel *acc_spline)
{
  if (phi < PARS.phi_freeze)
    {
      return eval(spline_cubic, phi, acc_spline)/PARS.Pi_ast;
    }
  return eval(spline_cubic, PARS.phi_freeze, acc_spline)/PARS.Pi_ast;
}

