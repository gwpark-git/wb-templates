
#ifndef DISPERSION_PROPERTIES_H
#define DISPERSION_PROPERTIES_H

#include "parameters.h"
#include <memory>
#include "case_driver.h"

#include "transport_properties_SPHS.h"
#include "transport_properties_INTERPOLATE.h"
#include "transport_properties_SID.h"
#include "osmotic_pressure_CS.h"

extern "C"
{
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
}

struct SplineDeleter { void operator()(gsl_spline* p) const { gsl_spline_free(p); } };
struct AccelDeleter { void operator()(gsl_interp_accel* p) const { gsl_interp_accel_free(p); } };

#define True 1
#define False 0

class DISP_TYPE 
{

public:

  // Note: the varialbe "shear_rate_or_stress_nodim" can be either for the shear rate or shear stress
  // This depends on the selection of dispersion properties
  
  inline double fcn_Dc_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    if (transport_) return evaluate_live(phi).collective_diffusivity_over_D0;
    return (*ptr_Dc_div_D0)(phi, PARS, shear_rate_or_stress_nodim, additional_arg, spline_Dc_div_D0.get(), acc_spline_Dc_div_D0.get());
  }

  inline double fcn_eta_div_eta0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    if (transport_) return evaluate_live(phi).viscosity_over_solvent;
    return (*ptr_eta_div_eta0)(phi, PARS, shear_rate_or_stress_nodim, additional_arg, spline_eta_div_eta0.get(), acc_spline_eta_div_eta0.get());
  }

  inline double fcn_Pi_div_Piast(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return (*ptr_Pi_div_Piast)(phi, PARS, shear_rate_or_stress_nodim,  additional_arg, spline_Pi.get(), acc_spline_Pi.get());
  }
  
  
  DISP_TYPE() = delete;
  
  DISP_TYPE(PARAMETER_SET& PARS, mbla::EquilibriumTransport* transport = nullptr);
  DISP_TYPE(const DISP_TYPE&) = delete;
  DISP_TYPE& operator=(const DISP_TYPE&) = delete;
  ~DISP_TYPE();


private:
  mbla::TransportValues evaluate_live(double phi);
  mbla::EquilibriumTransport* transport_;

  double (*ptr_Dc_div_D0)
  (double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline_Dc_div_D0);
  
  double (*ptr_eta_div_eta0)
  (double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline_eta_div_eta0);
  
  double (*ptr_Pi_div_Piast)
  (double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline_Pi);



    // When there is given file for each of suspension properties
  // the following will be used
  bool IDENT_file_D = FALSE;
  bool IDENT_file_eta = FALSE;
  bool IDENT_file_Pi = FALSE;
  bool IDENT_any_file = FALSE;

  
  // this will give a interpolation function
  std::unique_ptr<gsl_spline, SplineDeleter> spline_Dc_div_D0; // when the data_given
  std::unique_ptr<gsl_spline, SplineDeleter> spline_eta_div_eta0;
  std::unique_ptr<gsl_spline, SplineDeleter> spline_Pi;

  // gsl_interp_accel *acc_spline;
  /*
    Essentially, we can use single object gsl_interp_accel for all of evaluation step.
    Which, however, is necessary to be reset this object whenever we did the evaluation.
    This means, we do not gain any benefit to use acceleration object.
    Therefore, here we allocate each of accerleation for each of spline object, then reset only when it is necessary.

   */
  std::unique_ptr<gsl_interp_accel, AccelDeleter> acc_spline_Dc_div_D0;
  std::unique_ptr<gsl_interp_accel, AccelDeleter> acc_spline_eta_div_eta0;
  std::unique_ptr<gsl_interp_accel, AccelDeleter> acc_spline_Pi;
  
  
};

// inlined wrapper functions
namespace WRP
{

  /*
    Below, the free parameter is specified as 1.0 without any meaning for equilibrium transport properties
    However, it will be importance for the case when shear-induced migration where the dimensionless_shear_rate should be provided.

    The default value for shear_rate_or_stress_nodim is NOT set here. 
    This is because the sequence of default value must be the end of parameter list.
    Rather than specify default value here, we put rigorous pointer_function whether shear-rate is accounted or not.
   */

  /*
    For SID namespace, (SID stands for shear-induced diffusivity) the original name comes from SID, but it should be replaced by SIM.
    For diffusivity, the current model uses dressed Peclet number, which in other word, a dimensionless shear rate.
    For viscosity, the Cross-like equation uses dimensionless shear stress instead of shear rate.
    Therefore, the additional arguments should pass information about viscosity (of previous iteration)
   */
  
  inline double SID_eta_div_eta0_test(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    // the given shear_rate_or_stress_nodim = shear_rate/characteristic_shear_rate
    // required shear_rate_or_stress_nodim inside the function is shear_rate * a_H^2 /D0.
    // Therefore, multiplication with PARS.Pe_a is required here.
    // return SID::eta_div_eta0(phi, PARS, shear_rate_or_stress_nodim*PARS.Pe_a);

    // additional_arg == relative_viscosity
    return SID::eta_div_eta0_test(phi, PARS, shear_rate_or_stress_nodim, additional_arg);
  }

  inline double SID_eta_div_eta0_test_EQ(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return SID::eta_div_eta0_test_EQ(phi, PARS, shear_rate_or_stress_nodim, additional_arg);
  }
  
  inline double SID_Dc_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return SID::Dc_short_div_D0(phi, PARS, shear_rate_or_stress_nodim);

    
  }
  
  inline double SPHS_fcn_unity(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return SPHS::fcn_unity(phi, PARS);
  }

  inline double SPHS_eta_div_eta0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return SPHS::eta_div_eta0(phi, PARS);
  }

  inline double HS_eta_div_eta0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return HS::eta_div_eta0(phi, PARS);
  }

  inline double SPHS_Dc_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return SPHS::Dc_short_div_D0(phi, PARS);
  }

  inline double CS_Pi_div_Piast(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return CS::get_Pi_div_Piast(phi, PARS);
  }

  inline double CS_fcn_zero(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim=0., double additional_arg=0., gsl_spline *spline_cubic=NULL, gsl_interp_accel *acc_spline=NULL)
  {
    return CS::fcn_zero(phi, PARS);
  }

  inline double INTERP_eta_div_eta0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline)
  {
    return INTERP::eta_div_eta0(phi, PARS, spline_cubic, acc_spline);
  }

  inline double INTERP_Dc_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline)
  {
    return INTERP::Dc_div_D0(phi, PARS, spline_cubic, acc_spline);
  }

  inline double INTERP_Pi_div_Piast(double phi, PARAMETER_SET& PARS, double shear_rate_or_stress_nodim, double additional_arg, gsl_spline *spline_cubic, gsl_interp_accel *acc_spline)
  {
    return INTERP::Pi_div_Piast(phi, PARS, spline_cubic, acc_spline);
  }

  

}

#endif
