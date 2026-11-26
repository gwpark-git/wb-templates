
#include "transport_properties_SPHS.h"

// using namespace SPHS;

double SPHS::fcn_unity(double phi, PARAMETER_SET& PARS)
{
  return 1.;
}

// Combine the high-frequency and structural viscosity terms, capped at the freeze fraction.
double SPHS::eta_div_eta0(double phi, PARAMETER_SET& PARS)
{

  if (phi < PARS.phi_freeze)
    {
      return SPHS::eta_inf_div_eta0(phi, PARS.gamma)*(1. + (1./SPHS::Gamma_S(phi, PARS.gamma))*SPHS::Del_eta_noHI_div_eta0(phi));
    }
  return SPHS::eta_inf_div_eta0(PARS.phi_freeze, PARS.gamma)*(1. + (1./SPHS::Gamma_S(PARS.phi_freeze, PARS.gamma))*SPHS::Del_eta_noHI_div_eta0(PARS.phi_freeze));
}

// Use sedimentation divided by the equilibrium structure factor for relative collective diffusion.
double SPHS::Dc_short_div_D0(double phi, PARAMETER_SET& PARS)
{
  if (phi < PARS.phi_freeze)
    {
      return SPHS::K(phi, PARS.gamma)/SPHS::S0_CS(phi);
    }
  return SPHS::K(PARS.phi_freeze, PARS.gamma)/SPHS::S0_CS(PARS.phi_freeze);
}

double SPHS::Ds_div_S0(double phi, double gamma)
{
  return 1 + SPHS::lambda_t(gamma)*phi*(1 + 0.12*phi - 0.70*pow(phi,2.0));
}

double SPHS::lambda_t(double gamma)
{
  return -1.8315 + 7.820 *(1-gamma) - 14.231 *pow(1-gamma, 2.0) + 14.908* pow(1-gamma, 3.0)- 9.383* pow(1-gamma, 4.0) + 2.717* pow(1-gamma, 5.0);
}

// Keep the inherited multiplicative high-frequency formula for the spherical-shell model.
double SPHS::eta_inf_div_eta0(double phi, double gamma)
{
  return 1 + 2.5*pow(gamma,3.0)*phi*(1 + SPHS::Saito_fcn(phi, gamma))*(1 - pow(gamma, 3.0) * phi * (1 + SPHS::Saito_fcn(phi, gamma)));
}

double SPHS::Saito_fcn(double phi, double gamma)
{
  // this is a compatible version for solvent permeable case
  // however, for the impermeable hard sphere, gamma = 1.0, the value may produe different
  return phi*(SPHS::lambda_V(gamma)/(2.5*pow(gamma, 3.0)) - pow(gamma, 3.0));
  
}



double SPHS::lambda_V(double gamma)
{
  return 5.0021 - 39.279* (1-gamma) + 143.179 *pow(1-gamma, 2.0) - 288.202* pow(1-gamma, 3.0) + 254.581 *pow(1- gamma, 4.0);  
}

double SPHS::Del_eta_noHI_div_eta0(double phi)
{
  return (12/5.)*pow(phi, 2.0)*(1 - 7.085*phi + 20.182*pow(phi, 2.0))/(1 - phi/0.64);

}

// Construct the relaxation factor from relative self-diffusion and high-frequency viscosity.
double SPHS::Gamma_S(double phi, double gamma)
{
  return SPHS::Ds_div_S0(phi, gamma)*SPHS::eta_inf_div_eta0(phi, gamma);
}

// Evaluate the sedimentation correlation at volume fraction phi and radius ratio gamma.
double SPHS::K(double phi, double gamma)
{
  double gp = gamma*phi;
  return 1 + SPHS::lambda_K(gamma)*phi*(1 - 3.348*(gp) + 7.426*pow(gp, 2.0) - 10.034*pow(gp, 3.0) + 5.882*pow(gp, 4.0));
  
}

double SPHS::lambda_K(double gamma)
{
  return -6.5464 + 8.592*(1-gamma) - 3.901*pow(1-gamma, 2.0) + 2.011*pow(1-gamma, 3.0) - 0.142*pow(1-gamma, 4.0);
}

double SPHS::S0_CS(double phi)
{
  return pow(1.-phi, 4.0)/(pow(1.+2.*phi, 2.0) + pow(phi, 3.0)*(phi - 4));
}

// below is SPHS::HS:: for the case of correct Saito function

double HS::Saito_fcn(double phi, double gamma)
{
  // this is an accurate version for hard spheres
  return phi*(1. + 0.95*phi - 2.15 * phi*phi);
}

double HS::eta_inf_div_eta0(double phi, double gamma)
{
  return 1 + 2.5*pow(gamma,3.0)*phi*(1 + HS::Saito_fcn(phi, gamma))*(1 - pow(gamma, 3.0) * phi * (1 + HS::Saito_fcn(phi, gamma)));
}

// Apply the hard-sphere Saito variant while retaining the inherited structural factor.
double HS::eta_div_eta0(double phi, PARAMETER_SET& PARS)
{

  if (phi < PARS.phi_freeze)
    {
      return HS::eta_inf_div_eta0(phi, PARS.gamma)*(1. + (1./SPHS::Gamma_S(phi, PARS.gamma))*SPHS::Del_eta_noHI_div_eta0(phi));
    }
  return HS::eta_inf_div_eta0(PARS.phi_freeze, PARS.gamma)*(1. + (1./SPHS::Gamma_S(PARS.phi_freeze, PARS.gamma))*SPHS::Del_eta_noHI_div_eta0(PARS.phi_freeze));
}
