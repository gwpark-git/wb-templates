
#include "transport_properties_SID.h"
#include "math.h" // for M_PI constant

// using namespace SID;

double SID::fcn_unity(double phi, PARAMETER_SET& PARS, double shear_rate_nodim)
{
  return 1.;
}

// Return the Krieger–Dougherty relative viscosity for the supplied packing fraction and exponent.
double SID::eta_div_eta0_KD(double phi, double phi_max, double q)
{
  return pow(1. - phi/phi_max, -q);
}

// Use the inherited TEST stress law with b = 1 and the raw phi in its two viscosity limits.
double SID::eta_div_eta0_test(double phi, PARAMETER_SET& PARS, double shear_stress_nodim, double temp)
{
    /*
      \tilde{\dot{gamma}} = \frac{a^2}{D_0}\dot{\gamma} = \frac{a^3}{kT}6\pi \eta_s \dot{\gamma}
      \tilde{\tau} = \frac{a^3}{kT}\tau \approx \frac{a^3}{kT}\eta \dot{\gamma}
      = \tilde{\eta}\tilde{\dot{gamma}}/(6\pi).
     */
  // double shear_stress_nodim = shear_rate_nodim * relative_viscosity / (6.* M_PI);
  double tmp_phi = phi;  
  if (phi > PARS.phi_freeze)
    {
      tmp_phi = PARS.phi_freeze;
    }
    
  double etar_0 = eta_div_eta0_KD(phi, 0.63, 2.0);
  double etar_inf = eta_div_eta0_KD(phi, 0.71, 2.0);
  return etar_inf + (etar_0 - etar_inf)/(1. + shear_stress_nodim);
}

// The TEST_EQ law uses the zero-stress viscosity limit based on the raw input volume fraction.
double SID::eta_div_eta0_test_EQ(double phi, PARAMETER_SET& PARS, double shear_rate_nodim, double relative_viscosity)
{
  // if (phi < PARS.phi_freeze)
  //   {
  //     return eta_inf_div_eta0(phi, PARS.gamma)*(1. + (1./Gamma_S(phi, PARS.gamma))*Del_eta_noHI_div_eta0(phi));
  //   }
  // return eta_inf_div_eta0(PARS.phi_freeze, PARS.gamma)*(1. + (1./Gamma_S(PARS.phi_freeze, PARS.gamma))*Del_eta_noHI_div_eta0(PARS.phi_freeze));

  /*
    This is still a test function. Note that from relation Eq. (22) and the following viscosity relation in manuscript by Jin et al. (2014), the shear-curvature viscosity plays a role. However, at this moment (21 JAN 2022), I have no clue the range of k. Instead of including this part, I simply ignore the term. Note that this is also the same for the gradient diffusion coefficient, although this is done after the dominant balance analysis.
   */
  double tmp_phi = phi;
  double phi_max_RCP = 0.63;
  double q = 2.0; // exponent in Krieger-Dougherty relation
  if (phi > PARS.phi_freeze)
    {
      tmp_phi = PARS.phi_freeze;
    }
  return pow(1. - phi/phi_max_RCP, -q);
}



// Compute the shear-dependent isotropic contact correction used in the diffusivity model.
double SID::delta_g_iso_contact(double phi, double shear_rate_nodim)
{
  double phi_max = 0.64; // RCP
  double s = 2.525;
  double m = 0.43 + 5.26*phi - 8.80*phi*phi;
  double A = 0.0140;
  double re = A*pow(shear_rate_nodim, m)*pow(1. - phi/phi_max, -s);
  return re;
}

// Add the nonnegative shear-induced contribution to equilibrium diffusivity above the shear threshold.
double SID::Dc_short_div_D0(double phi, PARAMETER_SET& PARS, double shear_rate_nodim, double relative_viscosity)
{
  /*
    Temporal approach at this moment (21 JAN 2022)
    Note that the shear-curavture diffusivity is ignored due to the fact that the effect within the thin CP layer is negligible.
    Therefore, the effective diffusivity only accounts the isotropic contributions from both of equilibrium and shear-induced parts.

    More importantly, we use EQUILIBRIUM sedimentation coefficient K(phi) for the both of equilibrium and isotropic parts of shear-induced diffusivity.

   */
  if (phi > PARS.phi_freeze)
    {
      phi = PARS.phi_freeze;
    }

  double K_eq_phi = SPHS::K(phi, 1.0); // gamma=1 since we do not consider the solvent-permeable hard spheres
  double D_div_D0_eq_phi = K_eq_phi / (SPHS::S0_CS(phi));
  
  double deltaD_div_D0_SID = 0.;

  if (shear_rate_nodim > 1e-2)
    {
      double phi_max = 0.64; // RCP
      double s = 2.525;
      double m = 0.43 + 5.26*phi - 8.80*phi*phi;
      double dm_dphi = 5.26 - 17.6*phi;
      if (dm_dphi < 0.)
	{
	  dm_dphi = 0.;
	}
  
      double A = 0.0140;

      double intermediate_term = dm_dphi*log(shear_rate_nodim) + s/(phi_max - phi);
      deltaD_div_D0_SID = K_eq_phi*(8*phi + 4.*phi*phi*intermediate_term)*SID::delta_g_iso_contact(phi, shear_rate_nodim);
      if(deltaD_div_D0_SID < 0)
	deltaD_div_D0_SID = 0.;
    }
  return D_div_D0_eq_phi + deltaD_div_D0_SID;
  
}
