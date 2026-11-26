
#include "dispersion_properties.h"

namespace {
// Validate a two-column interpolation table against the active concentration range.
void check_table(MATRIX& table, const PARAMETER_SET& pars, bool positive) {
  if (table.rows < 3 || table.cols != 2) throw std::runtime_error("cubic property table requires >=3 rows and exactly two columns");
  if (table(0,0) > 0 || table(table.rows-1,0) < pars.phi_freeze) throw std::runtime_error("property table must cover [0, phi_freeze]");
  for (long i=0; i<table.rows; ++i) {
    if (i && table(i,0) <= table(i-1,0)) throw std::runtime_error("property table concentration must increase strictly");
    if (positive && table(i,1) <= 0) throw std::runtime_error("diffusion/viscosity table values must be positive");
  }
}
}

// Select each named closure or own its cubic spline and lookup accelerator for the case.
DISP_TYPE::
DISP_TYPE(PARAMETER_SET& PARS, mbla::EquilibriumTransport* transport)
  : transport_(transport)
{
  const bool live_diffusion = PARS.str_Dc_div_D0 == "LIVE";
  const bool live_viscosity = PARS.str_eta_div_eta0 == "LIVE";
  if (live_diffusion != bool(transport_) || live_viscosity != bool(transport_)) {
    throw std::runtime_error("LIVE transport requires both selections and an injected equilibrium provider");
  }
  if (transport_) transport_->validate_case(PARS.phi_freeze, PARS.Pe_a);

  IDENT_file_D = False;
  IDENT_file_eta = False;
  IDENT_file_Pi = False;
    
  if (live_diffusion)
    ptr_Dc_div_D0 = nullptr;  // the paired provider owns evaluation in this mode
  else if (PARS.str_Dc_div_D0=="SPHS" || PARS.str_Dc_div_D0=="HS")
    {
      ptr_Dc_div_D0 = &WRP::SPHS_Dc_div_D0;
    }
  else if(PARS.str_Dc_div_D0=="SID")
    {
      ptr_Dc_div_D0 = &WRP::SID_Dc_div_D0;
    }
  else if(PARS.str_Dc_div_D0=="CT")
    {
      ptr_Dc_div_D0 = &WRP::SPHS_fcn_unity;
    }
  else
    {
      IDENT_file_D = True;
      IDENT_any_file = True;

      MATRIX data_phi_vs_Dc_div_D0(PARS.str_Dc_div_D0.c_str());
      check_table(data_phi_vs_Dc_div_D0, PARS, true);
      spline_Dc_div_D0.reset(gsl_spline_alloc(gsl_interp_cspline, data_phi_vs_Dc_div_D0.rows));
      
      if (!spline_Dc_div_D0) throw std::bad_alloc();
      INTERP::interpolate_spline(spline_Dc_div_D0.get(), data_phi_vs_Dc_div_D0);
      ptr_Dc_div_D0 = &WRP::INTERP_Dc_div_D0;
      acc_spline_Dc_div_D0.reset(gsl_interp_accel_alloc());
      if (!acc_spline_Dc_div_D0) throw std::bad_alloc();
    }

  if (live_viscosity)
    ptr_eta_div_eta0 = nullptr;
  else if (PARS.str_eta_div_eta0=="SPHS")
    {
      ptr_eta_div_eta0 = &WRP::SPHS_eta_div_eta0;
    }
  else if (PARS.str_eta_div_eta0=="HS")
    {
      ptr_eta_div_eta0 = &WRP::HS_eta_div_eta0;
    }
  else if(PARS.str_eta_div_eta0=="TEST")
    {
      ptr_eta_div_eta0 = &WRP::SID_eta_div_eta0_test;
    }
  else if(PARS.str_eta_div_eta0=="TEST_EQ")
    {
      ptr_eta_div_eta0 = &WRP::SID_eta_div_eta0_test_EQ;
    }
  else if(PARS.str_eta_div_eta0=="CT")
    {
      ptr_eta_div_eta0 = &WRP::SPHS_fcn_unity;
    }
  else
    {
      IDENT_file_eta = True;
      IDENT_any_file = True;
      MATRIX data_phi_vs_eta_div_eta0(PARS.str_eta_div_eta0.c_str());
      check_table(data_phi_vs_eta_div_eta0, PARS, true);
      spline_eta_div_eta0.reset(gsl_spline_alloc(gsl_interp_cspline, data_phi_vs_eta_div_eta0.rows));
      if (!spline_eta_div_eta0) throw std::bad_alloc();
      INTERP::interpolate_spline(spline_eta_div_eta0.get(), data_phi_vs_eta_div_eta0);
      ptr_eta_div_eta0 = &WRP::INTERP_eta_div_eta0;
      acc_spline_eta_div_eta0.reset(gsl_interp_accel_alloc());
      if (!acc_spline_eta_div_eta0) throw std::bad_alloc();
    }
  if (PARS.str_Pi=="CS")
      ptr_Pi_div_Piast = &WRP::CS_Pi_div_Piast;
  else if (PARS.str_Pi=="NO")
    ptr_Pi_div_Piast = &WRP::CS_fcn_zero;
  else
    {
      IDENT_file_Pi = True;
      IDENT_any_file = True;
      MATRIX data_phi_vs_Pi(PARS.str_Pi.c_str());
      check_table(data_phi_vs_Pi, PARS, false);
      spline_Pi.reset(gsl_spline_alloc(gsl_interp_cspline, data_phi_vs_Pi.rows));
      if (!spline_Pi) throw std::bad_alloc();
      INTERP::interpolate_spline(spline_Pi.get(), data_phi_vs_Pi);
      ptr_Pi_div_Piast = &WRP::INTERP_Pi_div_Piast;
      acc_spline_Pi.reset(gsl_interp_accel_alloc());
      if (!acc_spline_Pi) throw std::bad_alloc();
    }
}

// Each spline/accelerator has a unique owner, including constructor failure.
DISP_TYPE::~DISP_TYPE() = default;

// Validate the law's reduced values before the solver uses their reciprocals.
mbla::TransportValues DISP_TYPE::evaluate_live(double phi) {
  const auto values = transport_->evaluate(phi);
  const double d = values.collective_diffusivity_over_D0;
  const double eta = values.viscosity_over_solvent;
  if (!std::isfinite(d) || !std::isfinite(eta) || d <= 0. || eta <= 0.) {
    throw std::runtime_error("LIVE provider returned invalid relative transport values");
  }
  return values;
}
