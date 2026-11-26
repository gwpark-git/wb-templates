#include "ported_lib.h"
#include "fixed_point_controls.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok, const char* description) {
  if (!ok) throw std::runtime_error(description);
}
void near(double x, double y, double tolerance, const char* description) {
  require(std::isfinite(x) && std::abs(x-y) <= tolerance, description);
}
// This closed-form pressure reference is independent of the finite-difference
// assembly and LU implementation. It tests both boundaries and the sign of k².
double pressure_error(long nz, bool velocity=false, double cake_ratio=0.) {
  PARAMETER_SET p{};
  p.Nz=nz; p.L_channel=0.5; p.dz=p.L_channel/(nz-1);
  p.DLP=130; p.Pin=101455; p.Pout=101325; p.Pper=96400; p.k=0.4;
  p.inlet_BC=velocity ? "velocity" : "pressure";
  MATRIX pressure(nz,1), ln(nz,1,1./(1.+cake_ratio)), pi(nz,1,0.), a(nz,nz), y(nz,1), d(nz,1),
    r(nz,1), q(nz,1), f(nz,1), cake(nz,1);
  get_P_div_DLP_arr_with_fouling(pressure,p,ln,pi,a,y,d,r,q,f,cake);
  double error=0.;
  for (long i=0;i<nz;++i) {
    const double z=double(i)/(nz-1);
    const double rate=p.k/std::sqrt(1.+cake_ratio);
    double exact;
    if (velocity) {
      // q'(0)=-DLP and q(1)=Pout-Pper, q=P-Pper. This exact solution
      // exercises the fictitious inlet row independently of its assembly.
      const double c=(p.Pout-p.Pper+p.DLP*std::sinh(rate)/rate)/std::cosh(rate);
      exact=c*std::cosh(rate*z)-p.DLP*std::sinh(rate*z)/rate+p.Pper;
    } else {
      exact=((p.Pin-p.Pper)*std::sinh(rate*(1-z))+
             (p.Pout-p.Pper)*std::sinh(rate*z))/std::sinh(rate)+p.Pper;
    }
    error=std::max(error,std::abs(pressure(i)*p.DLP-exact));
  }
  if (!velocity) near(pressure(0)*p.DLP,p.Pin,2e-8,"inlet pressure");
  near(pressure(nz-1)*p.DLP,p.Pout,2e-8,"outlet pressure, zero osmotic source");
  return error;
}

// The exact gradients of a cubic test the production one-sided stencil
// without depending on pressure assembly or any frozen simulation output.
double inlet_velocity_error(long nz, bool legacy) {
  PARAMETER_SET p{};
  p.Nz = nz;
  p.L_channel = 0.5;
  p.dz = p.L_channel/(nz-1);
  p.inlet_BC = "pressure";
  p.legacy_inlet_velocity = legacy;
  p.CAKE_convert_unit = 0.;
  MATRIX pressure(nz,1), z(nz,1), uz(nz,1), cake(nz,1,0.);
  for (long i=0; i<nz; ++i) {
    const double x = double(i)/(nz-1);
    z(i) = x;
    pressure(i) = 2. - 3.*x + 4.*x*x + 5.*x*x*x;
  }
  get_uz_div_u0_arr_with_fouling(uz,z,pressure,cake,p);
  const double error = uz(0) - 3.;
  p.inlet_BC = "velocity";
  get_uz_div_u0_arr_with_fouling(uz,z,pressure,cake,p);
  near(uz(0),1.,0.,"prescribed inlet velocity must be used exactly");
  return error;
}
}
int main() {
  gsl_set_error_handler_off();
  try {
    MATRIX m(2,2,0.); m(0,0)=7; m(1,1)=-3;
    m=m; near(m(0,0),7,0,"self assignment must retain data");
    MATRIX copy(m); copy(0,0)=2;
    near(m(0,0),7,0,"matrix copy must not alias storage");
    MATRIX assigned; assigned=copy; copy.initial(3,1,8.);
    near(assigned(0,0),2,0,"assignment must own independent storage");
    bool rejected=false;try { m(2,0)=1; } catch(const std::out_of_range&) { rejected=true; }
    require(rejected,"bounds must reject invalid row");
    MATRIX a(2,2),b(2,1),x(2,1);
    a(0,0)=3; a(0,1)=2; a(1,0)=1; a(1,1)=2; b(0)=7; b(1)=5;
    gsl_sol_Ax_b(x,a,b);near(x(0),1,1e-13,"LU x");near(x(1),2,1e-13,"LU y");
    PARAMETER_SET p{};p.kT=2.;p.Va=4.;p.gamma=1.;p.phi_freeze=0.494;
    // At phi=1/5, Carnahan–Starling Z=2.40625; kT/Va is 1/2.
    near(CS::CS_Z(0.2),2.40625,1e-14,"CS compressibility factor");
    near(CS::get_Pi(0.2,p),0.240625,1e-14,"osmotic pressure unit conversion");
    near(SPHS::Dc_short_div_D0(0,p),1,1e-14,"dilute diffusivity");
    near(SPHS::eta_div_eta0(0,p),1,1e-14,"dilute viscosity");
    near(get_nodim_Lpf(0),1,0,"clean membrane permeability");
    near(get_nodim_Lpf(3),0.25,0,"series resistance law");
    require(get_nodim_Lpf(2)>get_nodim_Lpf(3),"increased cake must reduce permeability");
    const double coarse=pressure_error(25),fine=pressure_error(49);
    require(fine<0.001 && coarse/fine>3.8 && coarse/fine<4.2,"second order axial pressure convergence");
    const double cake_coarse=pressure_error(25,false,3.),cake_fine=pressure_error(49,false,3.);
    const double velocity_coarse=pressure_error(25,true),velocity_fine=pressure_error(49,true);
    require(cake_fine<0.001 && cake_coarse/cake_fine>3.8 && cake_coarse/cake_fine<4.2,
            "second order pressure convergence with prescribed uniform cake permeability");
    require(velocity_fine<0.01 && velocity_coarse/velocity_fine>3.8 && velocity_coarse/velocity_fine<4.2,
            "second order velocity-inlet pressure convergence");
    const double inlet_coarse = inlet_velocity_error(25,false);
    const double inlet_fine = inlet_velocity_error(49,false);
    near(inlet_coarse,10./(24.*24.),1e-12,"positive cubic truncation error for revised pressure inlet");
    near(inlet_fine,10./(48.*48.),1e-12,"halving axial spacing must quarter the inlet error");
    near(inlet_velocity_error(25,true),-4./24. - 5./(24.*24.),1e-12,
         "historical first-order inlet replay");
    const auto persistent = mbla::relax_wall(.5, 1., 1.25, 1., .5);
    require(persistent.used_upstream, "previous and upstream cake carry a below-freeze candidate");
    near(persistent.relaxed, 1.125, 0., "upstream substitution must precede relaxation");
    const auto no_upstream = mbla::relax_wall(.5, 1., .75, 1., .5);
    const auto no_history = mbla::relax_wall(.5, .75, 1.25, 1., .5);
    const auto exact_freeze = mbla::relax_wall(1., 1., 1.25, 1., .5);
    require(!no_upstream.used_upstream && !no_history.used_upstream && !exact_freeze.used_upstream,
            "each history guard operand matters and the trial comparison is strict");
    near(exact_freeze.relaxed, 1., 0., "threshold equality must survive wall relaxation");
    near(mbla::relax_wall(-1., .25, .125, 1., .5).relaxed, .25, 0.,
         "negative relaxed concentration falls back to the previous wall value");
    mbla::IterationStop stopping;
    require(!stopping.observe(.0625, 0., 0., .0625, .25), "wall tolerance is strict");
    require(!stopping.observe(0., .25, 0., .0625, .25), "new outlet cake waits for resistance");
    require(std::isinf(stopping.outlet_resistance_change), "first outlet resistance diagnostic");
    require(!stopping.observe(0., 1.25, 1., .0625, .25), "cake tolerance is strict");
    require(stopping.observe(0., 0., 1., .0625, .25), "cake-free convergence uses wall tolerance");
    near(stopping.outlet_resistance_change, .25, 0., "cake-free iteration retains the prior diagnostic");
    require(!stopping.should_stop(30), "no early break at the minimum index");
    require(!stopping.observe(.125, 0., 0., .0625, .25) && stopping.should_stop(31),
            "a failed subsequent check must not clear the latch");
    mbla::IterationStop settled_cake;
    require(settled_cake.observe(0., 1., 1., .0625, .25), "settled cake meets both tolerances");
    std::cout << "PASS matrix ownership/bounds, LU, CS, dilute closures, cake series law; pressure_errors_Pa="
              << coarse << ',' << fine << "; uniform_cake_errors_Pa=" << cake_coarse << ',' << cake_fine
              << "; velocity_inlet_errors_Pa=" << velocity_coarse << ',' << velocity_fine << '\n';
  } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
