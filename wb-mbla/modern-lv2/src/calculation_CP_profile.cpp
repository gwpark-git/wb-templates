#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "case_driver.h"
#include "fixed_point_controls.h"
#include "ported_lib.h"

namespace {

// Report early transients more often, keeping the historical iteration indices.
bool should_report_progress(long iteration) {
  return (iteration <= 50 && iteration % 10 == 0) ||
         (iteration <= 100 && iteration % 20 == 0) ||
         (iteration <= 500 && iteration % 50 == 0) ||
         iteration % 100 == 0;
}

// Capture the iteration's profiles before final analysis overwrites its work arrays.
void report_details(const std::string& output_path, AUX_MATRIX& flow, AUX_MATRIX& profiles,
                    DISP_TYPE& dispersion, PARAMETER_SET& parameters, MATRIX& axial_grid,
                    MATRIX& radial_grid) {
  const long outlet = parameters.Nz - 1;
  for (long j = 0; j < parameters.Ny; ++j) {
    const double phi = profiles.phi_yt_div_phib_arr(outlet, j) * parameters.phi_bulk;
    const double shear_rate = flow.shear_rate_nodim(outlet, j);
    const double shear_stress = flow.shear_stress_nodim(outlet, j);
    // The inherited .Rout coordinate is y/R². The reader multiplies it by R.
    flow.report_outlet_profile(j, 0) = radial_grid(j) / parameters.R_channel;
    flow.report_outlet_profile(j, 1) = phi;
    flow.report_outlet_profile(j, 2) = shear_rate;
    flow.report_outlet_profile(j, 3) = shear_stress;
    flow.report_outlet_profile(j, 4) = dispersion.fcn_Dc_div_D0(phi, parameters, shear_rate);
    flow.report_outlet_profile(j, 5) = dispersion.fcn_eta_div_eta0(phi, parameters, shear_stress);
  }
  for (long i = 0; i < parameters.Nz; ++i) {
    flow.F_info_arr(i, 0) = axial_grid(i);
  }

  const auto outlet_path = output_path + ".Rout";
  std::ofstream outlet_report(outlet_path);
  if (!outlet_report) {
    throw std::runtime_error("cannot create outlet report");
  }
  flow.report_outlet_profile.fprint_basic(outlet_report);
  finish_output(outlet_report, outlet_path);

  const auto integral_path = output_path + ".Fval";
  std::ofstream integral_report(integral_path);
  if (!integral_report) {
    throw std::runtime_error("cannot create integral report");
  }
  flow.F_info_arr.fprint_basic(integral_report);
  finish_output(integral_report, integral_path);
}

// Report the current wall RMS and the saved outlet resistance under their legacy labels.
void report_progress(long iteration, long peak, long cake_onset, double wall_error,
                     double resistance_change, MATRIX& wall, MATRIX& axial_grid,
                     MATRIX& previous_resistance, AUX_MATRIX& flow, PARAMETER_SET& parameters) {
  const long outlet = parameters.Nz - 1;
  std::cout << "\nn_ITER: " << iteration << std::endl;
  std::cout << "cond_GT: Pin=" << parameters.Pin << ", Pout=" << parameters.Pout << std::endl;
  std::printf("n=%ld\n[CP] phiw(0)=%4.4f, phiw(L)=%4.4f, max:(phiw(%4.3f))=%4.4f\n",
              iteration, wall(0) * parameters.phi_bulk, wall(outlet) * parameters.phi_bulk,
              axial_grid(peak), wall(peak) * parameters.phi_bulk);
  std::printf("norm(1. - p2/p1) : %4.3e, weight : %4.3e, weight_f : %4.3e\n",
              wall_error, parameters.weight, parameters.weight_f);
  if (parameters.str_Dc_div_D0 == "SID") {
    const double outlet_phi = wall(outlet) * parameters.phi_bulk;
    const double shear_rate = flow.shear_rate_nodim(outlet, 0);
    const double sid = SID::Dc_short_div_D0(outlet_phi, parameters, shear_rate);
    const double hs = SPHS::Dc_short_div_D0(outlet_phi, parameters);
    std::printf("[SID] at outlet wall: (actual) Pe_a_wL = %4.3e (Pa_a_ast = %4.3e, Pe_a_wL/Pe_a_ast = %4.3e), D_SID = %4.3e, D_HS = %4.3e, D_SID/D_HS = %4.3e\n",
                shear_rate, parameters.Pe_a, shear_rate / parameters.Pe_a, sid, hs, sid / hs);
    if (parameters.str_eta_div_eta0 == "TEST") {
      const double shear_stress = flow.shear_stress_nodim(outlet, 0);
      const double sid_viscosity = SID::eta_div_eta0_test(outlet_phi, parameters, shear_stress);
      const double equilibrium_viscosity = SID::eta_div_eta0_test(outlet_phi, parameters, 0.);
      std::printf("                      tau_wL/tau_ast  = %4.3e (tau_ast = %4.3e Pa), etar_SID = %4.3e, etar_SID/etar_EQ = %4.3e\n",
                  shear_stress, parameters.tau_ast, sid_viscosity,
                  sid_viscosity / equilibrium_viscosity);
    }
  }
  if (cake_onset <= parameters.Nz) {
    std::printf("[CAKE] (z/L)[%ld] = %0.3f\n", cake_onset, axial_grid(cake_onset));
    std::printf("\tind=%ld, Rcake(L)/Rmem=%4.3e, abs(1-dR_2/dR_1)=%4.3e\n",
                cake_onset, previous_resistance(outlet), resistance_change);
    std::printf("\td_cake/R=%4.3e, uL_div_u0=%4.3e\n",
                parameters.CAKE_convert_unit * previous_resistance(outlet), flow.uz_div_u0_arr(outlet));
  }
}

// Seed the wall from column 2, in volume fraction. The inlet keeps the imposed bulk value.
void load_continuation(COND& condition, PARAMETER_SET& parameters, MATRIX& wall) {
  if (condition("CONT") != "TRUE") {
    return;
  }
  MATRIX seed(condition("FN_CONT").c_str());
  if (seed.rows != parameters.Nz || seed.cols < 2) {
    throw std::runtime_error("continuation table requires Nz rows and at least two columns");
  }
  for (long i = 1; i < parameters.Nz; ++i) {
    wall(i) = seed(i, 1) / parameters.phi_bulk;
    std::cout << i << '\t' << seed(i, 1) << std::endl;
  }
}

// Own one case's configuration, closure selection, work arrays and output lifecycle.
int run_case(const char* input_path, mbla::EquilibriumTransport* transport) {
  COND condition(input_path);
  condition.cond_print();
  PARAMETER_SET parameters;
  parameters.init(condition);
  parameters.print();
  DISP_TYPE dispersion(parameters, transport);

  MATRIX axial_grid(parameters.Nz, 1, 0.);
  gen_z_div_L_arr(axial_grid, parameters);
  MATRIX previous_wall(parameters.Nz, 1, parameters.phiw_init / parameters.phi_bulk);
  previous_wall(0) = parameters.phi_bulk / parameters.phi_bulk;
  load_continuation(condition, parameters, previous_wall);
  MATRIX wall(previous_wall);

  MATRIX previous_resistance(parameters.Nz, 1, 0.);
  MATRIX resistance(parameters.Nz, 1, 0.);
  MATRIX osmotic_pressure(parameters.Nz, 1, 0.);
  MATRIX radial_grid;
  gen_yt_arr(radial_grid, parameters);
  AUX_MATRIX flow(parameters), profiles(parameters);
  shear_rate_initialization(flow.shear_rate_nodim, radial_grid, parameters);
  shear_stress_initialization(flow.shear_stress_nodim, flow.shear_rate_nodim,
                              radial_grid, parameters, dispersion);

  const long outlet = parameters.Nz - 1;
  long cake_onset = parameters.Nz + 1;
  long last_reported_iteration = 0;
  mbla::IterationStop stopping;
  for (long iteration = 0; iteration < parameters.N_iter; ++iteration) {
    const long previous_cake_onset = cake_onset;
    // Osmotic pressure and the hydraulic solve use the saved iteration.
    // The axial wall/cake sweep then updates stations in increasing z order.
    previous_wall.copy_obj(wall);
    previous_resistance.copy_obj(resistance);
    flow.shear_rate_nodim_pre.copy_obj(flow.shear_rate_nodim);
    flow.shear_stress_nodim_pre.copy_obj(flow.shear_stress_nodim);
    update_Pi_div_DLP_arr(osmotic_pressure, previous_wall, flow.shear_rate_nodim_pre,
                          dispersion, parameters);
    cake_onset = phiw_update(wall, resistance, axial_grid, previous_wall, previous_cake_onset,
                             previous_resistance, osmotic_pressure, radial_grid, iteration,
                             dispersion, parameters, flow, profiles);
    for (long i = 0; i < parameters.Nz; ++i) {
      if (!std::isfinite(wall(i)) || !std::isfinite(resistance(i))) {
        throw std::runtime_error("nonfinite fixed-point iterate; refine/review the numerical case");
      }
    }
    const long peak = wall.argmax_1d();
    const double wall_error = ABS_ERR_data(wall, previous_wall);
    const bool has_outlet_cake = resistance(outlet) > 0.;
    if (stopping.observe(wall_error, resistance(outlet), previous_resistance(outlet),
                         parameters.TOL_phiw, parameters.TOL_delta_f)) {
      std::printf(has_outlet_cake
                    ? "the iteration is converaged (with cake) at iter=%04ld\n"
                    : "the iteration is converaged (w/o cake) at iter=%04ld\n",
                  iteration);
    }

    if (stopping.should_stop(iteration)) {
      if (parameters.DETAILED_REPORT) {
        report_details(condition("fout"), flow, profiles, dispersion, parameters,
                       axial_grid, radial_grid);
      }
      break;
    }
    if (should_report_progress(iteration)) {
      report_progress(iteration, peak, cake_onset, wall_error, stopping.outlet_resistance_change,
                      wall, axial_grid, previous_resistance, flow, parameters);
    }
    last_reported_iteration = iteration;
  }

  // Sidecars describe the iteration above. The main table is recomputed here.
  std::printf("\nFINAL CHECK POINT (total iteration number = %ld):\n", last_reported_iteration);
  std::printf("COUNT force replace phiw : %ld\n", parameters.count_negative_phiw);
  fprint_analysis_with_fouling(axial_grid, wall, resistance, radial_grid, dispersion,
                               parameters, condition("fout").c_str(), flow);
  std::cout << "MBLA_STATUS=" << (stopping.latched ? "converged" : "iteration_limit") << std::endl;
  return stopping.latched ? 0 : 3;
}

} // namespace

// Dispatch positional inputs in order so later cases can use earlier outputs as seeds.
int mbla::run_cases(int argc, char* argv[], EquilibriumTransport* transport) {
  gsl_set_error_handler_off();
  if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--help")) {
    std::cout << "Usage: mbla case.inp [next_case.inp ...]\n"
                 "See docs/SCIENTIFIC_GUIDE.org for SI units, output columns and validation limits.\n";
    return 0;
  }
  try {
    // A continuation in a later input can use output from an earlier case.
    for (int i = 1; i < argc; ++i) {
      const int status = run_case(argv[i], transport);
      if (status != 0) {
        return status;
      }
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "mbla: " << error.what() << '\n';
    return 2;
  }
}
