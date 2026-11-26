#ifndef BDMC_HI_SD_H
#define BDMC_HI_SD_H

#include <vector>

// Deterministic equal-sphere Stokesian Dynamics in an unbounded solvent.
// All quantities use one caller-selected coherent length/time/force system.
// No Brownian noise, periodic images, walls or particle overlaps are supported.
namespace HI_SD {

struct Parameters {
  double radius = 1.0;
  double viscosity = 1.0;
  bool lubrication = true;
  double lubrication_cutoff = 3.99; // centre separation / radius, at most 3.99
};

struct Result {
  std::vector<double> translation; // 3N, includes caller's ambient velocity
  std::vector<double> rotation;    // 3N, includes caller's ambient angular velocity
  std::vector<double> stresslet;   // 9N, symmetric traceless row-major tensors
  double relative_residual = 0.0;
  double minimum_gap_over_radius = 0.0;
};

// q,F,T,Uinf,Omegainf: 3N each; strain: 9 uniform entries, or 9N entries.
// The strain tensors must be symmetric/traceless and use inverse time units.
// This is an instantaneous linear response; the caller owns time integration.
Result solve(const std::vector<double>& positions,
             const std::vector<double>& forces,
             const std::vector<double>& torques,
             const std::vector<double>& background_velocity,
             const std::vector<double>& background_angular_velocity,
             const std::vector<double>& strain,
             const Parameters& parameters = Parameters());

} // namespace HI_SD
#endif
