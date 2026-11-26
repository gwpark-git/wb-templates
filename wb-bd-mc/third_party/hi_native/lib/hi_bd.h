#ifndef BD_MC_HI_BD_H
#define BD_MC_HI_BD_H

#include <gsl/gsl_rng.h>
#include <vector>

// Native, equal-radius, unbounded three-dimensional translational Brownian HI.
// Units are caller-selected but must be consistent across radius, viscosity,
// force, energy, time and imposed background velocity. No periodic images,
// stresslet constraints, lubrication or wall Green functions are included.
namespace HI_BD {

enum class Model { Oseen, RPY };

struct Parameters {
    double radius;
    double viscosity;
    Model model;
};

// q is [x0,y0,z0,x1,y1,z1,...]. Returns a row-major (3N)x(3N) mobility.
// The self block is I/(6*pi*viscosity*radius). Oseen requires r>=2*radius;
// RPY uses the regularized overlapping-sphere branch. Matrix assembly does
// not certify positive definiteness; increment() factors and checks it.
std::vector<double> mobility(const std::vector<double>& q,
                             const Parameters& parameters);

// Euler--Maruyama dq = (u + M F)dt + sqrt(2*kBT*dt) L xi, L L^T=M.
// All three vectors have length 3N. For this equal-radius free-space mobility,
// div_q M=0, so there is no additional Ito thermal-drift term. This assumption
// must not be reused automatically for wall/constraint/periodic approximations.
// The matrix is checked even at kBT=0. rng may be null only when kBT=0.
// Nonpositive/numerically singular pivots raise std::domain_error; there is
// no matrix clipping, diagonal jitter or silently changed mobility.
std::vector<double> increment(const std::vector<double>& q,
                              const std::vector<double>& force,
                              const std::vector<double>& background,
                              double dt, double kBT, gsl_rng* rng,
                              const Parameters& parameters);

} // namespace HI_BD
#endif
