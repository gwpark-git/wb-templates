#pragma once

#include <cmath>

namespace mbla {

// The wall values in this control are dimensionless ratios phi_w/phi_b.
struct WallUpdate {
  double relaxed;
  bool used_upstream;
};

// The current upstream node has already been updated in this axial sweep.
inline WallUpdate relax_wall(double candidate, double previous, double upstream,
                             double freeze_ratio, double weight) {
  const bool used_upstream = candidate < freeze_ratio &&
                             previous >= freeze_ratio && upstream >= freeze_ratio;
  if (used_upstream) {
    candidate = upstream;
  }
  const double relaxed = (1. - weight) * previous + weight * candidate;
  return {relaxed < 0. ? previous : relaxed, used_upstream};
}

// The criterion can latch before the earliest allowed break at iteration 31.
struct IterationStop {
  bool latched = false;
  double outlet_resistance_change = 0.;

  // Preserve the latest positive-cake diagnostic, including its IEEE
  // infinity when positive resistance first appears after zero resistance.
  bool observe(double wall_error, double outlet_resistance, double previous_resistance,
               double wall_tolerance, double resistance_tolerance) {
    bool criteria_met;
    if (outlet_resistance > 0.) {
      outlet_resistance_change = std::fabs(1. - outlet_resistance / previous_resistance);
      criteria_met = wall_error < wall_tolerance &&
                     outlet_resistance_change < resistance_tolerance;
    } else {
      criteria_met = wall_error < wall_tolerance;
    }
    latched = latched || criteria_met;
    return criteria_met;
  }

  bool should_stop(long iteration) const {
    return latched && iteration > 30;
  }
};

} // namespace mbla
