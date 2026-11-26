#pragma once

#include "bdmc/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace stiff {

inline constexpr double cosine_pi = 3.14159265358979323846;
inline constexpr double cosine_amplitude = 100.;

// A cosine force term for the Level 3 BD-MC core. Lengths are in diameters
// and energies are in kBT. BD-MC owns the Gaussian noise, clock and step.
class CosinePairs {
public:
    explicit CosinePairs(double length) : length_(length), half_(0.5 * length) {
        if (!std::isfinite(length) || 2 * neighbor_cutoff >= length)
            throw std::invalid_argument("cosine neighbor radius must be inside half the box");
    }

    void evaluate(const bdmc::ParticleState& state, bdmc::ForceResult& result,
                  std::array<double, 6>* virial = nullptr) {
        refresh(state);
        if (virial) virial->fill(0.);
        for (const auto& [i, j] : pairs_) {
            const auto r = separation(state.positions, i, j);
            const double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            if (r2 >= 1.) continue;
            const double theta = cosine_pi * std::sqrt(r2);
            const double factor = cosine_amplitude * cosine_pi * cosine_pi *
                                  (theta == 0. ? 1. : std::sin(theta)/theta);
            result.potential_energy += cosine_amplitude * (1. + std::cos(theta));
            for (std::size_t axis = 0; axis < r.size(); ++axis) {
                const double force = factor * r[axis];
                result.values[i+axis] -= force;
                result.values[j+axis] += force;
            }
            if (virial) {
                // Positive pair virial, each pair counted once.
                (*virial)[0] += factor*r[0]*r[0];
                (*virial)[1] += factor*r[1]*r[1];
                (*virial)[2] += factor*r[2]*r[2];
                (*virial)[3] += factor*r[0]*r[1];
                (*virial)[4] += factor*r[0]*r[2];
                (*virial)[5] += factor*r[1]*r[2];
            }
        }
    }

    // The accelerator reuses the same double-precision neighbor list. Its
    // device representation needs updating only when the Verlet skin expires.
    bool refresh(const bdmc::ParticleState& state) {
        if (!needs_rebuild(state.unwrapped_positions)) return false;
        rebuild(state);
        return true;
    }

    const std::vector<std::pair<std::size_t, std::size_t>>& pairs() const { return pairs_; }

private:
    static constexpr double skin = 0.3;
    static constexpr double neighbor_cutoff = 1. + skin;

    std::array<double, 3> separation(const std::vector<double>& positions,
                                     std::size_t i, std::size_t j) const {
        std::array<double, 3> r{};
        for (std::size_t axis = 0; axis < r.size(); ++axis) {
            r[axis] = positions[j+axis] - positions[i+axis];
            // With wrapped inputs one correction implements BD-MC's minimum image.
            if (r[axis] >= half_) r[axis] -= length_;
            if (r[axis] < -half_) r[axis] += length_;
        }
        return r;
    }

    bool needs_rebuild(const std::vector<double>& unwrapped) const {
        if (reference_.size() != unwrapped.size()) return true;
        constexpr double limit = (skin/2) * (skin/2);
        for (std::size_t i = 0; i < unwrapped.size(); i += 3) {
            const double dx = unwrapped[i] - reference_[i];
            const double dy = unwrapped[i+1] - reference_[i+1];
            const double dz = unwrapped[i+2] - reference_[i+2];
            if (dx*dx + dy*dy + dz*dz >= limit) return true;
        }
        return false;
    }

    void rebuild(const bdmc::ParticleState& state) {
        // Each particle has moved less than half the skin between rebuilds.
        // A pair initially outside the list cannot cross the force cutoff.
        if (reference_.size() != state.unwrapped_positions.size())
            prepare_cells(state.positions.size()/3);
        reference_ = state.unwrapped_positions;
        pairs_.clear();
        if (!cells_.empty()) {
            rebuild_cells(state.positions);
            return;
        }
        for (std::size_t i = 0; i < state.positions.size(); i += 3) {
            for (std::size_t j = i+3; j < state.positions.size(); j += 3) {
                add_neighbor(state.positions, i, j);
            }
        }
    }

    void prepare_cells(std::size_t count) {
        // Cap occupancy metadata for very dilute systems. Slightly widen
        // the minimum bin width to cover floating-point cell-edge rounding.
        const double target_width = neighbor_cutoff + 32*std::numeric_limits<double>::epsilon()*length_;
        const auto n = static_cast<std::size_t>(std::min(length_/target_width, std::cbrt(double(count))));
        cells_per_axis_ = n;
        cell_scale_ = n/length_;
        // Three or fewer bins per axis already visit every bin.
        cells_.resize(n >= 4 ? n*n*n : 0);
    }

    void add_neighbor(const std::vector<double>& positions, std::size_t i, std::size_t j) {
        const auto r = separation(positions, i, j);
        if (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] < neighbor_cutoff*neighbor_cutoff)
            pairs_.emplace_back(i, j);
    }

    void rebuild_cells(const std::vector<double>& positions) {
        const auto n = cells_per_axis_;
        const auto cell_id = [n](std::size_t x, std::size_t y, std::size_t z) {
            return (x*n + y)*n + z;
        };
        for (auto& cell : cells_) cell.clear();
        particle_cells_.resize(positions.size()/3);
        for (std::size_t i = 0; i < positions.size(); i += 3) {
            auto& index = particle_cells_[i/3];
            for (std::size_t a = 0; a < index.size(); ++a)
                index[a] = std::min(n-1, static_cast<std::size_t>(positions[i+a]*cell_scale_));
            cells_[cell_id(index[0], index[1], index[2])].push_back(i);
        }
        const auto neighbors = [n](std::size_t index) {
            return std::array<std::size_t, 3>{index == 0 ? n-1 : index-1, index,
                                              index+1 == n ? 0 : index+1};
        };
        // Bin width is at least the list radius. Periodic nearest bins span
        // all candidates. Preserve the original i,j force accumulation order.
        for (std::size_t i = 0; i < positions.size(); i += 3) {
            const auto& index = particle_cells_[i/3];
            const auto first = pairs_.size();
            for (auto x : neighbors(index[0]))
                for (auto y : neighbors(index[1]))
                    for (auto z : neighbors(index[2]))
                        for (auto j : cells_[cell_id(x, y, z)])
                            if (j > i) add_neighbor(positions, i, j);
            std::sort(pairs_.begin()+first, pairs_.end());
        }
    }

    const double length_;
    const double half_;
    std::size_t cells_per_axis_{};
    double cell_scale_{};
    std::vector<std::vector<std::size_t>> cells_;
    std::vector<std::array<std::size_t, 3>> particle_cells_;
    std::vector<double> reference_;
    std::vector<std::pair<std::size_t, std::size_t>> pairs_;
};

} // namespace stiff
