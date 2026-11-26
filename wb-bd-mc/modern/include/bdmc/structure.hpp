#pragma once

#include "bdmc/spatial.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace bdmc {

// Counts are ordered distinct-particle pairs; both +/- displacements are counted.
// Frames have fixed N and box. This is an observer, independent of topology/RNG.
struct PairCorrelationSpec {
    std::vector<double> box_lengths;
    std::vector<double> radial_edges;
    SpatialBackend backend = SpatialBackend::CellList;
    bool finite_population = true;
    // Flow lattice shift on crossing one gradient box, in length units.
    // Nonzero shear currently requires the explicit brute-force backend.
    double shear_displacement = 0.0;
    std::size_t flow_axis = 0;
    std::size_t gradient_axis = 1;
    void validate() const;
};

// Keep raw counts and their expected ideal-gas measure together. The density
// used to normalize finite-N g can be (N-1)/V, while the system density is N/V;
// downstream coordination must use the recorded convention, not silently N/V.
struct PairCorrelationResult {
    std::vector<std::uint64_t> counts;
    std::vector<double> bin_volumes, expected_counts, g;
    std::vector<double> cumulative_coordination;
    std::uint64_t candidate_pairs = 0;
    std::size_t n_frames = 0, n_particles = 0;
    double volume = 0.0, number_density = 0.0, normalization_density = 0.0;
};

PairCorrelationResult pair_correlation(const std::vector<double>& positions,
    std::size_t n_frames, std::size_t n_particles, const PairCorrelationSpec& spec);

// 3D solid-angle sectors: mu=cos(theta) about polar_axis; phi in [0,2*pi].
// azimuth_reference defines phi=0 by its projection normal to polar_axis.
struct DirectionalCorrelationSpec {
    PairCorrelationSpec radial;
    std::vector<double> mu_edges{-1.0, 1.0};
    std::vector<double> phi_edges{0.0, 6.283185307179586476925286766559};
    std::array<double, 3> polar_axis{{0.0, 0.0, 1.0}};
    std::array<double, 3> azimuth_reference{{1.0, 0.0, 0.0}};
    void validate() const;
};

struct DirectionalCorrelationResult {
    PairCorrelationResult histogram;
    std::array<std::size_t, 3> shape{};
    std::uint64_t zero_distance_pairs = 0;
};

DirectionalCorrelationResult directional_pair_correlation(const std::vector<double>& positions,
    std::size_t n_frames, std::size_t n_particles, const DirectionalCorrelationSpec& spec);

// A genuine 3D pair-displacement projection, not a 2D simulation RDF.
// The unplotted coordinate is integrated over the full box by default, or a
// symmetric slab |dr_normal|<=slab_half_width. Bins refer to signed dr_flow,dr_gradient.
struct PlanarCorrelationSpec {
    PairCorrelationSpec geometry;
    std::vector<double> flow_edges, gradient_edges;
    double slab_half_width = -1.0;
    void validate() const;
};

struct PlanarCorrelationResult {
    PairCorrelationResult histogram;
    std::array<std::size_t, 2> shape{};
    double slab_width = 0.0;
};

PlanarCorrelationResult planar_pair_correlation(const std::vector<double>& positions,
    std::size_t n_frames, std::size_t n_particles, const PlanarCorrelationSpec& spec);

enum class CorrelationWindow { AllOrigins, LegacyFixedWindow };
struct AutocorrelationResult {
    std::vector<double> values;
    std::vector<std::uint64_t> origin_counts;
    std::size_t n_components = 0;
};
// Row-major (time,component); result is (lag,component). LegacyFixedWindow
// intentionally retains floor(Nt/2)-1 origins for compatibility.
AutocorrelationResult autocorrelation(const std::vector<double>& values,
    std::size_t n_samples, std::size_t n_components, std::size_t max_lag,
    bool center = false, CorrelationWindow window = CorrelationWindow::AllOrigins);

} // namespace bdmc
