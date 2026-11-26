#pragma once

#include "bdmc/core.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace bdmc {

enum class SpatialBackend { CellList, BruteForce };
std::string spatial_backend_name(SpatialBackend backend);
SpatialBackend parse_spatial_backend(const std::string& name);

struct SpatialSearchSpec {
    std::vector<double> box_lengths;
    double cutoff = 0.0;
    SpatialBackend backend = SpatialBackend::CellList;
    void validate() const;
};

// Storage indices address today's coordinate buffer; IDs carry identity to
// downstream physics. A cell/bin index is never a replacement for a particle
// or a chain-end label when a neighbor structure is rebuilt.
struct NeighborPair {
    std::size_t first_index = 0;
    std::size_t second_index = 0;
    ParticleId first_id = 0;
    ParticleId second_id = 0;
    double distance = 0.0;
    // Vector from first_index to second_index; unused coordinates are zero.
    std::array<double, 3> displacement{{0.0, 0.0, 0.0}};
};

struct SpatialSearchStats {
    SpatialBackend backend = SpatialBackend::CellList;
    std::size_t n_particles = 0;
    std::size_t dimension = 0;
    std::size_t occupied_cells = 0;
    std::vector<std::uint64_t> grid_shape;
    std::uint64_t candidate_pairs = 0;
    std::size_t stored_particle_indices = 0;
};

struct NeighborResult {
    std::vector<NeighborPair> pairs;
    SpatialSearchStats stats;
};

// Independent spatial candidates, not permanent/dynamic physical bonds.
// Storage is sparse in occupied cells; accepted pairs are retained, candidate
// pairs are evaluated as a stream. The brute-force route does not use bins.
class SpatialSearch {
public:
    explicit SpatialSearch(SpatialSearchSpec spec);
    const SpatialSearchSpec& spec() const noexcept { return spec_; }
    NeighborResult operator()(const std::vector<double>& positions,
                              std::vector<ParticleId> particle_ids = {}) const;

private:
    const SpatialSearchSpec spec_;
};

}  // namespace bdmc
