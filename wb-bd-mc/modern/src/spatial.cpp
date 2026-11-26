#include "bdmc/spatial.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace bdmc {
namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::invalid_argument(message);
}

SpatialSearchSpec validated(SpatialSearchSpec spec) {
    spec.validate();
    return spec;
}

struct CellKey {
    std::array<std::uint64_t, 3> value{{0, 0, 0}};
    bool operator==(const CellKey& other) const noexcept { return value == other.value; }
    bool operator<(const CellKey& other) const noexcept { return value < other.value; }
};

struct CellHash {
    std::size_t operator()(const CellKey& key) const noexcept {
        std::size_t hash = 0;
        for (const auto coordinate : key.value) {
            // Unsigned hash arithmetic intentionally wraps; grid coordinates
            // themselves are never flattened into an overflowing product.
            hash ^= std::hash<std::uint64_t>{}(coordinate) +
                    std::size_t{0x9e3779b9U} + (hash << 6U) + (hash >> 2U);
        }
        return hash;
    }
};

std::vector<std::uint64_t> grid_shape(const SpatialSearchSpec& spec) {
    std::vector<std::uint64_t> shape;
    for (const double length : spec.box_lengths) {
        // A conservative floating-point margin makes bins slightly wider than
        // the physical cutoff. It covers arithmetic at cell/PBC boundaries;
        // the actual pair acceptance still uses exactly distance <= cutoff.
        const long double target = static_cast<long double>(spec.cutoff) +
            64.0L * std::numeric_limits<double>::epsilon() * static_cast<long double>(length);
        const long double count = std::floor(static_cast<long double>(length) / target);
        constexpr std::uint64_t max_axis_cells = std::uint64_t{1} << 52U;
        require(std::isfinite(count) && count <= static_cast<long double>(max_axis_cells),
                "cell coordinate count exceeds exact integer grid limit");
        auto cells = count < 2.0L ? std::uint64_t{2} : static_cast<std::uint64_t>(count);
        // Strict slack away from half-box cases; two-cell axes already visit
        // every periodic cell and therefore require no numerical margin.
        if (cells > 2) --cells;
        shape.push_back(cells);
    }
    return shape;
}

CellKey cell_of(const std::vector<double>& wrapped, std::size_t particle,
                 const std::vector<double>& lengths,
                 const std::vector<std::uint64_t>& shape) {
    CellKey key;
    for (std::size_t c = 0; c < lengths.size(); ++c) {
        const long double coordinate = static_cast<long double>(wrapped[particle * lengths.size() + c]);
        const long double fraction = coordinate / static_cast<long double>(lengths[c]);
        const long double bin = std::floor(fraction * static_cast<long double>(shape[c]));
        require(bin >= 0.0L && std::isfinite(bin), "invalid wrapped cell coordinate");
        // The wrapped coordinate is <L; a rounded upper-edge value maps to the
        // final cell instead of overflowing the periodic integer range.
        key.value[c] = bin >= static_cast<long double>(shape[c])
            ? shape[c] - 1 : static_cast<std::uint64_t>(bin);
    }
    return key;
}

std::uint64_t offset_cell(std::uint64_t coordinate, int offset, std::uint64_t count) {
    if (offset == -1) return coordinate == 0 ? count - 1 : coordinate - 1;
    if (offset == 1) return coordinate == count - 1 ? 0 : coordinate + 1;
    return coordinate;
}

std::vector<CellKey> adjacent_cells(const CellKey& key,
                                    const std::vector<std::uint64_t>& shape) {
    std::vector<CellKey> adjacent{key};
    for (std::size_t c = 0; c < shape.size(); ++c) {
        std::vector<CellKey> expanded;
        expanded.reserve(adjacent.size() * 3);
        for (const auto& current : adjacent) {
            for (const int offset : {-1, 0, 1}) {
                auto next = current;
                next.value[c] = offset_cell(key.value[c], offset, shape[c]);
                expanded.push_back(next);
            }
        }
        adjacent.swap(expanded);
    }
    std::sort(adjacent.begin(), adjacent.end());
    // With two cells along an axis, -1 and +1 identify the same periodic cell.
    adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
    return adjacent;
}

void consider_pair(std::size_t i, std::size_t j,
                    const std::vector<double>& wrapped,
                    const std::vector<ParticleId>& ids,
                    const OrthorhombicBox& box, double cutoff,
                    NeighborResult& result) {
    if (j < i) std::swap(i, j);
    if (result.stats.candidate_pairs == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("candidate pair count overflows uint64");
    }
    ++result.stats.candidate_pairs;
    NeighborPair pair;
    pair.first_index = i;
    pair.second_index = j;
    pair.first_id = ids[i];
    pair.second_id = ids[j];
    for (std::size_t c = 0; c < box.dimension(); ++c) {
        pair.displacement[c] = box.minimum_image(
            wrapped[j * box.dimension() + c] - wrapped[i * box.dimension() + c], c);
        // hypot avoids an overflowing squared distance for large finite boxes.
        pair.distance = std::hypot(pair.distance, pair.displacement[c]);
    }
    if (pair.distance <= cutoff) result.pairs.push_back(pair);
}

}  // namespace

std::string spatial_backend_name(SpatialBackend backend) {
    switch (backend) {
        case SpatialBackend::CellList: return "cell_list";
        case SpatialBackend::BruteForce: return "brute_force";
    }
    throw std::invalid_argument("unknown spatial search enum");
}

SpatialBackend parse_spatial_backend(const std::string& name) {
    if (name == "cell_list") return SpatialBackend::CellList;
    if (name == "brute_force") return SpatialBackend::BruteForce;
    throw std::invalid_argument("spatial backend must be 'cell_list' or 'brute_force'");
}

void SpatialSearchSpec::validate() const {
    (void)spatial_backend_name(backend);
    const OrthorhombicBox box(box_lengths);
    require(std::isfinite(cutoff) && cutoff > 0.0, "cutoff must be finite and positive");
    const double shortest = *std::min_element(box_lengths.begin(), box_lengths.end());
    require(cutoff <= 0.5 * shortest, "cutoff must not exceed half the shortest periodic box length");
}

SpatialSearch::SpatialSearch(SpatialSearchSpec spec) : spec_(validated(std::move(spec))) {}

NeighborResult SpatialSearch::operator()(const std::vector<double>& positions,
                                         std::vector<ParticleId> ids) const {
    const OrthorhombicBox box(spec_.box_lengths);
    const std::size_t dimension = box.dimension();
    require(positions.size() % dimension == 0, "spatial positions must have N*dimension values");
    const std::size_t n = positions.size() / dimension;
    require(n <= std::numeric_limits<ParticleId>::max(), "particle count exceeds uint64 IDs");
    if (ids.empty()) {
        ids.resize(n);
        std::iota(ids.begin(), ids.end(), ParticleId{0});
    }
    require(ids.size() == n, "spatial particle ID count mismatch");
    std::unordered_set<ParticleId> unique_ids;
    for (const auto id : ids) require(unique_ids.insert(id).second, "duplicate spatial particle ID");
    std::vector<double> wrapped(positions.size());
    for (std::size_t k = 0; k < positions.size(); ++k) wrapped[k] = box.wrap(positions[k], k % dimension);

    NeighborResult result;
    result.stats.backend = spec_.backend;
    result.stats.n_particles = n;
    result.stats.dimension = dimension;
    if (spec_.backend == SpatialBackend::BruteForce) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) consider_pair(i, j, wrapped, ids, box, spec_.cutoff, result);
        }
    } else {
        result.stats.grid_shape = grid_shape(spec_);
        std::unordered_map<CellKey, std::vector<std::size_t>, CellHash> cells;
        // Reserve by particle count, never by a product of cell counts.
        cells.reserve(n);
        for (std::size_t i = 0; i < n; ++i) cells[cell_of(wrapped, i, spec_.box_lengths, result.stats.grid_shape)].push_back(i);
        result.stats.occupied_cells = cells.size();
        result.stats.stored_particle_indices = n;
        for (const auto& entry : cells) {
            for (const auto& neighbor : adjacent_cells(entry.first, result.stats.grid_shape)) {
                if (neighbor < entry.first) continue;  // each unordered cell pair once
                const auto found = cells.find(neighbor);
                if (found == cells.end()) continue;
                if (neighbor == entry.first) {
                    for (std::size_t a = 0; a < entry.second.size(); ++a) {
                        for (std::size_t b = a + 1; b < entry.second.size(); ++b) {
                            consider_pair(entry.second[a], entry.second[b], wrapped, ids, box, spec_.cutoff, result);
                        }
                    }
                } else {
                    for (const auto i : entry.second) {
                        for (const auto j : found->second) consider_pair(i, j, wrapped, ids, box, spec_.cutoff, result);
                    }
                }
            }
        }
        std::sort(result.pairs.begin(), result.pairs.end(), [](const NeighborPair& a, const NeighborPair& b) {
            // Hash-table iteration is not a physical order. Return canonical
            // pairs so a cache implementation does not silently reorder a
            // downstream force reduction; IDs remain separate from these slots.
            return std::make_pair(a.first_index, a.second_index) < std::make_pair(b.first_index, b.second_index);
        });
    }
    return result;
}

}  // namespace bdmc
