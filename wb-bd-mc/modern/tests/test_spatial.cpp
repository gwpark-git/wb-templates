#include "bdmc/spatial.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void rejects(Function function, const std::string& message) {
    bool rejected = false;
    try { function(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, message);
}

bdmc::NeighborResult search(const std::vector<double>& positions,
                            const std::vector<double>& lengths, double cutoff,
                            bdmc::SpatialBackend backend = bdmc::SpatialBackend::CellList,
                            const std::vector<bdmc::ParticleId>& ids = {}) {
    return bdmc::SpatialSearch({lengths, cutoff, backend})(positions, ids);
}

void compare(const std::vector<double>& positions,
              const std::vector<double>& lengths, double cutoff,
              const std::vector<bdmc::ParticleId>& ids = {}) {
    const auto cell = search(positions, lengths, cutoff, bdmc::SpatialBackend::CellList, ids);
    const auto brute = search(positions, lengths, cutoff, bdmc::SpatialBackend::BruteForce, ids);
    check(cell.pairs.size() == brute.pairs.size(), "cell/brute pair counts differ");
    const auto n = positions.size() / lengths.size();
    const auto all_pairs = n == 0 ? std::size_t{0} : n * (n - 1) / 2;
    check(cell.stats.candidate_pairs <= all_pairs, "duplicate cell-pair evaluation");
    check(brute.stats.candidate_pairs == all_pairs, "brute-force candidate count");
    check(cell.stats.occupied_cells <= n, "occupied storage exceeds particle count");
    check(cell.stats.stored_particle_indices == n, "one stored index per particle");
    std::pair<std::size_t, std::size_t> previous{0, 0};
    for (std::size_t k = 0; k < cell.pairs.size(); ++k) {
        const auto& a = cell.pairs[k];
        const auto& b = brute.pairs[k];
        check(a.first_index == b.first_index && a.second_index == b.second_index,
              "cell/brute pair indices differ");
        check(a.first_id == b.first_id && a.second_id == b.second_id, "stable pair IDs differ");
        check(a.first_index < a.second_index, "pair must be oriented by storage index");
        const auto current = std::make_pair(a.first_index, a.second_index);
        check(k == 0 || previous < current, "pairs sorted strictly with no duplicates");
        previous = current;
        check(a.distance == b.distance && a.displacement == b.displacement,
              "cell/brute geometry differs");
        check(a.distance <= cutoff, "retained pair exceeds cutoff");
    }
}

void analytic_pairs() {
    const auto result = search({0.0, 0.5, 1.0, 5.0, 9.0}, {10.0}, 1.0);
    const std::vector<std::pair<std::size_t, std::size_t>> expected{{0, 1}, {0, 2}, {0, 4}, {1, 2}};
    check(result.pairs.size() == expected.size(), "analytic periodic pair count");
    for (std::size_t k = 0; k < expected.size(); ++k) {
        check(std::make_pair(result.pairs[k].first_index, result.pairs[k].second_index) == expected[k],
              "analytic periodic pair identity");
    }
    check(result.pairs[1].distance == 1.0, "cutoff equality is included");
    check(result.pairs[2].displacement[0] == -1.0, "minimum-image vector orientation");
    const auto half = search({0.0, 1.0}, {2.0}, 1.0);
    check(half.pairs.size() == 1 && half.pairs[0].displacement[0] == -1.0,
          "half-box tie follows declared convention");
}

void empty_single_and_overlap() {
    for (std::size_t dimension = 1; dimension <= 3; ++dimension) {
        const std::vector<double> lengths(dimension, 2.0);
        compare({}, lengths, 1.0);
        compare(std::vector<double>(dimension, 0.0), lengths, 1.0);
        compare(std::vector<double>(5 * dimension, 0.0), lengths, 1.0);
    }
    const auto coincident = search({0.0, 0.0, 1.0}, {1.0},
                                   std::numeric_limits<double>::denorm_min());
    check(coincident.pairs.size() == 3, "tiny positive cutoff retains coincident wrapped particles");
}

void few_cells_and_boundary_particles() {
    for (std::size_t dimension = 1; dimension <= 3; ++dimension) {
        const std::vector<double> lengths(dimension, 2.0);
        std::vector<double> points;
        const std::vector<double> coordinate{
            -2.0, -std::numeric_limits<double>::epsilon(), 0.0,
            std::nextafter(0.0, 1.0), 0.5, 1.0,
            std::nextafter(2.0, 0.0), 2.0, 4.0};
        std::size_t total = 1;
        for (std::size_t c = 0; c < dimension; ++c) total *= coordinate.size();
        for (std::size_t k = 0; k < total; ++k) {
            auto index = k;
            for (std::size_t c = 0; c < dimension; ++c) {
                points.push_back(coordinate[index % coordinate.size()]);
                index /= coordinate.size();
            }
        }
        compare(points, lengths, 1.0);
        compare(points, lengths, std::nextafter(1.0, 0.0));
        compare(points, lengths, 0.1);

        // Exercise internal bin faces as well as physical box faces. Include
        // immediate floating-point neighbors and exact-cutoff translations.
        const auto shape = search({}, lengths, 0.1).stats.grid_shape;
        std::vector<double> bin_boundary_points;
        for (std::size_t axis = 0; axis < dimension; ++axis) {
            const double width = lengths[axis] / static_cast<double>(shape[axis]);
            for (std::uint64_t cell = 0; cell < shape[axis]; ++cell) {
                const double edge = static_cast<double>(cell) * width;
                for (const double value : {
                        std::nextafter(edge, -std::numeric_limits<double>::infinity()), edge,
                        std::nextafter(edge, std::numeric_limits<double>::infinity()),
                        std::nextafter(edge + 0.1, -std::numeric_limits<double>::infinity()), edge + 0.1,
                        std::nextafter(edge + 0.1, std::numeric_limits<double>::infinity())}) {
                    const auto start = bin_boundary_points.size();
                    bin_boundary_points.resize(start + dimension, 0.0);
                    bin_boundary_points[start + axis] = value;
                }
            }
        }
        compare(bin_boundary_points, lengths, 0.1);
    }
}

void randomized_reference_comparison() {
    std::mt19937_64 generator(19051);
    std::uniform_real_distribution<double> distribution(-3.0, 4.0);
    for (std::size_t dimension = 1; dimension <= 3; ++dimension) {
        std::vector<double> lengths;
        for (std::size_t c = 0; c < dimension; ++c) lengths.push_back(3.0 + c);
        for (std::size_t replicate = 0; replicate < 12; ++replicate) {
            std::vector<double> points(120 * dimension);
            for (std::size_t k = 0; k < points.size(); ++k) points[k] = distribution(generator) * lengths[k % dimension];
            for (const double cutoff : {0.04, 0.3, 0.9, 1.5}) compare(points, lengths, cutoff);
        }
    }
}

void sparse_huge_box() {
    const std::vector<double> lengths{1e12, 2e12, 3e12};
    const std::vector<double> points{
        0.0, 0.0, 0.0, 0.004, 0.002, 0.0,
        1e12 - 0.003, 0.0, 0.0, 5e11, 1e12, 1.5e12,
        -1e12, 4e12, -3e12, 0.0, 0.0, 0.0};
    compare(points, lengths, 0.01);
    const auto result = search(points, lengths, 0.01);
    check(result.stats.occupied_cells <= 6, "huge empty grid stores only occupied cells");
    long double theoretical = 1.0L;
    for (const auto side : result.stats.grid_shape) theoretical *= side;
    check(theoretical > static_cast<long double>(std::numeric_limits<std::uint64_t>::max()),
          "test grid product deliberately exceeds native uint64");
    std::cout << "  sparse grid theoretical cells=" << theoretical
              << " occupied=" << result.stats.occupied_cells
              << " candidates=" << result.stats.candidate_pairs << '\n';
}

void stable_ids_and_permutation() {
    const std::vector<double> points{0.1, 0.2, 1.0, 4.9, 9.9};
    const std::vector<bdmc::ParticleId> ids{91, 4, 72, 6, 1000};
    compare(points, {10.0}, 1.0, ids);
    const auto original = search(points, {10.0}, 1.0, bdmc::SpatialBackend::CellList, ids);
    const std::vector<std::size_t> permutation{3, 0, 4, 1, 2};
    std::vector<double> reordered;
    std::vector<bdmc::ParticleId> reordered_ids;
    for (const auto k : permutation) {
        reordered.push_back(points[k]);
        reordered_ids.push_back(ids[k]);
    }
    const auto changed = search(reordered, {10.0}, 1.0, bdmc::SpatialBackend::CellList, reordered_ids);
    auto keyed = [](const bdmc::NeighborResult& result) {
        std::map<std::pair<bdmc::ParticleId, bdmc::ParticleId>, std::pair<double, double>> map;
        for (const auto& pair : result.pairs) {
            const bool forward = pair.first_id < pair.second_id;
            const auto key = std::minmax(pair.first_id, pair.second_id);
            map[{key.first, key.second}] = {pair.distance, forward ? pair.displacement[0] : -pair.displacement[0]};
        }
        return map;
    };
    check(keyed(original) == keyed(changed), "stable-ID neighbor geometry survives storage permutation");
}

void invalid_search_contracts() {
    for (const double cutoff : {0.0, -1.0, 0.500001,
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::quiet_NaN()}) {
        rejects([&] { search({0.0}, {1.0}, cutoff); }, "invalid cutoff rejected");
    }
    rejects([] { search({}, {}, 1.0); }, "zero-dimensional box rejected");
    rejects([] { search({}, {1, 1, 1, 1}, 0.1); }, "four-dimensional box rejected");
    rejects([] { search({0.0}, {1.0, 1.0}, 0.1); }, "coordinate shape mismatch rejected");
    rejects([] { search({0.0}, {0.0}, 0.1); }, "invalid box length rejected");
    rejects([] { search({std::numeric_limits<double>::infinity()}, {1.0}, 0.1); }, "nonfinite positions rejected");
    rejects([] { search({0.0, 0.1}, {1.0}, 0.1, bdmc::SpatialBackend::CellList, {4, 4}); }, "duplicate IDs rejected");
    rejects([] { search({0.0, 0.1}, {1.0}, 0.1, bdmc::SpatialBackend::CellList, {4}); }, "ID count mismatch rejected");
    rejects([] { bdmc::parse_spatial_backend("automatic"); }, "unsupported backend rejected");
}

void large_geometry_avoids_squared_overflow() {
    const auto result = search({0.0, 0.0, 0.0, 0.5e307, 0.5e307, 0.5e307},
                               {1e308, 1e308, 1e308}, 1e307);
    check(result.pairs.size() == 1, "large finite geometry accepted without squared-distance overflow");
    check(std::isfinite(result.pairs[0].distance), "large pair distance remains finite");
    compare({0.0, 0.0, 0.0, 0.5e307, 0.5e307, 0.5e307}, {1e308, 1e308, 1e308}, 1e307);
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"analytic_pairs", analytic_pairs},
        {"empty_single_and_overlap", empty_single_and_overlap},
        {"few_cells_and_boundary_particles", few_cells_and_boundary_particles},
        {"randomized_reference_comparison", randomized_reference_comparison},
        {"sparse_huge_box", sparse_huge_box},
        {"stable_ids_and_permutation", stable_ids_and_permutation},
        {"invalid_search_contracts", invalid_search_contracts},
        {"large_geometry_avoids_squared_overflow", large_geometry_avoids_squared_overflow},
    };
    int failures = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "PASS " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    std::cout << tests.size() - failures << '/' << tests.size() << " spatial check groups passed\n";
    return failures == 0 ? 0 : 1;
}
