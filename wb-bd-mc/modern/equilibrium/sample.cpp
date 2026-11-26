#include "hs_mc.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
std::uint64_t integer(std::string_view text) {
    std::uint64_t value;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size())
        throw std::invalid_argument("unsigned integer required");
    return value;
}

double number(const std::string& text) {
    std::size_t parsed;
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value))
        throw std::invalid_argument("finite number required");
    return value;
}

// Use more FCC sites than particles and select vacancies independently of the MC stream.
std::vector<double> fcc_start(int particles, double length, std::uint64_t seed) {
    int side = 1;
    while (4 * side * side * side < particles) ++side;
    const double spacing = length / side;
    if (spacing / std::sqrt(2.) < 1)
        throw std::invalid_argument("FCC lattice at these parameters can overlap");
    using Site = std::array<double, 3>;
    constexpr std::array<Site, 4> basis{{{{0, 0, 0}}, {{0, .5, .5}},
                                         {{.5, 0, .5}}, {{.5, .5, 0}}}};
    std::vector<Site> sites;
    for (int i = 0; i < side; ++i)
        for (int j = 0; j < side; ++j)
            for (int k = 0; k < side; ++k)
                for (const auto& shift : basis)
                    sites.push_back({(i + shift[0]) * spacing, (j + shift[1]) * spacing,
                                     (k + shift[2]) * spacing});
    std::mt19937_64 initial_rng(seed ^ 0x97ba31ce85e0b517ULL);
    std::shuffle(sites.begin(), sites.end(), initial_rng);
    std::vector<double> positions;
    positions.reserve(3 * particles);
    for (int i = 0; i < particles; ++i)
        positions.insert(positions.end(), sites[i].begin(), sites[i].end());
    return positions;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 9 && argc != 10)
            throw std::invalid_argument("usage: bdmc_hard_sphere_mc phi nside step seed stream equilibrate production stride [fcc]");
        if (argc == 10 && std::string_view(argv[9]) != "fcc")
            throw std::invalid_argument("unknown initial geometry");
        const double phi = number(argv[1]);
        const double step = number(argv[3]);
        const auto side = integer(argv[2]), seed = integer(argv[4]), stream = integer(argv[5]);
        const auto equilibration = integer(argv[6]), production = integer(argv[7]), stride = integer(argv[8]);
        const double pi = std::acos(-1.);
        if (!(phi > 0 && phi < pi / 6) || side == 0 || side > 16 ||
            stream > UINT64_MAX / 2 || equilibration > 1000000 || production > 1000000 ||
            stride == 0 || production == 0 || production % stride != 0)
            throw std::invalid_argument("invalid simple-cubic or sampling parameters");
        const int n = static_cast<int>(side * side * side);
        const double length = std::pow(n * pi / (6 * phi), 1. / 3.);
        std::vector<double> initial;
        if (argc == 10) {
            initial = fcc_start(n, length, seed);
        } else {
            initial.reserve(3 * n);
            for (unsigned i = 0; i < side; ++i)
                for (unsigned j = 0; j < side; ++j)
                    for (unsigned k = 0; k < side; ++k)
                        for (unsigned coordinate : {i, j, k})
                            initial.push_back((coordinate + 0.5) * length / side);
        }

        // One sweep makes N independent random selections with replacement.
        // Equilibration consumes the same stream and contributes no output rows.
        const std::unique_ptr<HS, decltype(&hs_destroy)> state(
            hs_create(n, length, seed, stream, step, 0, initial.data()), hs_destroy);
        if (!state) throw std::runtime_error("invalid or overlapping initial state, unsupported box or allocation failure");
        hs_sweeps(state.get(), static_cast<int>(equilibration));
        std::cout << std::setprecision(17)
                  << "sweep,min_r2_capped,pairs_0025,pairs_005,pairs_01,pairs_02,pairs_04,mc_msd,acceptance\n";
        for (std::uint64_t done = 0; done < production; done += stride) {
            const auto accepted = hs_sweeps(state.get(), static_cast<int>(stride));
            double observed[7];
            hs_measure(state.get(), observed);
            if (observed[0] < 1) throw std::runtime_error("sample contains overlapping spheres");
            std::cout << equilibration + done + stride;
            for (double value : observed) std::cout << ',' << value;
            std::cout << ',' << static_cast<double>(accepted) / (stride * n) << '\n';
        }
        return std::cout ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
