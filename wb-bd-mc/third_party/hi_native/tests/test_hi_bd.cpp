#include "../lib/hi_bd.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
const double pi = 3.141592653589793238462643383279502884;
int checks_run = 0;
bool passed = true;

void check(const char* name, bool ok) {
    if (checks_run++) std::cout << ",\n";
    std::cout << "    \"" << name << "\": " << (ok ? "true" : "false");
    passed = passed && ok;
}

bool near(double a, double b, double tolerance = 2e-13) {
    return std::abs(a - b) <= tolerance;
}

bool throws(const std::function<void()>& f) {
    try { f(); } catch (const std::exception&) { return true; }
    return false;
}

double divergence(const std::vector<double>& q, const HI_BD::Parameters& parameters) {
    std::vector<double> div(q.size(), 0.0);
    const double h = 1e-4;
    for (std::size_t j = 0; j < q.size(); ++j) {
        std::vector<double> plus(q), minus(q);
        plus[j] += h;
        minus[j] -= h;
        const std::vector<double> mp = HI_BD::mobility(plus, parameters);
        const std::vector<double> mm = HI_BD::mobility(minus, parameters);
        for (std::size_t i = 0; i < q.size(); ++i)
            div[i] += (mp[i * q.size() + j] - mm[i * q.size() + j]) / (2.0 * h);
    }
    double maximum = 0.0;
    for (std::size_t i = 0; i < div.size(); ++i)
        maximum = std::max(maximum, std::abs(div[i]));
    return maximum;
}
} // namespace

int main() {
    std::cout << std::setprecision(17) << "{\n  \"checks\": {\n";
    const HI_BD::Parameters rpy = {1.0, 1.0 / (6.0 * pi), HI_BD::Model::RPY};
    const HI_BD::Parameters oseen = {1.0, 1.0 / (6.0 * pi), HI_BD::Model::Oseen};
    const std::vector<double> pair = {0., 0., 0., 4., 0., 0.};
    const std::vector<double> zeros(6, 0.0);
    const std::vector<double> mr = HI_BD::mobility(pair, rpy);
    const std::vector<double> mo = HI_BD::mobility(pair, oseen);
    check("rpy_pair_at_4a_literal", near(mr[3], .359375)
          && near(mr[1 * 6 + 4], .1953125) && near(mr[2 * 6 + 5], .1953125));
    check("oseen_pair_at_4a_literal", near(mo[3], .375)
          && near(mo[1 * 6 + 4], .1875) && near(mo[2 * 6 + 5], .1875));
    const std::vector<double> overlap = {0., 0., 0., 1., 0., 0.};
    const std::vector<double> m_overlap = HI_BD::mobility(overlap, rpy);
    check("rpy_overlap_at_a_literal", near(m_overlap[3], .8125)
          && near(m_overlap[1 * 6 + 4], .71875));

    const HI_BD::Parameters dimensional = {2.0, 3.0, HI_BD::Model::RPY};
    const std::vector<double> single = HI_BD::mobility({0., 0., 0.}, dimensional);
    check("single_sphere_stokes_drag", near(single[0], 1.0 / (36.0 * pi))
          && near(single[4], single[0]) && near(single[8], single[0])
          && single[1] == 0.0 && single[2] == 0.0);

    const double s = std::sqrt(8.0);
    const std::vector<double> rotated = HI_BD::mobility({0.,0.,0.,s,s,0.}, rpy);
    bool symmetric = true;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            symmetric = symmetric && near(rotated[i * 6 + j], rotated[j * 6 + i]);
    check("rotation_and_lorentz_reciprocity", symmetric
          && near(rotated[3], .27734375) && near(rotated[4], .08203125)
          && near(rotated[1 * 6 + 4], .27734375));

    const double div_rpy = std::max(divergence(pair, rpy), divergence(overlap, rpy));
    const double div_oseen = divergence(pair, oseen);
    check("free_space_configuration_divergence", std::max(div_rpy, div_oseen) < 1e-8);

    std::vector<double> force(6, 0.0), background(6, 0.0);
    force[0] = 1.0;
    background[3] = .2;
    const std::vector<double> deterministic = HI_BD::increment(pair, force, background,
                                                               .01, 0.0, 0, rpy);
    check("deterministic_force_and_background", near(deterministic[0], .01)
          && near(deterministic[3], .00559375));

    // A known-zero-mean second moment at one fixed configuration. The expected
    // covariance uses literal pair references, not another call to mobility().
    const std::size_t samples = 160000;
    const double dt = .002, energy = 1.7;
    const double c0 = 2 * energy * dt;
    std::vector<double> expected(36, 0.0), sums(6, 0.0), moments(36, 0.0);
    for (std::size_t i = 0; i < 6; ++i) expected[i * 6 + i] = c0;
    expected[3] = expected[3 * 6] = c0 * .359375;
    expected[1 * 6 + 4] = expected[4 * 6 + 1] = c0 * .1953125;
    expected[2 * 6 + 5] = expected[5 * 6 + 2] = c0 * .1953125;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 810729UL);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const std::vector<double> dq = HI_BD::increment(pair, zeros, zeros, dt, energy, rng, rpy);
        for (std::size_t i = 0; i < 6; ++i) {
            sums[i] += dq[i];
            for (std::size_t j = 0; j < 6; ++j) moments[i * 6 + j] += dq[i] * dq[j];
        }
    }
    double covariance_z = 0.0, mean_z = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        mean_z = std::max(mean_z, std::abs(sums[i] / samples) / std::sqrt(c0 / samples));
        for (std::size_t j = 0; j < 6; ++j) {
            const double target = expected[i * 6 + j];
            const double se = std::sqrt((c0 * c0 + target * target) / samples);
            covariance_z = std::max(covariance_z,
                std::abs(moments[i * 6 + j] / samples - target) / se);
        }
    }
    check("correlated_thermal_covariance_6SE", covariance_z < 6.0);
    check("zero_thermal_noise_mean_6SE", mean_z < 6.0);

    // Validate covariances for the distinct Oseen approximation as well.
    gsl_rng_set(rng, 3519UL);
    const std::size_t oseen_samples = 80000;
    double ox = 0.0;
    for (std::size_t sample = 0; sample < oseen_samples; ++sample) {
        const std::vector<double> dq = HI_BD::increment(pair, zeros, zeros, dt, energy, rng, oseen);
        ox += dq[0] * dq[3];
    }
    const double oseen_target = c0 * .375;
    const double oseen_z = std::abs(ox / oseen_samples - oseen_target)
        / std::sqrt((c0 * c0 + oseen_target * oseen_target) / oseen_samples);
    check("oseen_cross_particle_noise_6SE", oseen_z < 6.0);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const HI_BD::Parameters bad_radius = {-1., 1., HI_BD::Model::RPY};
    bool invalid = throws([&] { HI_BD::mobility({0.,0.}, rpy); })
        && throws([&] { HI_BD::mobility({nan,0.,0.}, rpy); })
        && throws([&] { HI_BD::mobility(pair, bad_radius); })
        && throws([&] { HI_BD::increment(pair, zeros, zeros, 0., 1., rng, rpy); })
        && throws([&] { HI_BD::increment(pair, zeros, zeros, dt, -1., rng, rpy); })
        && throws([&] { HI_BD::increment(pair, zeros, zeros, dt, 1., 0, rpy); })
        && throws([&] { HI_BD::increment(pair, {}, zeros, dt, 0., 0, rpy); })
        && throws([&] { HI_BD::mobility(overlap, oseen); });
    check("invalid_inputs_and_oseen_overlap_rejected", invalid);
    check("coincident_rpy_singular_matrix_rejected", throws([&] {
        HI_BD::increment({0.,0.,0.,0.,0.,0.}, zeros, zeros, dt, 0., 0, rpy);
    }));
    gsl_rng_free(rng);
    std::cout << "\n  },\n  \"status\": \"" << (passed ? "passed" : "failed") << "\",\n"
              << "  \"check_count\": " << checks_run << ",\n"
              << "  \"RPY_covariance_samples\": " << samples << ",\n"
              << "  \"RPY_maximum_covariance_z\": " << covariance_z << ",\n"
              << "  \"RPY_maximum_mean_z\": " << mean_z << ",\n"
              << "  \"Oseen_cross_covariance_z\": " << oseen_z << ",\n"
              << "  \"free_space_divergence_maximum\": " << std::max(div_rpy, div_oseen) << ",\n"
              << "  \"scope\": \"Native free-space translational small-N HI; conditional short-time tests, not long-time suspension or wall validation\"\n}\n";
    return passed ? 0 : 1;
}
