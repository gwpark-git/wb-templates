#include "hi_bd.h"

#include <gsl/gsl_randist.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

const double pi = 3.141592653589793238462643383279502884;

void finite_vector(const std::vector<double>& values, const char* label) {
    for (std::size_t i = 0; i < values.size(); ++i)
        if (!std::isfinite(values[i]))
            throw std::invalid_argument(std::string(label) + " must be finite");
}

// A small dense reference factorization avoids changing GSL's process-global
// error handler. Its tolerance rejects unresolved singular matrices rather
// than modifying the requested physical mobility.
std::vector<double> cholesky(const std::vector<double>& matrix, std::size_t n) {
    std::vector<double> lower(n * n, 0.0);
    double max_diagonal = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        max_diagonal = std::max(max_diagonal, matrix[i * n + i]);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                           * max_diagonal;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double entry = matrix[i * n + j];
            for (std::size_t k = 0; k < j; ++k)
                entry -= lower[i * n + k] * lower[j * n + k];
            if (i == j) {
                if (!std::isfinite(entry) || !(entry > tolerance))
                    throw std::domain_error(
                        "HI mobility is not numerically positive definite; "
                        "no regularisation was applied");
                lower[i * n + i] = std::sqrt(entry);
            } else {
                lower[i * n + j] = entry / lower[j * n + j];
                if (!std::isfinite(lower[i * n + j]))
                    throw std::domain_error("Nonfinite HI Cholesky factor");
            }
        }
    }
    return lower;
}

} // namespace

std::vector<double> HI_BD::mobility(const std::vector<double>& q,
                                   const Parameters& parameters) {
    if (q.empty() || q.size() % 3 != 0)
        throw std::invalid_argument("HI positions must contain N>=1 3D particles");
    finite_vector(q, "HI positions");
    if (!std::isfinite(parameters.radius) || !(parameters.radius > 0.0)
        || !std::isfinite(parameters.viscosity) || !(parameters.viscosity > 0.0))
        throw std::invalid_argument("HI radius and solvent viscosity must be finite and positive");
    if (parameters.model != Model::Oseen && parameters.model != Model::RPY)
        throw std::invalid_argument("Unknown HI mobility model");

    const std::size_t dimension = q.size();
    if (dimension > std::numeric_limits<std::size_t>::max() / dimension)
        throw std::length_error("HI mobility matrix dimension overflow");
    const double m0 = 1.0 / (6.0 * pi * parameters.viscosity * parameters.radius);
    if (!std::isfinite(m0) || !(m0 > 0.0))
        throw std::invalid_argument("HI isolated mobility is outside finite positive range");
    std::vector<double> matrix(dimension * dimension, 0.0);
    for (std::size_t d = 0; d < dimension; ++d)
        matrix[d * dimension + d] = m0;

    const std::size_t particles = dimension / 3;
    for (std::size_t i = 0; i < particles; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            double separation[3];
            for (std::size_t d = 0; d < 3; ++d)
                separation[d] = q[3 * i + d] - q[3 * j + d];
            const double distance = std::hypot(std::hypot(separation[0], separation[1]),
                                               separation[2]);
            if (!std::isfinite(distance))
                throw std::invalid_argument("HI pair separation is nonfinite");
            if (parameters.model == Model::Oseen
                && distance / parameters.radius < 2.0)
                throw std::domain_error(
                    "Oseen HI requires nonoverlapping far-field spheres (r>=2a); "
                    "use regularized RPY for soft overlap");

            double isotropic = 0.0;
            double parallel = 0.0;
            if (parameters.model == Model::Oseen) {
                isotropic = 0.75 * parameters.radius / distance;
                parallel = isotropic;
            } else if (distance / parameters.radius >= 2.0) {
                const double ar = parameters.radius / distance;
                isotropic = 0.75 * ar + 0.5 * ar * ar * ar;
                parallel = 0.75 * ar - 1.5 * ar * ar * ar;
            } else {
                const double ra = distance / parameters.radius;
                isotropic = 1.0 - 9.0 * ra / 32.0;
                parallel = 3.0 * ra / 32.0;
            }
            for (std::size_t a = 0; a < 3; ++a) {
                for (std::size_t b = 0; b < 3; ++b) {
                    const double dyad = distance > 0.0
                        ? (separation[a] / distance) * (separation[b] / distance) : 0.0;
                    const double entry = m0 * (isotropic * (a == b ? 1.0 : 0.0)
                                               + parallel * dyad);
                    matrix[(3 * i + a) * dimension + 3 * j + b] = entry;
                    matrix[(3 * j + b) * dimension + 3 * i + a] = entry;
                }
            }
        }
    }
    finite_vector(matrix, "HI mobility");
    return matrix;
}

std::vector<double> HI_BD::increment(const std::vector<double>& q,
                                    const std::vector<double>& force,
                                    const std::vector<double>& background,
                                    double dt, double kBT, gsl_rng* rng,
                                    const Parameters& parameters) {
    if (!std::isfinite(dt) || !(dt > 0.0) || !std::isfinite(kBT) || kBT < 0.0)
        throw std::invalid_argument("HI step requires finite dt>0 and kBT>=0");
    if (force.size() != q.size() || background.size() != q.size())
        throw std::invalid_argument("HI force/background shape must equal position shape");
    finite_vector(force, "HI force");
    finite_vector(background, "HI background velocity");
    if (kBT > 0.0 && rng == 0)
        throw std::invalid_argument("Thermal HI increment requires an RNG");
    const std::vector<double> matrix = mobility(q, parameters);
    const std::size_t dimension = q.size();
    const std::vector<double> lower = cholesky(matrix, dimension);
    std::vector<double> result(dimension, 0.0);
    std::vector<double> gaussian(dimension, 0.0);
    const double noise_scale = std::sqrt(2.0 * kBT) * std::sqrt(dt);
    if (!std::isfinite(noise_scale))
        throw std::invalid_argument("HI stochastic scale is nonfinite");
    if (kBT > 0.0)
        for (std::size_t d = 0; d < dimension; ++d)
            gaussian[d] = gsl_ran_gaussian_ziggurat(rng, 1.0);

    for (std::size_t i = 0; i < dimension; ++i) {
        double velocity = background[i];
        for (std::size_t j = 0; j < dimension; ++j)
            velocity += matrix[i * dimension + j] * force[j];
        result[i] = dt * velocity;
        if (kBT > 0.0) {
            double fluctuation = 0.0;
            for (std::size_t j = 0; j <= i; ++j)
                fluctuation += lower[i * dimension + j] * gaussian[j];
            result[i] += noise_scale * fluctuation;
        }
    }
    finite_vector(result, "HI displacement");
    return result;
}
