#include "bdmc/execution.hpp"

#include <Accelerate/Accelerate.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace bdmc {
namespace {

// Accelerate deliberately uses the serial force oracle and CBLAS for the bulk
// Euler update. This is an Apple CPU backend, not a GPU implementation.
class AccelerateKernels final : public ExecutionKernels {
public:
    AccelerateKernels() : physics_(make_execution_kernels("serial")) {}
    std::string name() const override { return "accelerate"; }
    std::string precision() const override { return "float64"; }

    KernelResult repulsive(const std::vector<double>& positions, double box,
                           double epsilon, double sigma) override {
        return physics_->repulsive(positions, box, epsilon, sigma);
    }
    KernelResult springs(const std::vector<double>& connectors, double h,
                         double q0, bool fene) override {
        return physics_->springs(connectors, h, q0, fene);
    }
    std::vector<double> linear_step(const std::vector<double>& old,
                                   const std::vector<double>& forces,
                                   const std::vector<double>& noise,
                                   double drift_scale, double noise_scale) override {
        if (old.empty() || old.size() % 3 || forces.size() != old.size() ||
            noise.size() != old.size())
            throw std::invalid_argument("linear_step requires matching nonempty N*3 arrays");
        if (!std::isfinite(drift_scale) || !std::isfinite(noise_scale))
            throw std::invalid_argument("linear_step scales must be finite");
        for (std::size_t i = 0; i < old.size(); ++i)
            if (!std::isfinite(old[i]) || !std::isfinite(forces[i]) || !std::isfinite(noise[i]))
                throw std::invalid_argument("linear_step inputs must be finite");
        std::vector<double> result(old.size());
        // Use the default LP64 CBLAS interface, with bounded chunks even when
        // the host array exceeds the signed BLAS integer range.
        const auto max_chunk = static_cast<std::size_t>(std::numeric_limits<int>::max());
        for (std::size_t offset = 0; offset < old.size();) {
            const auto count = static_cast<int>(std::min(old.size() - offset, max_chunk));
            cblas_dcopy(count, old.data() + offset, 1, result.data() + offset, 1);
            cblas_daxpy(count, drift_scale, forces.data() + offset, 1, result.data() + offset, 1);
            cblas_daxpy(count, noise_scale, noise.data() + offset, 1, result.data() + offset, 1);
            offset += static_cast<std::size_t>(count);
        }
        for (double x : result)
            if (!std::isfinite(x)) throw std::invalid_argument("linear_step output overflows");
        return result;
    }
private:
    std::unique_ptr<ExecutionKernels> physics_;
};

} // namespace

BackendCapability accelerate_capability() {
    return {"accelerate", "float64", "Apple Accelerate CBLAS CPU updates; serial force kernels", true, true};
}
std::unique_ptr<ExecutionKernels> make_accelerate_kernels() {
    return std::make_unique<AccelerateKernels>();
}

} // namespace bdmc
