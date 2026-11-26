#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace bdmc {

struct BackendCapability {
    std::string name, precision, reason;
    bool compiled = false, available = false;
};

// Bulk host/device boundary. No per-particle virtual calls and no device RNG.
// Forces and per-item energies are returned in double host storage; Metal's
// arithmetic is explicitly float32. Pair virials have shape (N,3,3), counted
// half per directed neighbor. Spring virials are reconstructed from q and F.
struct KernelResult {
    std::vector<double> forces, energies, pair_virials;
};

// A backend replaces bulk arithmetic, not the physical controller or random
// process. The caller supplies the same ordered noise values and commits state.
// A compiled/available backend still needs its own numerical acceptance checks;
// host float64 output does not upgrade Metal's float32 arithmetic.
class ExecutionKernels {
public:
    virtual ~ExecutionKernels() = default;
    virtual std::string name() const = 0;
    virtual std::string precision() const = 0;
    virtual KernelResult repulsive(const std::vector<double>& positions,
                                  double box, double epsilon, double sigma) = 0;
    virtual KernelResult springs(const std::vector<double>& connectors,
                                double h, double q0, bool fene) = 0;
    // old + drift_scale*force + noise_scale*z; caller owns the RNG transaction.
    virtual std::vector<double> linear_step(const std::vector<double>& old,
                                           const std::vector<double>& forces,
                                           const std::vector<double>& noise,
                                           double drift_scale, double noise_scale) = 0;
};

std::vector<BackendCapability> execution_backends();
std::unique_ptr<ExecutionKernels> make_execution_kernels(const std::string& backend,
                                                       int threads = 1);
void validate_execution_selection(const std::string& backend, int threads);

// Optional implementations, linked only when their CMake option is enabled.
BackendCapability cuda_capability();
std::unique_ptr<ExecutionKernels> make_cuda_kernels();
BackendCapability metal_capability();
std::unique_ptr<ExecutionKernels> make_metal_kernels();
BackendCapability accelerate_capability();
std::unique_ptr<ExecutionKernels> make_accelerate_kernels();

} // namespace bdmc
