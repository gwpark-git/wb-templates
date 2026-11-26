#include "bdmc/execution.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace bdmc {
namespace {

void check_cuda(cudaError_t status, const char* operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string("CUDA ") + operation + ": " + cudaGetErrorString(status));
}
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
void finite_array(const std::vector<double>& values) {
    for (double x : values) require(std::isfinite(x), "kernel input must be finite");
}
std::size_t items(const std::vector<double>& values, bool virials = false) {
    require(!values.empty() && values.size() % 3 == 0, "kernel input requires nonempty N*3 storage");
    const auto n = values.size()/3;
    require(!virials || n <= std::vector<double>().max_size()/9, "N*9 virial storage overflows");
    finite_array(values);
    return n;
}

class Stream {
public:
    Stream() { check_cuda(cudaStreamCreateWithFlags(&value, cudaStreamNonBlocking), "create stream"); }
    ~Stream() { if (value) (void)cudaStreamDestroy(value); }
    Stream(const Stream&) = delete;
    Stream& operator=(const Stream&) = delete;
    cudaStream_t value = nullptr;
};
template<class T> class Buffer {
public:
    ~Buffer() { if (data_) (void)cudaFree(data_); }
    Buffer() = default;
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    void reserve(std::size_t count) {
        if (count <= capacity_) return;
        require(count <= std::numeric_limits<std::size_t>::max()/sizeof(T), "device buffer size overflows");
        T* replacement = nullptr;
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&replacement), count*sizeof(T)), "allocate device buffer");
        if (data_) {
            const auto status = cudaFree(data_);
            if (status != cudaSuccess) { (void)cudaFree(replacement); check_cuda(status, "release device buffer"); }
        }
        data_ = replacement;
        capacity_ = count;
    }
    T* data() { return data_; }
private:
    T* data_ = nullptr;
    std::size_t capacity_ = 0;
};

enum DeviceError : unsigned int { nonfinite = 1, coincident = 2, fene_domain = 4 };
__device__ void report(unsigned int* error, unsigned int value) { atomicOr(error, value); }

__global__ void repulsive_kernel(const double* positions, std::size_t n,
                                double box, double epsilon, double sigma,
                                double* forces, double* energies, double* virials,
                                unsigned int* error) {
    const auto stride = static_cast<std::size_t>(blockDim.x)*gridDim.x;
    for (std::size_t i = static_cast<std::size_t>(blockIdx.x)*blockDim.x + threadIdx.x; i < n; i += stride) {
        double force[3] = {0,0,0}, virial[9] = {}, energy = 0;
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            double dx[3];
            for (int a = 0; a < 3; ++a) {
                dx[a] = positions[3*i+a] - positions[3*j+a];
                dx[a] -= box*floor(dx[a]/box + 0.5);
            }
            const double r = hypot(hypot(dx[0],dx[1]),dx[2]);
            if (r == 0) { report(error, coincident); continue; }
            if (!isfinite(r)) { report(error, nonfinite); continue; }
            if (r >= sigma) continue;
            const double ratio = r/sigma;
            const double magnitude = (epsilon/sigma)*(1-ratio*ratio);
            const double u = (epsilon/3)*(1-ratio)*(1-ratio)*(2+ratio);
            double f[3];
            energy += 0.5*u;
            for (int a = 0; a < 3; ++a) { f[a] = magnitude*(dx[a]/r); force[a] += f[a]; }
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) virial[3*a+b] += (0.5*dx[a])*f[b];
        }
        energies[i] = energy;
        if (!isfinite(energy)) report(error, nonfinite);
        for (int a = 0; a < 3; ++a) {
            forces[3*i+a] = force[a];
            if (!isfinite(force[a])) report(error, nonfinite);
        }
        for (int a = 0; a < 9; ++a) {
            virials[9*i+a] = virial[a];
            if (!isfinite(virial[a])) report(error, nonfinite);
        }
    }
}

__global__ void springs_kernel(const double* q, std::size_t n, double h, double q0, bool fene,
                              double* forces, double* energies, unsigned int* error) {
    const auto stride = static_cast<std::size_t>(blockDim.x)*gridDim.x;
    for (std::size_t i = static_cast<std::size_t>(blockIdx.x)*blockDim.x + threadIdx.x; i < n; i += stride) {
        const double norm = hypot(hypot(q[3*i],q[3*i+1]),q[3*i+2]);
        double ratio = 0;
        if (!isfinite(norm)) { report(error, nonfinite); continue; }
        if (fene) {
            if (!(norm < q0)) { report(error, fene_domain); continue; }
            ratio = norm/q0; ratio *= ratio;
            if (!(ratio < 1)) { report(error, fene_domain); continue; }
        }
        const double factor = h/(1-ratio);
        const double energy = fene ? -0.5*h*q0*q0*log1p(-ratio)
                                  : 0.5*h*(q[3*i]*q[3*i]+q[3*i+1]*q[3*i+1]+q[3*i+2]*q[3*i+2]);
        energies[i] = energy;
        if (!isfinite(energy)) report(error, nonfinite);
        for (int a = 0; a < 3; ++a) {
            const double f = factor*q[3*i+a];
            forces[3*i+a] = f;
            if (!isfinite(f)) report(error, nonfinite);
        }
    }
}

__global__ void step_kernel(const double* old, const double* forces, const double* noise,
                           std::size_t count, double drift, double amplitude,
                           double* output, unsigned int* error) {
    const auto stride = static_cast<std::size_t>(blockDim.x)*gridDim.x;
    for (std::size_t k = static_cast<std::size_t>(blockIdx.x)*blockDim.x + threadIdx.x; k < count; k += stride) {
        const double x = old[k] + drift*forces[k] + amplitude*noise[k];
        output[k] = x;
        if (!isfinite(x)) report(error, nonfinite);
    }
}

class CudaKernels final : public ExecutionKernels {
public:
    CudaKernels() {
        check_cuda(cudaGetDevice(&device_), "get current device");
        error_.reserve(1);
    }
    ~CudaKernels() override {
        // Restore this instance's device before member buffers/stream release.
        (void)cudaSetDevice(device_);
    }
    std::string name() const override { return "cuda"; }
    std::string precision() const override { return "float64"; }

    KernelResult repulsive(const std::vector<double>& positions, double box,
                           double epsilon, double sigma) override {
        const auto n = items(positions, true);
        require(std::isfinite(box) && box > 0, "box must be finite and positive");
        require(std::isfinite(sigma) && sigma > 0 && sigma <= 0.5*box, "sigma must satisfy 0<sigma<=box/2");
        require(std::isfinite(epsilon) && epsilon >= 0 && std::isfinite(epsilon/sigma), "epsilon/force scale is invalid");
        for (double x : positions) require(x >= 0 && x < box, "repulsive positions must be wrapped to [0,box)");
        prepare(); upload(input_, positions);
        output_.reserve(positions.size()); energies_.reserve(n); virials_.reserve(9*n);
        repulsive_kernel<<<blocks(n),256,0,stream_.value>>>(input_.data(),n,box,epsilon,sigma,
            output_.data(),energies_.data(),virials_.data(),error_.data());
        check_cuda(cudaGetLastError(), "launch repulsive kernel");
        finish();
        KernelResult result; result.forces.resize(positions.size()); result.energies.resize(n); result.pair_virials.resize(9*n);
        download(result.forces,output_); download(result.energies,energies_); download(result.pair_virials,virials_);
        return result;
    }

    KernelResult springs(const std::vector<double>& connectors, double h, double q0, bool fene) override {
        const auto n = items(connectors);
        require(std::isfinite(h) && h > 0, "spring constant must be finite and positive");
        require(std::isfinite(q0) && q0 > 0, "q0 must be finite and positive");
        for (std::size_t i = 0; i < n; ++i) {
            const double norm = std::hypot(connectors[3*i],connectors[3*i+1],connectors[3*i+2]);
            require(std::isfinite(norm), "connector norm overflows");
            require(!fene || norm < q0, "FENE requires |q|<q0 before CUDA execution");
        }
        prepare(); upload(input_,connectors); output_.reserve(connectors.size()); energies_.reserve(n);
        springs_kernel<<<blocks(n),256,0,stream_.value>>>(input_.data(),n,h,q0,fene,output_.data(),energies_.data(),error_.data());
        check_cuda(cudaGetLastError(), "launch spring kernel");
        finish();
        KernelResult result; result.forces.resize(connectors.size()); result.energies.resize(n);
        download(result.forces,output_); download(result.energies,energies_);
        return result;
    }

    std::vector<double> linear_step(const std::vector<double>& old, const std::vector<double>& forces,
                                   const std::vector<double>& noise, double drift, double amplitude) override {
        (void)items(old);
        require(forces.size() == old.size() && noise.size() == old.size(), "linear_step array shapes differ");
        finite_array(forces); finite_array(noise);
        require(std::isfinite(drift) && std::isfinite(amplitude), "linear_step scales must be finite");
        prepare(); upload(input_,old); upload(forces_,forces); upload(noise_,noise); output_.reserve(old.size());
        step_kernel<<<blocks(old.size()),256,0,stream_.value>>>(input_.data(),forces_.data(),noise_.data(),old.size(),drift,amplitude,output_.data(),error_.data());
        check_cuda(cudaGetLastError(), "launch linear step kernel");
        finish();
        std::vector<double> result(old.size()); download(result,output_);
        return result;
    }
private:
    static unsigned int blocks(std::size_t count) {
        return static_cast<unsigned int>(std::min<std::size_t>(1+(count-1)/256,65535));
    }
    void prepare() {
        check_cuda(cudaSetDevice(device_), "select instance device");
        check_cuda(cudaMemsetAsync(error_.data(),0,sizeof(unsigned int),stream_.value), "reset device error");
    }
    void upload(Buffer<double>& buffer, const std::vector<double>& values) {
        buffer.reserve(values.size());
        check_cuda(cudaMemcpy(buffer.data(),values.data(),values.size()*sizeof(double),cudaMemcpyHostToDevice), "upload input");
    }
    void download(std::vector<double>& values, Buffer<double>& buffer) {
        check_cuda(cudaMemcpy(values.data(),buffer.data(),values.size()*sizeof(double),cudaMemcpyDeviceToHost), "download output");
    }
    void finish() {
        unsigned int error = 0;
        check_cuda(cudaStreamSynchronize(stream_.value), "synchronize kernels");
        check_cuda(cudaMemcpy(&error,error_.data(),sizeof(error),cudaMemcpyDeviceToHost), "download device error");
        if (error & coincident) throw std::invalid_argument("coincident periodic particles have undefined repulsive gradient");
        if (error & fene_domain) throw std::invalid_argument("FENE requires |q|<q0 in CUDA float64 arithmetic");
        if (error) throw std::invalid_argument("CUDA physics/update output overflows");
    }
    int device_ = 0;
    Stream stream_;
    Buffer<double> input_, forces_, noise_, output_, energies_, virials_;
    Buffer<unsigned int> error_;
};

} // namespace

BackendCapability cuda_capability() {
    BackendCapability result{"cuda","float64","",true,false};
    int count = 0;
    const auto count_status = cudaGetDeviceCount(&count);
    if (count_status != cudaSuccess) { result.reason = cudaGetErrorString(count_status); (void)cudaGetLastError(); return result; }
    if (count == 0) { result.reason = "no CUDA devices"; return result; }
    int device = 0;
    const auto current_status = cudaGetDevice(&device);
    if (current_status != cudaSuccess) { result.reason = cudaGetErrorString(current_status); return result; }
    cudaDeviceProp property{};
    const auto property_status = cudaGetDeviceProperties(&property,device);
    if (property_status != cudaSuccess) { result.reason = cudaGetErrorString(property_status); return result; }
    if (property.major < 1 || (property.major == 1 && property.minor < 3)) {
        result.reason = "current CUDA device does not support float64"; return result;
    }
    result.available = true;
    result.reason = std::string(property.name) + "; CUDA float64 kernels; host-owned RNG";
    return result;
}
std::unique_ptr<ExecutionKernels> make_cuda_kernels() {
    const auto capability = cuda_capability();
    if (!capability.available) throw std::runtime_error("CUDA backend unavailable: " + capability.reason);
    return std::make_unique<CudaKernels>();
}

} // namespace bdmc
