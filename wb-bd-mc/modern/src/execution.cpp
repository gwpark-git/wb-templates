#include "bdmc/execution.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace bdmc {
namespace {
void require(bool value, const char* message) {
    if (!value) throw std::invalid_argument(message);
}
std::size_t coordinates(const std::vector<double>& values) {
    require(!values.empty() && values.size()%3==0, "backend expects nonempty Nx3 coordinates");
    require(values.size()/3 <= std::vector<double>().max_size()/9, "backend tensor storage overflows");
    for(double x:values) require(std::isfinite(x), "backend coordinates must be finite");
    return values.size()/3;
}
void finite_result(const KernelResult& r) {
    for(const auto* values:{&r.forces,&r.energies,&r.pair_virials})
        for(double x:*values) require(std::isfinite(x), "backend force/energy/virial overflows");
}
class CpuKernels final:public ExecutionKernels {
public:
    CpuKernels(std::string name,int threads):name_(std::move(name)),threads_(threads){}
    std::string name() const override{return name_;}
    std::string precision() const override{return "float64";}
    KernelResult repulsive(const std::vector<double>& x,double box,double epsilon,double sigma) override {
        const auto n=coordinates(x);
        require(std::isfinite(box)&&box>0&&std::isfinite(sigma)&&sigma>0&&sigma<=box/2,
                "backend requires positive box and 0<sigma<=box/2");
        require(std::isfinite(epsilon)&&epsilon>=0&&std::isfinite(epsilon/sigma),"backend repulsion parameters invalid");
        for(double v:x) require(v>=0&&v<box,"backend pair positions must be wrapped into [0,L)");
        KernelResult out{std::vector<double>(n*3),std::vector<double>(n),std::vector<double>(n*9)};
        std::atomic<bool> coincidence{false};
        // One output row per worker: no atomics on forces and no concurrent
        // reduction. Summation order within a row is independent of threads.
#ifdef BDMC_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(threads_) if(threads_>1)
#endif
        for(std::size_t i=0;i<n;++i) {
            for(std::size_t j=0;j<n;++j) if(i!=j) {
                double dx[3];
                for(std::size_t c=0;c<3;++c) {dx[c]=x[3*i+c]-x[3*j+c];dx[c]-=box*std::floor(dx[c]/box+.5);}
                const double r=std::hypot(dx[0],dx[1],dx[2]);
                if(r==0) {coincidence.store(true,std::memory_order_relaxed);continue;}
                if(r>=sigma) continue;
                const double u=r/sigma, radial=(epsilon/sigma)*(1-u*u), energy=(epsilon/3)*(1-u)*(1-u)*(2+u);
                out.energies[i]+=.5*energy;
                for(std::size_t a=0;a<3;++a) {
                    const double f=radial*(dx[a]/r);out.forces[3*i+a]+=f;
                    for(std::size_t b=0;b<3;++b) out.pair_virials[9*i+3*a+b]+=.5*dx[b]*f;
                }
            }
        }
        require(!coincidence.load(),"coincident periodic particles have undefined repulsive gradient");
        finite_result(out);return out;
    }
    KernelResult springs(const std::vector<double>& q,double h,double q0,bool fene) override {
        const auto n=coordinates(q);
        require(std::isfinite(h)&&h>0&&std::isfinite(q0)&&q0>0,"backend spring parameters invalid");
        for(std::size_t i=0;i<n;++i) {
            const double r=std::hypot(q[3*i],q[3*i+1],q[3*i+2]);
            require(std::isfinite(r)&&(!fene||r<q0),"FENE requires |q|<q0; backend has not clamped the state");
        }
        KernelResult out{std::vector<double>(n*3),std::vector<double>(n),{}};
#ifdef BDMC_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(threads_) if(threads_>1)
#endif
        for(std::size_t i=0;i<n;++i) {
            const double r=std::hypot(q[3*i],q[3*i+1],q[3*i+2]);
            const double ratio=fene?(r/q0)*(r/q0):0., factor=h/(1-ratio);
            out.energies[i]=fene?-.5*h*q0*q0*std::log1p(-ratio):.5*h*r*r;
            for(std::size_t c=0;c<3;++c) out.forces[3*i+c]=factor*q[3*i+c];
        }
        finite_result(out);return out;
    }
    std::vector<double> linear_step(const std::vector<double>& old,const std::vector<double>& forces,
                                    const std::vector<double>& noise,double drift,double amplitude) override {
        coordinates(old);
        require(forces.size()==old.size()&&noise.size()==old.size(),"backend Euler arrays must have equal lengths");
        require(std::isfinite(drift)&&std::isfinite(amplitude),"backend Euler scales must be finite");
        for(const auto* values:{&forces,&noise}) for(double x:*values) require(std::isfinite(x),"backend Euler inputs must be finite");
        std::vector<double> result(old.size());
#ifdef BDMC_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(threads_) if(threads_>1)
#endif
        for(std::size_t k=0;k<old.size();++k) result[k]=old[k]+drift*forces[k]+amplitude*noise[k];
        for(double x:result) require(std::isfinite(x),"backend Euler result overflows");
        return result;
    }
private:std::string name_;int threads_;
};
} // namespace

std::vector<BackendCapability> execution_backends() {
    std::vector<BackendCapability> result{{"serial","float64","CPU reference",true,true}};
#ifdef BDMC_HAVE_OPENMP
    result.push_back({"openmp","float64","parallel bulk CPU kernels",true,true});
#else
    result.push_back({"openmp","float64","configure BDMC_ENABLE_OPENMP=ON",false,false});
#endif
#ifdef BDMC_HAVE_CUDA
    result.push_back(cuda_capability());
#else
    result.push_back({"cuda","float64","configure BDMC_ENABLE_CUDA=ON with CUDA toolkit",false,false});
#endif
#ifdef BDMC_HAVE_METAL
    result.push_back(metal_capability());
#else
    result.push_back({"metal","float32","requires macOS and BDMC_ENABLE_METAL=ON",false,false});
#endif
#ifdef BDMC_HAVE_ACCELERATE
    result.push_back(accelerate_capability());
#else
    result.push_back({"accelerate","float64","requires macOS and BDMC_ENABLE_ACCELERATE=ON",false,false});
#endif
    return result;
}
void validate_execution_selection(const std::string& backend,int threads) {
    // An unavailable requested backend is an error, not an automatic fallback:
    // silently switching precision or execution would invalidate provenance and
    // make a claimed GPU/parallel comparison impossible to interpret.
    require(threads>0&&threads<=256,"backend_threads must be in 1..256");
    require(backend=="openmp"||threads==1,"backend_threads>1 requires the explicit openmp backend");
    if(backend=="serial") return;
    for(const auto& c:execution_backends()) if(c.name==backend) {
        if(!c.available) throw std::invalid_argument("requested execution backend "+backend+" unavailable: "+c.reason);
        return;
    }
    throw std::invalid_argument("execution_backend must be serial, openmp, cuda, metal or accelerate");
}
std::unique_ptr<ExecutionKernels> make_execution_kernels(const std::string& backend,int threads) {
    validate_execution_selection(backend,threads);
    if(backend=="serial"||backend=="openmp") return std::make_unique<CpuKernels>(backend,threads);
#ifdef BDMC_HAVE_CUDA
    if(backend=="cuda") return make_cuda_kernels();
#endif
#ifdef BDMC_HAVE_METAL
    if(backend=="metal") return make_metal_kernels();
#endif
#ifdef BDMC_HAVE_ACCELERATE
    if(backend=="accelerate") return make_accelerate_kernels();
#endif
    throw std::invalid_argument("requested execution backend was not built");
}
} // namespace bdmc
