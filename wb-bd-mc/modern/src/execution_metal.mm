#include "bdmc/execution.hpp"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace bdmc {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
std::string description(NSError* error) {
    return error ? std::string([[error localizedDescription] UTF8String]) : "unspecified Metal error";
}
float float_value(double value) {
    require(std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max(),
            "Metal float32 input is nonfinite or outside its finite range");
    const float converted = static_cast<float>(value);
    require(value == 0 || (converted != 0 && std::abs(converted) >= std::numeric_limits<float>::min()),
            "Metal float32 input underflows or is subnormal; input rejected without clamping");
    return converted;
}
std::size_t items(const std::vector<double>& values, bool virials = false) {
    require(!values.empty() && values.size()%3 == 0, "kernel input requires nonempty N*3 storage");
    const auto n = values.size()/3;
    require(n <= std::numeric_limits<std::uint32_t>::max()/9,
            "Metal item count exceeds checked uint32 indexing range");
    require(!virials || n <= std::vector<double>().max_size()/9, "N*9 virial storage overflows");
    return n;
}

// Input/output arrays are packed scalar buffers: Metal float3 has 16-byte
// alignment and must never be used as the host N*3 storage representation.
const char* shader_source = R"metal(
#include <metal_stdlib>
using namespace metal;
struct Parameters { uint n; float a; float b; float c; uint fene; };
void report(device atomic_uint* error, uint code) {
    atomic_fetch_or_explicit(error,code,memory_order_relaxed);
}
float norm3(float3 q) {
    float largest = max(abs(q.x),max(abs(q.y),abs(q.z)));
    if (largest == 0) return 0;
    float3 scaled = q/largest;
    return largest*sqrt(dot(scaled,scaled));
}
kernel void repulsive_kernel(device const float* p [[buffer(0)]],
                             device float* forces [[buffer(1)]],
                             device float* energies [[buffer(2)]],
                             device float* virials [[buffer(3)]],
                             device atomic_uint* error [[buffer(4)]],
                             constant Parameters& parameters [[buffer(5)]],
                             uint i [[thread_position_in_grid]]) {
    if (i >= parameters.n) return;
    float box=parameters.a, epsilon=parameters.b, sigma=parameters.c;
    float3 force=0; float energy=0; float virial[9]={0,0,0,0,0,0,0,0,0};
    float3 x=float3(p[3*i],p[3*i+1],p[3*i+2]);
    for (uint j=0;j<parameters.n;++j) {
        if (j==i) continue;
        float3 dx=x-float3(p[3*j],p[3*j+1],p[3*j+2]);
        dx-=box*floor(dx/box+0.5f);
        float r=norm3(dx);
        if (r==0) { report(error,2); continue; }
        if (!isfinite(r)) { report(error,1); continue; }
        if (r>=sigma) continue;
        float ratio=r/sigma;
        float3 f=((epsilon/sigma)*(1-ratio*ratio))*(dx/r);
        force+=f;
        energy+=0.5f*(epsilon/3)*(1-ratio)*(1-ratio)*(2+ratio);
        for(uint a=0;a<3;++a) for(uint b=0;b<3;++b) virial[3*a+b]+=(0.5f*dx[a])*f[b];
    }
    energies[i]=energy;
    if(!isfinite(energy)) report(error,1);
    for(uint a=0;a<3;++a) { forces[3*i+a]=force[a]; if(!isfinite(force[a])) report(error,1); }
    for(uint a=0;a<9;++a) { virials[9*i+a]=virial[a]; if(!isfinite(virial[a])) report(error,1); }
}
kernel void springs_kernel(device const float* q [[buffer(0)]],
                           device float* forces [[buffer(1)]],
                           device float* energies [[buffer(2)]],
                           device atomic_uint* error [[buffer(4)]],
                           constant Parameters& parameters [[buffer(5)]],
                           uint i [[thread_position_in_grid]]) {
    if(i>=parameters.n) return;
    float h=parameters.a,q0=parameters.b;
    float3 connector=float3(q[3*i],q[3*i+1],q[3*i+2]);
    float norm=norm3(connector),ratio=0;
    if(!isfinite(norm)) { report(error,1); return; }
    if(parameters.fene) {
        if(!(norm<q0)) { report(error,4); return; }
        ratio=norm/q0; ratio*=ratio;
        if(!(ratio<1)) { report(error,4); return; }
    }
    float3 force=(h/(1-ratio))*connector;
    float energy;
    if(parameters.fene) {
        // Compensated log1p(-ratio), using Metal's standard log primitive.
        float remainder=1-ratio;
        float logarithm=remainder==1 ? -ratio : log(remainder)*(-ratio)/(remainder-1);
        energy=-0.5f*h*q0*q0*logarithm;
    } else energy=0.5f*h*dot(connector,connector);
    energies[i]=energy;
    if(!isfinite(energy)) report(error,1);
    for(uint a=0;a<3;++a) { forces[3*i+a]=force[a]; if(!isfinite(force[a])) report(error,1); }
}
kernel void step_kernel(device const float* old [[buffer(0)]],
                        device const float* forces [[buffer(1)]],
                        device const float* noise [[buffer(2)]],
                        device float* output [[buffer(3)]],
                        device atomic_uint* error [[buffer(4)]],
                        constant Parameters& parameters [[buffer(5)]],
                        uint k [[thread_position_in_grid]]) {
    if(k>=parameters.n) return;
    float value=old[k]+parameters.a*forces[k]+parameters.b*noise[k];
    output[k]=value;
    if(!isfinite(value)) report(error,1);
}
)metal";

struct Parameters { std::uint32_t n; float a, b, c; std::uint32_t fene; };
static_assert(sizeof(Parameters)==20, "Metal parameter layout must match shader");

struct Buffer {
    id<MTLBuffer> storage = nil;
    std::size_t capacity = 0;
    void reserve(id<MTLDevice> device, std::size_t bytes) {
        if (bytes <= capacity) return;
        if (bytes > [device maxBufferLength]) throw std::runtime_error("Metal buffer exceeds device maxBufferLength");
        id<MTLBuffer> replacement = [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
        if (!replacement) throw std::runtime_error("Metal shared buffer allocation failed");
        storage = replacement;
        capacity = bytes;
    }
};

class MetalKernels final : public ExecutionKernels {
public:
    MetalKernels() {
        @autoreleasepool {
            device_ = MTLCreateSystemDefaultDevice();
            if (!device_) throw std::runtime_error("Metal backend unavailable: no default Metal device");
            queue_ = [device_ newCommandQueue];
            if (!queue_) throw std::runtime_error("Metal command queue creation failed");
            MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
            options.fastMathEnabled = NO;
            NSError* error = nil;
            id<MTLLibrary> library = [device_ newLibraryWithSource:[NSString stringWithUTF8String:shader_source]
                                                        options:options error:&error];
            if (!library) throw std::runtime_error("Metal shader compilation failed: " + description(error));
            repulsive_ = pipeline(library,@"repulsive_kernel");
            springs_ = pipeline(library,@"springs_kernel");
            step_ = pipeline(library,@"step_kernel");
            error_.reserve(device_,sizeof(std::uint32_t));
        }
    }
    std::string name() const override { return "metal"; }
    std::string precision() const override { return "float32"; }

    KernelResult repulsive(const std::vector<double>& positions, double box,
                           double epsilon, double sigma) override {
        @autoreleasepool {
            const auto n = items(positions,true);
            require(std::isfinite(box) && box>0, "box must be finite and positive");
            require(std::isfinite(sigma) && sigma>0 && sigma<=0.5*box, "sigma must satisfy 0<sigma<=box/2");
            require(std::isfinite(epsilon) && epsilon>=0 && std::isfinite(epsilon/sigma), "epsilon/force scale is invalid");
            for(double x:positions) require(std::isfinite(x) && x>=0 && x<box,"repulsive positions must be wrapped to [0,box)");
            Parameters p{static_cast<std::uint32_t>(n),float_value(box),float_value(epsilon),float_value(sigma),0};
            require(p.c<=0.5f*p.a && std::isfinite(p.b/p.c),"Metal float32 repulsive parameters leave valid numerical range");
            upload(input_,positions);
            const auto* converted = static_cast<const float*>([input_.storage contents]);
            for(std::size_t k=0;k<positions.size();++k)
                require(converted[k]>=0 && converted[k]<p.a,"Metal float32 coordinate rounds outside wrapped [0,box)");
            reserve(output_,positions.size()); reserve(energies_,n); reserve(virials_,9*n);
            run(repulsive_,p,{input_.storage,output_.storage,energies_.storage,virials_.storage});
            KernelResult result; result.forces=download(output_,positions.size());
            result.energies=download(energies_,n); result.pair_virials=download(virials_,9*n);
            return result;
        }
    }

    KernelResult springs(const std::vector<double>& connectors,double h,double q0,bool fene) override {
        @autoreleasepool {
            const auto n=items(connectors);
            require(std::isfinite(h) && h>0,"spring constant must be finite and positive");
            require(std::isfinite(q0) && q0>0,"q0 must be finite and positive");
            // Original double-domain validation precedes float32 conversion:
            // rounding must never turn an invalid original connector into a valid one.
            for(std::size_t i=0;i<n;++i) {
                const double norm=std::hypot(connectors[3*i],connectors[3*i+1],connectors[3*i+2]);
                require(std::isfinite(norm),"connector norm is nonfinite");
                require(!fene || norm<q0,"FENE requires |q|<q0 before Metal float32 conversion");
            }
            Parameters p{static_cast<std::uint32_t>(n),float_value(h),float_value(q0),0,fene ? 1u:0u};
            upload(input_,connectors); reserve(output_,connectors.size()); reserve(energies_,n);
            run(springs_,p,{input_.storage,output_.storage,energies_.storage,nil});
            KernelResult result; result.forces=download(output_,connectors.size()); result.energies=download(energies_,n);
            return result;
        }
    }

    std::vector<double> linear_step(const std::vector<double>& old,const std::vector<double>& forces,
                                   const std::vector<double>& noise,double drift,double amplitude) override {
        @autoreleasepool {
            (void)items(old);
            require(forces.size()==old.size() && noise.size()==old.size(),"linear_step array shapes differ");
            Parameters p{static_cast<std::uint32_t>(old.size()),float_value(drift),float_value(amplitude),0,0};
            upload(input_,old); upload(forces_,forces); upload(noise_,noise); reserve(output_,old.size());
            run(step_,p,{input_.storage,forces_.storage,noise_.storage,output_.storage});
            return download(output_,old.size());
        }
    }
private:
    id<MTLComputePipelineState> pipeline(id<MTLLibrary> library,NSString* name) {
        id<MTLFunction> function=[library newFunctionWithName:name];
        if(!function) throw std::runtime_error("Metal kernel function missing");
        NSError* error=nil;
        id<MTLComputePipelineState> result=[device_ newComputePipelineStateWithFunction:function error:&error];
        if(!result) throw std::runtime_error("Metal pipeline creation failed: "+description(error));
        return result;
    }
    void reserve(Buffer& buffer,std::size_t count) {
        require(count<=std::numeric_limits<std::size_t>::max()/sizeof(float),"Metal byte count overflows");
        buffer.reserve(device_,count*sizeof(float));
    }
    void upload(Buffer& buffer,const std::vector<double>& values) {
        reserve(buffer,values.size());
        auto* converted=static_cast<float*>([buffer.storage contents]);
        for(std::size_t k=0;k<values.size();++k) converted[k]=float_value(values[k]);
    }
    std::vector<double> download(const Buffer& buffer,std::size_t count) {
        const auto* source=static_cast<const float*>([buffer.storage contents]);
        std::vector<double> result(count);
        for(std::size_t k=0;k<count;++k) {
            require(std::isfinite(source[k]),"Metal physics/update output overflows");
            result[k]=source[k];
        }
        return result;
    }
    void run(id<MTLComputePipelineState> pipeline,const Parameters& parameters,
             std::initializer_list<id<MTLBuffer>> buffers) {
        *static_cast<std::uint32_t*>([error_.storage contents])=0;
        id<MTLCommandBuffer> command=[queue_ commandBuffer];
        if(!command) throw std::runtime_error("Metal command buffer creation failed");
        id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
        if(!encoder) throw std::runtime_error("Metal compute encoder creation failed");
        [encoder setComputePipelineState:pipeline];
        NSUInteger index=0;
        for(id<MTLBuffer> buffer:buffers) [encoder setBuffer:buffer offset:0 atIndex:index++];
        [encoder setBuffer:error_.storage offset:0 atIndex:4];
        [encoder setBytes:&parameters length:sizeof(parameters) atIndex:5];
        const NSUInteger width=std::min<NSUInteger>(256,[pipeline maxTotalThreadsPerThreadgroup]);
        [encoder dispatchThreadgroups:MTLSizeMake(1+(parameters.n-1)/width,1,1)
               threadsPerThreadgroup:MTLSizeMake(width,1,1)];
        [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
        if([command status]!=MTLCommandBufferStatusCompleted)
            throw std::runtime_error("Metal command failed: "+description([command error]));
        const auto error=*static_cast<const std::uint32_t*>([error_.storage contents]);
        if(error&2) throw std::invalid_argument("coincident periodic particles in Metal float32 have undefined repulsive gradient");
        if(error&4) throw std::invalid_argument("FENE requires |q|<q0 in Metal float32 arithmetic; no clamping");
        if(error) throw std::invalid_argument("Metal float32 physics/update output overflows");
    }
    id<MTLDevice> device_ = nil;
    id<MTLCommandQueue> queue_ = nil;
    id<MTLComputePipelineState> repulsive_ = nil, springs_ = nil, step_ = nil;
    Buffer input_, forces_, noise_, output_, energies_, virials_, error_;
};

} // namespace

BackendCapability metal_capability() {
    @autoreleasepool {
        id<MTLDevice> device=MTLCreateSystemDefaultDevice();
        if(!device) return {"metal","float32","no default Metal device",true,false};
        return {"metal","float32",std::string([[device name] UTF8String])+"; float32 only; host-owned RNG",true,true};
    }
}
std::unique_ptr<ExecutionKernels> make_metal_kernels() { return std::make_unique<MetalKernels>(); }

} // namespace bdmc
