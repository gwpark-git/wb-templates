#include "bdmc/reference.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace bdmc {
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
void finite(double x, const char* message) { require(std::isfinite(x), message); }
template<class Spec> Spec validated(Spec spec) { spec.validate(); return spec; }
void common(std::size_t n, double dt, double mobility, double diffusion, NoisePolicy noise) {
    require(n > 0 && n <= std::vector<double>().max_size()/3, "count must be positive and fit N*3 storage");
    require(std::isfinite(dt) && dt > 0, "dt must be finite and positive");
    require(std::isfinite(mobility) && mobility > 0, "bead mobility must be finite and positive");
    require(std::isfinite(diffusion) && diffusion >= 0, "bead diffusion must be finite and nonnegative");
    finite(diffusion/mobility, "thermal energy D/mobility overflows");
    finite(4*diffusion*dt, "noise increment variance overflows");
    (void)noise_name(noise);
}
void horizon(std::uint64_t current, std::uint64_t steps, double dt) {
    require(steps <= std::numeric_limits<std::uint64_t>::max()-current, "step counter overflow");
    finite(static_cast<double>(current+steps)*dt, "physical time overflows");
}
void check_values(const std::vector<double>& values, std::size_t n) {
    require(values.size() == n*3, "initial state must have N*3 values");
    for (auto x : values) finite(x, "coordinates must be finite");
}
void pair_parameters(double r, double epsilon, double sigma) {
    require(std::isfinite(r) && r >= 0, "distance must be finite and nonnegative");
    require(std::isfinite(epsilon) && epsilon >= 0, "epsilon must be finite and nonnegative");
    require(std::isfinite(sigma) && sigma > 0, "sigma must be finite and positive");
}
double spring_ratio(const std::array<double, 3>& q, const DumbbellSpec& spec) {
    double norm = 0;
    for (double x : q) { finite(x, "connector must be finite"); norm = std::hypot(norm, x); }
    finite(norm, "connector norm overflows");
    if (spec.spring_law == SpringLaw::FENE) {
        require(norm < spec.q0, "FENE requires |q|<q0; explicit Euler step rejected without state/RNG change");
        const double ratio = norm/spec.q0;
        return ratio*ratio;
    }
    return 0;
}
std::vector<double> regular_grid(std::size_t n, double length) {
    auto side = static_cast<std::size_t>(std::ceil(std::cbrt(static_cast<double>(n))));
    // ceil(cbrt) may round down at a cube boundary. Long double avoids integer products overflowing.
    while (static_cast<long double>(side)*side*side < n) ++side;
    std::vector<double> result(n*3);
    for (std::size_t i=0; i<n; ++i) {
        const std::array<std::size_t,3> cell{{i%side, (i/side)%side, (i/side)/side}};
        for (std::size_t c=0; c<3; ++c) result[i*3+c] = (static_cast<double>(cell[c])+0.5)*length/side;
    }
    return result;
}
} // namespace

void RepulsiveSpec::validate() const {
    common(n_particles,dt,mobility,diffusion,noise_policy);
    require(std::isfinite(box_length) && box_length>0, "box_length must be finite and positive");
    pair_parameters(0,epsilon,sigma);
    finite(epsilon/sigma, "repulsive force scale overflows");
    SpatialSearchSpec{{box_length,box_length,box_length},sigma,force_search}.validate();
    validate_execution_selection(execution_backend,backend_threads);
    require(execution_backend=="serial"||force_search==SpatialBackend::BruteForce,
            "nonserial repulsive backends require explicit force_search='brute_force' (all-pairs kernels)");
}
double repulsive_pair_energy(double r, double epsilon, double sigma) {
    pair_parameters(r,epsilon,sigma);
    if (r>=sigma) return 0;
    const double x = r/sigma;
    return (epsilon/3)*(1-x)*(1-x)*(2+x);
}
double repulsive_radial_force(double r, double epsilon, double sigma) {
    pair_parameters(r,epsilon,sigma);
    if (r>=sigma) return 0;
    const double x = r/sigma;
    return (epsilon/sigma)*(1-x*x);
}
RepulsiveSimulation::RepulsiveSimulation(RepulsiveSpec spec, std::vector<double> positions)
    : spec_(validated(std::move(spec))), explicit_initial_state_(!positions.empty()),
      kernels_(spec_.execution_backend=="serial"?nullptr:make_execution_kernels(spec_.execution_backend,spec_.backend_threads)),
      box_({spec_.box_length,spec_.box_length,spec_.box_length}),
      search_({box_.lengths(),spec_.sigma,spec_.force_search}), noise_(spec_.seed,spec_.noise_policy) {
    if (positions.empty()) positions=regular_grid(spec_.n_particles,spec_.box_length);
    check_values(positions,spec_.n_particles);
    auto unwrapped=positions;
    for (std::size_t k=0;k<positions.size();++k) positions[k]=box_.wrap(positions[k],k%3);
    state_=evaluate(std::move(positions),std::move(unwrapped),0);
}
RepulsiveFrame RepulsiveSimulation::evaluate(std::vector<double> positions, std::vector<double> unwrapped,
                                             std::uint64_t step) const {
    RepulsiveFrame frame;
    frame.positions=std::move(positions); frame.unwrapped_positions=std::move(unwrapped);
    frame.forces.assign(frame.positions.size(),0); frame.step=step; frame.time=static_cast<double>(step)*spec_.dt;
    if(kernels_) {
        const auto result=kernels_->repulsive(frame.positions,spec_.box_length,spec_.epsilon,spec_.sigma);
        require(result.forces.size()==frame.positions.size()&&result.energies.size()==spec_.n_particles&&
                result.pair_virials.size()==spec_.n_particles*9,"backend repulsive output shape invalid");
        frame.forces=result.forces;
        for(std::size_t i=0;i<spec_.n_particles;++i) {
            frame.potential_energy+=result.energies[i];
            for(std::size_t a=0;a<9;++a) frame.pair_virial[a]+=result.pair_virials[9*i+a];
        }
    } else {
      const auto neighbors=search_(frame.positions);
      for (const auto& pair:neighbors.pairs) {
        require(pair.distance>0, "coincident distinct periodic particles have an undefined repulsive gradient");
        const double force=repulsive_radial_force(pair.distance,spec_.epsilon,spec_.sigma);
        frame.potential_energy+=repulsive_pair_energy(pair.distance,spec_.epsilon,spec_.sigma);
        for (std::size_t a=0;a<3;++a) {
            const double f=force*(pair.displacement[a]/pair.distance);
            frame.forces[pair.first_index*3+a]-=f;
            frame.forces[pair.second_index*3+a]+=f;
            for (std::size_t b=0;b<3;++b) frame.pair_virial[a*3+b]+=pair.displacement[b]*f;
        }
      }
    }
    finite(frame.potential_energy,"potential energy overflow");
    for (auto f:frame.forces) finite(f,"repulsive force overflow");
    for (auto w:frame.pair_virial) finite(w,"pair virial overflow");
    return frame;
}
void RepulsiveSimulation::advance_unlocked(std::uint64_t steps) {
    horizon(state_.step,steps,spec_.dt);
    const double amplitude=std::sqrt(2*spec_.diffusion*spec_.dt);
    for (std::uint64_t i=0;i<steps;++i) {
        // Speculative state and RNG advance together. A rejected/nonfinite
        // candidate must not consume a different random path on the next call.
        // This per-step guarantee differs from rollback of an entire batch.
        auto next_noise=noise_;
        auto positions=state_.positions, unwrapped=state_.unwrapped_positions;
        std::vector<double> displacement;
        if(kernels_) {
            std::vector<double> z(positions.size(),0),zero(positions.size(),0);
            if(amplitude>0) for(auto& x:z) x=next_noise.standardized();
            displacement=kernels_->linear_step(zero,state_.forces,z,spec_.mobility*spec_.dt,amplitude);
            require(displacement.size()==positions.size(),"backend displacement output shape invalid");
        }
        for (std::size_t k=0;k<positions.size();++k) {
            const double dx=kernels_?displacement[k]:spec_.mobility*state_.forces[k]*spec_.dt+
                (amplitude>0 ? amplitude*next_noise.standardized():0);
            finite(dx,"repulsive displacement overflow");
            unwrapped[k]+=dx; finite(unwrapped[k],"unwrapped position overflow");
            positions[k]=box_.wrap(positions[k]+dx,k%3);
        }
        auto candidate=evaluate(std::move(positions),std::move(unwrapped),state_.step+1);
        state_=std::move(candidate); noise_=std::move(next_noise);
    }
}
void RepulsiveSimulation::advance(std::uint64_t steps) { std::lock_guard<std::mutex> lock(mutex_); advance_unlocked(steps); }
RepulsiveFrame RepulsiveSimulation::observe() const { std::lock_guard<std::mutex> lock(mutex_); return state_; }
std::vector<RepulsiveFrame> RepulsiveSimulation::sample(std::uint64_t steps,std::uint64_t stride,bool initial) {
    std::lock_guard<std::mutex> lock(mutex_); require(stride>0,"stride must be positive"); horizon(state_.step,steps,spec_.dt);
    std::vector<RepulsiveFrame> result; if(initial) result.push_back(state_);
    for(std::uint64_t done=0;done<steps;) { const auto chunk=std::min(stride,steps-done); advance_unlocked(chunk); done+=chunk; result.push_back(state_); }
    return result;
}

std::string spring_law_name(SpringLaw law) {
    if(law==SpringLaw::Hookean) return "hookean";
    if(law==SpringLaw::FENE) return "fene";
    throw std::invalid_argument("unknown spring law");
}
SpringLaw parse_spring_law(const std::string& name) {
    if(name=="hookean") return SpringLaw::Hookean;
    if(name=="fene") return SpringLaw::FENE;
    throw std::invalid_argument("spring_law must be 'hookean' or 'fene'");
}
std::string flow_protocol_name(FlowProtocol protocol) {
    switch(protocol) {
    case FlowProtocol::NoFlow:return "no_flow";
    case FlowProtocol::SimpleShear:return "simple_shear";
    case FlowProtocol::StepShear:return "step_shear";
    case FlowProtocol::OscillatoryShear:return "oscillatory_shear";
    }
    throw std::invalid_argument("unknown flow protocol");
}
FlowProtocol parse_flow_protocol(const std::string& name) {
    if(name=="no_flow") return FlowProtocol::NoFlow;
    if(name=="simple_shear") return FlowProtocol::SimpleShear;
    if(name=="step_shear") return FlowProtocol::StepShear;
    if(name=="oscillatory_shear") return FlowProtocol::OscillatoryShear;
    throw std::invalid_argument("unsupported dumbbell flow protocol");
}
void DumbbellSpec::validate() const {
    common(n_dumbbells,dt,mobility,diffusion,noise_policy);
    require(std::isfinite(spring_constant)&&spring_constant>0,"spring_constant must be finite and positive");
    require(std::isfinite(q0)&&q0>0,"q0 must be finite and positive");
    require(mobility*spring_constant*dt<1,"requires mobility*spring_constant*dt<1 for the Hookean Euler limit");
    (void)spring_law_name(spring_law); (void)flow_protocol_name(protocol);
    finite(shear_rate,"shear_rate must be finite"); finite(step_strain,"step_strain must be finite");
    finite(strain_amplitude,"strain_amplitude must be finite");
    require(std::isfinite(angular_frequency)&&angular_frequency>=0,"angular_frequency must be finite and nonnegative");
    finite(strain_amplitude*angular_frequency,"oscillatory rate amplitude overflows");
    validate_execution_selection(execution_backend,backend_threads);
}
namespace {
double spring_energy_validated(const std::array<double,3>& q,const DumbbellSpec& spec) {
    const double ratio=spring_ratio(q,spec);
    double energy;
    if(spec.spring_law==SpringLaw::FENE) energy=-0.5*spec.spring_constant*spec.q0*spec.q0*std::log1p(-ratio);
    else energy=0.5*spec.spring_constant*(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
    finite(energy,"spring energy overflow"); return energy;
}
std::array<double,3> spring_force_validated(const std::array<double,3>& q,const DumbbellSpec& spec) {
    const double ratio=spring_ratio(q,spec);
    const double factor=spec.spring_constant/(1-ratio);
    std::array<double,3> force{};
    for(std::size_t c=0;c<3;++c) { force[c]=factor*q[c]; finite(force[c],"spring force overflow"); }
    return force;
}
} // namespace
double dumbbell_spring_energy(const std::array<double,3>& q,const DumbbellSpec& spec) {
    spec.validate(); return spring_energy_validated(q,spec);
}
std::array<double,3> dumbbell_spring_force(const std::array<double,3>& q,const DumbbellSpec& spec) {
    spec.validate(); return spring_force_validated(q,spec);
}
double protocol_shear_rate(const DumbbellSpec& spec,double time) {
    finite(time,"protocol time must be finite");
    if(spec.protocol==FlowProtocol::SimpleShear) return spec.shear_rate;
    if(spec.protocol==FlowProtocol::OscillatoryShear) {
        const double phase=spec.angular_frequency*time; finite(phase,"oscillatory phase overflow");
        return spec.strain_amplitude*spec.angular_frequency*std::cos(phase);
    }
    return 0;
}
DumbbellSimulation::DumbbellSimulation(DumbbellSpec spec,std::vector<double> connectors)
    : spec_(validated(std::move(spec))),explicit_initial_state_(!connectors.empty()),
      kernels_(spec_.execution_backend=="serial"?nullptr:make_execution_kernels(spec_.execution_backend,spec_.backend_threads)),
      noise_(spec_.seed,spec_.noise_policy) {
    if(connectors.empty()) connectors.assign(spec_.n_dumbbells*3,0);
    check_values(connectors,spec_.n_dumbbells);
    // Validate the undeformed physical state as well as the affinely deformed state.
    (void)evaluate(connectors,0);
    if(spec_.protocol==FlowProtocol::StepShear) for(std::size_t i=0;i<spec_.n_dumbbells;++i) connectors[i*3]+=spec_.step_strain*connectors[i*3+1];
    state_=evaluate(std::move(connectors),0);
}
DumbbellFrame DumbbellSimulation::evaluate(std::vector<double> connectors,std::uint64_t step) const {
    DumbbellFrame frame; frame.connectors=std::move(connectors); frame.spring_forces.resize(frame.connectors.size());
    frame.step=step; frame.time=static_cast<double>(step)*spec_.dt; frame.shear_rate=protocol_shear_rate(spec_,frame.time);
    const double inverse=1.0/static_cast<double>(spec_.n_dumbbells);
    KernelResult bulk;
    if(kernels_) {
        // Validate the physical float64 domain independently of a device's
        // representation before accepting any candidate state.
        for(std::size_t i=0;i<spec_.n_dumbbells;++i)
            (void)spring_ratio({{frame.connectors[3*i],frame.connectors[3*i+1],frame.connectors[3*i+2]}},spec_);
        bulk=kernels_->springs(frame.connectors,spec_.spring_constant,spec_.q0,spec_.spring_law==SpringLaw::FENE);
        require(bulk.forces.size()==frame.connectors.size()&&bulk.energies.size()==spec_.n_dumbbells,
                "backend spring output shape invalid");
    }
    for(std::size_t i=0;i<spec_.n_dumbbells;++i) {
        const std::array<double,3> q{{frame.connectors[3*i],frame.connectors[3*i+1],frame.connectors[3*i+2]}};
        const auto force=kernels_?std::array<double,3>{{bulk.forces[3*i],bulk.forces[3*i+1],bulk.forces[3*i+2]}}:spring_force_validated(q,spec_);
        frame.mean_energy+=(kernels_?bulk.energies[i]:spring_energy_validated(q,spec_))*inverse;
        for(std::size_t a=0;a<3;++a) {
            frame.spring_forces[i*3+a]=force[a];
            for(std::size_t b=0;b<3;++b) {
                frame.conformation[a*3+b]+=q[a]*q[b]*inverse;
                frame.polymer_stress[a*3+b]+=q[a]*force[b]*inverse;
            }
        }
    }
    const double kT=spec_.diffusion/spec_.mobility;
    // This is extra stress per dumbbell, with its isotropic thermal part
    // removed. A caller still needs number density and any solvent contribution
    // before interpreting it as stress/viscosity of a suspension.
    for(std::size_t a=0;a<3;++a) frame.polymer_stress[a*3+a]-=kT;
    finite(frame.mean_energy,"mean spring energy overflow");
    for(double m:frame.conformation) finite(m,"conformation overflow");
    for(double s:frame.polymer_stress) finite(s,"polymer stress overflow");
    return frame;
}
void DumbbellSimulation::advance_unlocked(std::uint64_t steps) {
    horizon(state_.step,steps,spec_.dt);
    // q is the difference of two independently diffusing bead positions:
    // Var(dq_random)=4*D*dt and relaxation drift=-2*mu*F. Reusing the single-bead
    // factors would give the wrong connector equilibrium and relaxation time.
    const double amplitude=std::sqrt(4*spec_.diffusion*spec_.dt);
    for(std::uint64_t i=0;i<steps;++i) {
        auto next_noise=noise_; auto connectors=state_.connectors;
        const double rate=protocol_shear_rate(spec_,state_.time);
        if(kernels_) {
            std::vector<double> z(connectors.size(),0);
            if(amplitude>0) for(auto& x:z) x=next_noise.standardized();
            connectors=kernels_->linear_step(state_.connectors,state_.spring_forces,z,
                                              -2*spec_.mobility*spec_.dt,amplitude);
            require(connectors.size()==state_.connectors.size(),"backend connector output shape invalid");
        }
        for(std::size_t j=0;j<spec_.n_dumbbells;++j) for(std::size_t c=0;c<3;++c) {
            const auto k=j*3+c;
            const double affine=c==0 ? rate*state_.connectors[j*3+1]*spec_.dt:0;
            connectors[k]+=kernels_?affine:affine-2*spec_.mobility*state_.spring_forces[k]*spec_.dt+
                (amplitude>0 ? amplitude*next_noise.standardized():0);
            finite(connectors[k],"connector Euler step overflow");
        }
        auto candidate=evaluate(std::move(connectors),state_.step+1);
        // This is a per-step transaction, including the Box-Muller spare. An
        // advance(n) call may retain its successfully completed earlier steps.
        state_=std::move(candidate); noise_=std::move(next_noise);
    }
}
void DumbbellSimulation::advance(std::uint64_t steps) { std::lock_guard<std::mutex> lock(mutex_); advance_unlocked(steps); }
DumbbellFrame DumbbellSimulation::observe() const { std::lock_guard<std::mutex> lock(mutex_); return state_; }
std::vector<DumbbellFrame> DumbbellSimulation::sample(std::uint64_t steps,std::uint64_t stride,bool initial) {
    std::lock_guard<std::mutex> lock(mutex_); require(stride>0,"stride must be positive"); horizon(state_.step,steps,spec_.dt);
    std::vector<DumbbellFrame> result; if(initial) result.push_back(state_);
    for(std::uint64_t done=0;done<steps;) { const auto chunk=std::min(stride,steps-done); advance_unlocked(chunk); done+=chunk; result.push_back(state_); }
    return result;
}
} // namespace bdmc
