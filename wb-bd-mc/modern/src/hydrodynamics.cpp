#include "bdmc/hydrodynamics.hpp"
#include "bdmc/reference.hpp"
#include "hi_bd.h"
#include "hi_sd.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <gsl/gsl_rng.h>

namespace bdmc {
namespace {
void require(bool condition,const char* message) {
    if(!condition) throw std::invalid_argument(message);
}
void finite(const std::vector<double>& values,const char* message) {
    for(double x:values) require(std::isfinite(x),message);
}
HydrodynamicSpec validated(HydrodynamicSpec s) { s.validate(); return s; }
HI_BD::Parameters native(const HydrodynamicSpec& s) {
    return {s.radius,s.viscosity,s.mobility==HydrodynamicMobility::Oseen ? HI_BD::Model::Oseen : HI_BD::Model::RPY};
}
void positions_shape(const std::vector<double>& q,std::size_t n) {
    require(q.size()==3*n,"HI positions must contain n_particles*3 entries");
    finite(q,"HI positions must be finite");
}
void horizon(std::uint64_t step,std::uint64_t count,double dt) {
    require(count<=std::numeric_limits<std::uint64_t>::max()-step,"HI step count overflow");
    require(std::isfinite(static_cast<double>(step+count)*dt),"HI observation time overflow");
}
}

std::string hydrodynamic_mobility_name(HydrodynamicMobility m) {
    if(m==HydrodynamicMobility::RPY) return "rpy";
    if(m==HydrodynamicMobility::Oseen) return "oseen";
    throw std::invalid_argument("Unknown hydrodynamic mobility");
}
HydrodynamicMobility parse_hydrodynamic_mobility(const std::string& m) {
    if(m=="rpy") return HydrodynamicMobility::RPY;
    if(m=="oseen") return HydrodynamicMobility::Oseen;
    throw std::invalid_argument("Hydrodynamic mobility must be rpy or oseen");
}
void HydrodynamicSpec::validate() const {
    require(n_particles>0&&n_particles<=128,"Dense Brownian HI requires 1..128 particles");
    require(std::isfinite(radius)&&radius>0&&std::isfinite(viscosity)&&viscosity>0,"HI radius/viscosity must be finite positive");
    require(std::isfinite(dt)&&dt>0&&std::isfinite(thermal_energy)&&thermal_energy>=0,"HI dt>0 and thermal_energy>=0 required");
    require(std::isfinite(epsilon)&&epsilon>=0&&std::isfinite(sigma)&&sigma>0&&std::isfinite(epsilon/sigma),"HI soft repulsion requires finite epsilon>=0 and sigma>0");
    require(std::isfinite(shear_rate),"HI shear_rate must be finite");
    require(seed<=std::numeric_limits<std::uint32_t>::max(),"GSL mt19937 seed must fit uint32 to avoid silent truncation");
    hydrodynamic_mobility_name(mobility);
    require(external_forces.empty()||external_forces.size()==3*n_particles,"HI external force shape mismatch");
    require(background_velocity.empty()||background_velocity.size()==3*n_particles,"HI background velocity shape mismatch");
    finite(external_forces,"HI external force must be finite");
    finite(background_velocity,"HI background velocity must be finite");
    // Validate representable isolated mobility and stochastic scaling now.
    HI_BD::mobility({0,0,0},native(*this));
    require(std::isfinite(std::sqrt(2*thermal_energy)*std::sqrt(dt)),"HI stochastic scale overflow");
}

std::vector<double> hydrodynamic_mobility(const std::vector<double>& q,const HydrodynamicSpec& s) {
    s.validate(); positions_shape(q,s.n_particles);
    return HI_BD::mobility(q,native(s));
}

struct HydrodynamicSimulation::RngState {
    gsl_rng* value;
    explicit RngState(std::uint64_t seed):value(gsl_rng_alloc(gsl_rng_mt19937)) {
        if(!value) throw std::bad_alloc();
        gsl_rng_set(value,static_cast<unsigned long>(seed));
    }
    RngState(const RngState& r):value(gsl_rng_clone(r.value)) { if(!value) throw std::bad_alloc(); }
    ~RngState() { gsl_rng_free(value); }
};

HydrodynamicSimulation::HydrodynamicSimulation(HydrodynamicSpec s,std::vector<double> q)
    :spec_(validated(std::move(s))),rng_(std::make_unique<RngState>(spec_.seed)),state_(evaluate(std::move(q),0)) {}
HydrodynamicSimulation::~HydrodynamicSimulation()=default;

HydrodynamicFrame HydrodynamicSimulation::evaluate(std::vector<double> q,std::uint64_t step) const {
    positions_shape(q,spec_.n_particles);
    const std::vector<double> zero(q.size(),0);
    // Validate every committed configuration, including final candidate states.
    // This extra factorization is included in end-to-end performance timings.
    HI_BD::increment(q,zero,zero,spec_.dt,0,nullptr,native(spec_));
    HydrodynamicFrame f;f.positions=std::move(q);f.step=step;f.time=static_cast<double>(step)*spec_.dt;
    f.forces=spec_.external_forces.empty()?zero:spec_.external_forces;
    if(spec_.epsilon>0) for(std::size_t i=0;i<spec_.n_particles;++i) for(std::size_t j=i+1;j<spec_.n_particles;++j) {
        std::array<double,3> delta{};
        double r=0;
        for(std::size_t k=0;k<3;++k) { delta[k]=f.positions[3*j+k]-f.positions[3*i+k];r=std::hypot(r,delta[k]); }
        require(std::isfinite(r)&&r>0,"HI pair gradient requires finite noncoincident centres");
        if(r>=spec_.sigma) continue;
        const double radial=repulsive_radial_force(r,spec_.epsilon,spec_.sigma);
        f.potential_energy+=repulsive_pair_energy(r,spec_.epsilon,spec_.sigma);
        for(std::size_t a=0;a<3;++a) {
            const double force=radial*(delta[a]/r);
            f.forces[3*i+a]-=force;f.forces[3*j+a]+=force;
            for(std::size_t b=0;b<3;++b) f.pair_virial[3*a+b]+=force*delta[b];
        }
    }
    finite(f.forces,"HI force overflow");
    require(std::isfinite(f.potential_energy),"HI potential energy overflow");
    for(double x:f.pair_virial) require(std::isfinite(x),"HI pair virial overflow");
    return f;
}
void HydrodynamicSimulation::advance_unlocked(std::uint64_t steps) {
    horizon(state_.step,steps,spec_.dt);
    for(std::uint64_t step=0;step<steps;++step) {
        // Mobility/noise uses the full unwrapped free-space configuration.
        // Clone GSL state until the candidate also passes mobility/domain
        // checks; a failed step must preserve both geometry and random stream.
        auto rng=std::make_unique<RngState>(*rng_);
        auto background=spec_.background_velocity;
        if(background.empty()) background.assign(3*spec_.n_particles,0);
        for(std::size_t i=0;i<spec_.n_particles;++i) background[3*i]+=spec_.shear_rate*state_.positions[3*i+1];
        const auto dq=HI_BD::increment(state_.positions,state_.forces,background,spec_.dt,spec_.thermal_energy,rng->value,native(spec_));
        auto positions=state_.positions;
        for(std::size_t k=0;k<positions.size();++k) positions[k]+=dq[k];
        auto candidate=evaluate(std::move(positions),state_.step+1);
        state_=std::move(candidate);rng_=std::move(rng);
    }
}
void HydrodynamicSimulation::advance(std::uint64_t steps) {
    std::lock_guard<std::mutex> lock(mutex_);advance_unlocked(steps);
}
HydrodynamicFrame HydrodynamicSimulation::observe() const {
    std::lock_guard<std::mutex> lock(mutex_);return state_;
}
std::vector<HydrodynamicFrame> HydrodynamicSimulation::sample(std::uint64_t steps,std::uint64_t stride,bool initial) {
    std::lock_guard<std::mutex> lock(mutex_);
    require(stride>0,"HI sample stride must be positive");horizon(state_.step,steps,spec_.dt);
    const std::uint64_t intervals=steps/stride+(steps%stride!=0);
    require(!initial||intervals<std::numeric_limits<std::uint64_t>::max(),"HI sample frame count overflow");
    const std::uint64_t count=intervals+static_cast<unsigned>(initial);
    const std::uint64_t frame_bytes=6*spec_.n_particles*sizeof(double)+sizeof(HydrodynamicFrame);
    require(count<=256*1024*1024/frame_bytes,"HI sample exceeds 256MiB in-memory frame bound; sample in batches");
    std::vector<HydrodynamicFrame> frames;frames.reserve(static_cast<std::size_t>(count));
    if(initial) frames.push_back(state_);
    while(steps>0) { const auto n=std::min(steps,stride);advance_unlocked(n);steps-=n;frames.push_back(state_); }
    return frames;
}

void StokesianSpec::validate() const {
    require(std::isfinite(radius)&&radius>0&&std::isfinite(viscosity)&&viscosity>0,"SD radius/viscosity must be finite positive");
    require(std::isfinite(lubrication_cutoff)&&lubrication_cutoff>2&&lubrication_cutoff<=3.99,"SD lubrication cutoff must be in (2,3.99]");
}
StokesianResponse solve_stokesian(const std::vector<double>& q,const std::vector<double>& f,
    const std::vector<double>& t,const std::vector<double>& u,const std::vector<double>& w,
    const std::vector<double>& e,const StokesianSpec& s) {
    s.validate();const auto r=HI_SD::solve(q,f,t,u,w,e,{s.radius,s.viscosity,s.lubrication,s.lubrication_cutoff});
    return {r.translation,r.rotation,r.stresslet,r.relative_residual,r.minimum_gap_over_radius};
}
} // namespace bdmc
