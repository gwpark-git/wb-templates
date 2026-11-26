#include "bdmc/reference.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void check(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
void close(double a,double b,double tol,const char* message) { if(!std::isfinite(a)||std::abs(a-b)>tol) throw std::runtime_error(std::string(message)+" actual="+std::to_string(a)+" expected="+std::to_string(b)); }
template<class F> void rejects(F f,const char* message) { bool yes=false;try{f();}catch(const std::invalid_argument&){yes=true;}check(yes,message); }
void same(const std::vector<double>& a,const std::vector<double>& b,double tol,const char* message) {check(a.size()==b.size(),message);for(std::size_t i=0;i<a.size();++i)close(a[i],b[i],tol,message);}

void repulsive_potential_and_gradient() {
    constexpr double epsilon=2.3,sigma=1.7,h=1e-6;
    for(double x:{0.03,0.2,0.7,0.98}) {
        const double r=x*sigma;
        const double derivative=(bdmc::repulsive_pair_energy(r+h,epsilon,sigma)-bdmc::repulsive_pair_energy(r-h,epsilon,sigma))/(2*h);
        close(-derivative,bdmc::repulsive_radial_force(r,epsilon,sigma),1e-9,"radial energy derivative");
    }
    close(bdmc::repulsive_pair_energy(0,epsilon,sigma),2*epsilon/3,1e-15,"contact energy");
    for(double r:{sigma,2*sigma}) {
        close(bdmc::repulsive_pair_energy(r,epsilon,sigma),0,0,"energy cutoff");
        close(bdmc::repulsive_radial_force(r,epsilon,sigma),0,0,"force cutoff");
    }
    bdmc::RepulsiveSpec s;s.n_particles=2;s.epsilon=epsilon;s.sigma=sigma;s.diffusion=0;
    bdmc::RepulsiveSimulation sim(s,{9.8,0,0,0.3,0,0});const auto f=sim.observe();
    const double force=bdmc::repulsive_radial_force(.5,epsilon,sigma);
    close(f.forces[0],-force,1e-13,"minimum image force orientation");
    close(f.forces[3],force,1e-13,"Newton antisymmetry");
    close(f.pair_virial[0],.5*force,1e-13,"positive pair virial");
    const std::vector<double> q{1,1,1,1.4,1.3,1.2};
    bdmc::RepulsiveSimulation baseline(s,q);const auto native=baseline.observe();
    for(std::size_t k=0;k<q.size();++k) {
        auto plus=q,minus=q;plus[k]+=h;minus[k]-=h;
        const double derivative=(bdmc::RepulsiveSimulation(s,plus).observe().potential_energy-bdmc::RepulsiveSimulation(s,minus).observe().potential_energy)/(2*h);
        close(native.forces[k],-derivative,1e-9,"Cartesian force gradient");
    }
    rejects([&]{bdmc::RepulsiveSimulation bad(s,{0,0,0,10,0,0});},"coincident periodic particles rejected");
}
void repulsive_cell_oracle() {
    std::mt19937_64 rng(1042);std::uniform_real_distribution<double> uniform(-20,20);
    for(std::size_t count:{2,17,83}) for(double cutoff:{.3,1.0,5.0}) {
        bdmc::RepulsiveSpec a;a.n_particles=count;a.sigma=cutoff;a.epsilon=2.1;
        auto b=a;b.force_search=bdmc::SpatialBackend::BruteForce;
        std::vector<double> positions(count*3);for(auto&x:positions)x=uniform(rng);
        bdmc::RepulsiveSimulation cells(a,positions), brute(b,positions);
        const auto ca=cells.observe(),cb=brute.observe();
        same(ca.forces,cb.forces,0,"cells vs allpairs forces");
        close(ca.potential_energy,cb.potential_energy,0,"cells vs allpairs energy");
        for(std::size_t k=0;k<9;++k)close(ca.pair_virial[k],cb.pair_virial[k],0,"cells vs allpairs virial");
        cells.advance(3);brute.advance(3);
        same(cells.observe().positions,brute.observe().positions,0,"seeded cells vs allpairs trajectory");
    }
}
void repulsive_deterministic_and_wrap() {
    bdmc::RepulsiveSpec s;s.n_particles=2;s.diffusion=0;s.mobility=.75;s.dt=.05;s.epsilon=2;
    bdmc::RepulsiveSimulation sim(s,{.001,0,0,.501,0,0});const auto old=sim.observe();sim.advance(1);const auto next=sim.observe();
    for(std::size_t k=0;k<6;++k)close(next.unwrapped_positions[k],old.unwrapped_positions[k]+s.mobility*s.dt*old.forces[k],1e-14,"simultaneous deterministic drift");
    check(next.positions[0]>9,"negative x wrapped");
    bdmc::RepulsiveSpec free;free.n_particles=1;free.dt=.2;free.diffusion=3;free.seed=8;
    bdmc::RepulsiveSimulation walk(free,{-20.2,10.1,30.4});
    walk.advance(50);const auto f=walk.observe();bdmc::OrthorhombicBox box({10,10,10});
    for(std::size_t c=0;c<3;++c) {check(f.positions[c]>=0&&f.positions[c]<10,"wrapped interval");close(f.positions[c],box.wrap(f.unwrapped_positions[c],c),2e-13,"unwrapped PBC correspondence");}
}
void free_particle_fdt() {
    for(auto policy:{bdmc::NoisePolicy::Gaussian,bdmc::NoisePolicy::LegacyUniform}) {
        bdmc::RepulsiveSpec s;s.n_particles=18000;s.box_length=100;s.sigma=.05;s.epsilon=0;s.dt=.02;s.diffusion=.7;s.mobility=.35;s.noise_policy=policy;s.seed=791;
        bdmc::RepulsiveSimulation sim(s);const auto before=sim.observe();sim.advance(1);const auto after=sim.observe();
        double mean=0,second=0,fourth=0;const auto count=before.positions.size();
        const double variance=2*s.diffusion*s.dt;
        for(std::size_t k=0;k<count;++k) {const double dx=after.unwrapped_positions[k]-before.unwrapped_positions[k];mean+=dx;second+=dx*dx;fourth+=dx*dx*dx*dx;}
        close(mean/count,0,5*std::sqrt(variance/count),"free Brownian mean");
        close(second/count,variance,.03*variance,"free FDT variance");
        close(fourth/count/(variance*variance),policy==bdmc::NoisePolicy::Gaussian?3:1.8,.17,"finite step noise fourth moment");
    }
}
void dumbbell_force_energy_stress() {
    bdmc::DumbbellSpec s;s.n_dumbbells=2;s.spring_constant=2.1;s.diffusion=.7;s.mobility=.4;s.q0=3;
    const std::array<double,3> q{{.3,-.4,.2}};const double h=1e-6;
    for(auto law:{bdmc::SpringLaw::Hookean,bdmc::SpringLaw::FENE}) {
        s.spring_law=law;const auto f=bdmc::dumbbell_spring_force(q,s);
        for(std::size_t c=0;c<3;++c) {auto p=q,m=q;p[c]+=h;m[c]-=h;const double d=(bdmc::dumbbell_spring_energy(p,s)-bdmc::dumbbell_spring_energy(m,s))/(2*h);close(f[c],d,2e-10,"spring tension is positive gradient");}
        bdmc::DumbbellSimulation sim(s,{q[0],q[1],q[2],-q[0],-q[1],-q[2]});const auto frame=sim.observe();
        for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b) {
            close(frame.conformation[a*3+b],q[a]*q[b],1e-15,"conformation outer product");
            close(frame.polymer_stress[a*3+b],q[a]*f[b]-(a==b?s.diffusion/s.mobility:0),1e-15,"positive-Cauchy stress sign");
        }
        close(frame.mean_energy,bdmc::dumbbell_spring_energy(q,s),1e-15,"energy mean per dumbbell");
    }
}
void dumbbell_protocols() {
    bdmc::DumbbellSpec s;s.n_dumbbells=1;s.dt=.01;s.diffusion=0;s.spring_constant=1.7;s.mobility=.6;
    const double a=1-2*s.mobility*s.spring_constant*s.dt;
    for(double rate:{-120.,-2.,0.,2.,120.}) {
        s.protocol=bdmc::FlowProtocol::SimpleShear;s.shear_rate=rate;
        bdmc::DumbbellSimulation sim(s,{.3,.7,-.2});sim.advance(1);const auto f=sim.observe();
        same(f.connectors,{a*.3+s.dt*rate*.7,a*.7,a*(-.2)},1e-15,"signed shear old-state Euler");
    }
    s.protocol=bdmc::FlowProtocol::StepShear;s.step_strain=-.8;
    bdmc::DumbbellSimulation step(s,{.3,.7,-.2});
    close(step.observe().connectors[0],.3-.8*.7,1e-15,"step affine initialization");
    step.advance(3);close(step.observe().connectors[0],(.3-.8*.7)*std::pow(a,3),1e-15,"step affine not repeated");
    s.protocol=bdmc::FlowProtocol::OscillatoryShear;s.strain_amplitude=-.9;s.angular_frequency=2.3;
    bdmc::DumbbellSimulation oscillatory(s,{.3,.7,-.2});oscillatory.advance(1);
    const auto one=oscillatory.observe();oscillatory.advance(1);const auto two=oscillatory.observe();
    close(two.connectors[0],a*one.connectors[0]+s.dt*s.strain_amplitude*s.angular_frequency*std::cos(s.angular_frequency*s.dt)*one.connectors[1],1e-15,"oscillatory rate at old time");
    s.protocol=bdmc::FlowProtocol::NoFlow;bdmc::DumbbellSimulation rest(s,{.3,.7,-.2});rest.advance(5);
    close(rest.observe().connectors[1],.7*std::pow(a,5),1e-15,"no flow relaxation");
}
void hookean_finite_step_moments() {
    for(auto policy:{bdmc::NoisePolicy::Gaussian,bdmc::NoisePolicy::LegacyUniform}) {
        bdmc::DumbbellSpec s;s.n_dumbbells=32000;s.dt=.012;s.mobility=.6;s.diffusion=.7;s.spring_constant=1.4;s.noise_policy=policy;s.seed=1204;
        std::vector<double> initial(s.n_dumbbells*3,.8);bdmc::DumbbellSimulation sim(s,initial);constexpr unsigned steps=25;sim.advance(steps);const auto frame=sim.observe();
        const double a=1-2*s.mobility*s.spring_constant*s.dt,mean=.8*std::pow(a,steps);
        const double var=4*s.diffusion*s.dt*(1-std::pow(a,2*steps))/(1-a*a);
        double average=0,variance=0;for(auto q:frame.connectors) {average+=q;variance+=(q-mean)*(q-mean);}const double count=frame.connectors.size();
        close(average/count,mean,5*std::sqrt(var/count),"Hookean finite step mean both noises");
        close(variance/count,var,.025*var,"Hookean finite step covariance both noises");
    }
}
void fene_domain_transaction() {
    bdmc::DumbbellSpec s;s.n_dumbbells=1;s.spring_law=bdmc::SpringLaw::FENE;s.q0=.7;s.dt=.05;
    rejects([&]{bdmc::DumbbellSimulation bad(s,{.7,0,0});},"FENE boundary initial state rejected");
    rejects([&]{bdmc::DumbbellSimulation bad(s,{.71,0,0});},"FENE exterior initial state rejected");
    s.protocol=bdmc::FlowProtocol::StepShear;s.step_strain=2;
    rejects([&]{bdmc::DumbbellSimulation bad(s,{0,.5,0});},"step shear FENE domain validated");
    s.protocol=bdmc::FlowProtocol::NoFlow;
    // Select a fixed seed for which a wrongly consumed failed-step stream
    // would make the next attempt succeed. This checks the Gaussian spare too.
    bool found=false;for(std::uint64_t seed=0;seed<1000&&!found;++seed) {
        bdmc::NoiseSource noise(seed,bdmc::NoisePolicy::Gaussian);double n1=0,n2=0;
        for(int c=0;c<3;++c)n1=std::hypot(n1,std::sqrt(4*s.diffusion*s.dt)*noise.standardized());
        for(int c=0;c<3;++c)n2=std::hypot(n2,std::sqrt(4*s.diffusion*s.dt)*noise.standardized());
        if(n1>s.q0&&n2<s.q0) {s.seed=seed;found=true;}
    }
    check(found,"fixed failed/successful trial seed exists");bdmc::DumbbellSimulation sim(s);const auto before=sim.observe();
    rejects([&]{sim.advance(1);},"FENE proposed outside state rejected");
    rejects([&]{sim.advance(1);},"failed FENE step retains RNG and cached spare");
    const auto after=sim.observe();same(after.connectors,before.connectors,0,"failed FENE state unchanged");check(after.step==0&&after.time==0,"failed FENE clock unchanged");
    s.diffusion=0;s.q0=2;s.dt=.001;bdmc::DumbbellSimulation bounded(s,{.3,.2,.1});bounded.advance(100);check(bounded.observe().mean_energy<.1,"bounded deterministic FENE relaxation");
}
void sampling_seed_and_invalids() {
    bdmc::DumbbellSpec d;d.n_dumbbells=3;d.seed=17;bdmc::DumbbellSimulation a(d),b(d);
    const auto frames=a.sample(7,3);check(frames.size()==4,"sample initial and final partial stride");
    check(frames[0].step==0&&frames[1].step==3&&frames[2].step==6&&frames[3].step==7,"sample exact steps");
    b.advance(7);same(a.observe().connectors,b.observe().connectors,0,"sampling no RNG consumption");
    check(a.sample(0,1,false).empty(),"empty sample supported");
    rejects([&]{a.sample(1,0);},"zero stride rejected");
    auto bad=d;bad.mobility=0;rejects([&]{bad.validate();},"zero mobility rejects ambiguous kT");
    bad=d;bad.dt=1;rejects([&]{bad.validate();},"unstable Hookean timestep rejected");
    bad=d;bad.diffusion=-1;rejects([&]{bad.validate();},"negative diffusion rejected");
    bad=d;bad.shear_rate=std::numeric_limits<double>::infinity();rejects([&]{bad.validate();},"nonfinite rate rejected");
    rejects([&]{bdmc::DumbbellSimulation wrong(d,{0,0});},"wrong connector shape rejected");
    bdmc::RepulsiveSpec p;p.sigma=6;rejects([&]{p.validate();},"cutoff exceeds half box rejected");
    p=bdmc::RepulsiveSpec{};p.epsilon=-1;rejects([&]{p.validate();},"negative epsilon rejected");
    rejects([]{bdmc::parse_flow_protocol("periodic_shear");},"unsupported protocol rejected");
}
} // namespace
int main() {
    const std::vector<std::pair<const char*,std::function<void()>>> tests{
        {"repulsive_potential_and_gradient",repulsive_potential_and_gradient},
        {"repulsive_cell_oracle",repulsive_cell_oracle},
        {"repulsive_deterministic_and_wrap",repulsive_deterministic_and_wrap},
        {"free_particle_fdt",free_particle_fdt},
        {"dumbbell_force_energy_stress",dumbbell_force_energy_stress},
        {"dumbbell_protocols",dumbbell_protocols},
        {"hookean_finite_step_moments",hookean_finite_step_moments},
        {"fene_domain_transaction",fene_domain_transaction},
        {"sampling_seed_and_invalids",sampling_seed_and_invalids}};
    unsigned failures=0;for(const auto& test:tests)try{test.second();std::cout<<"PASS "<<test.first<<'\n';}catch(const std::exception& e){++failures;std::cerr<<"FAIL "<<test.first<<": "<<e.what()<<'\n';}
    std::cout<<tests.size()-failures<<'/'<<tests.size()<<" reference check groups passed\n";return failures?1:0;
}
