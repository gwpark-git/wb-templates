#include "bdmc/hydrodynamics.hpp"
#include "hi_bd.h"
#include "hi_sd.h"
#include <gsl/gsl_rng.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void check(bool ok,const char* label) { if(!ok) throw std::runtime_error(label); }
bool near(double a,double b,double tol=2e-13) {return std::abs(a-b)<=tol;}
void rejects(const std::function<void()>& f) {bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw,"expected rejection");}
}
int main() try {
    using D=std::vector<double>;
    const double pi=std::acos(-1.);
    bdmc::HydrodynamicSpec s;s.n_particles=2;s.radius=1;s.viscosity=1/(6*pi);s.dt=.002;s.thermal_energy=1.7;s.seed=810729;
    const D q{0,0,0,4,0,0},zero(6,0);
    const auto matrix=bdmc::hydrodynamic_mobility(q,s);
    check(near(matrix[3],.359375)&&near(matrix[10],.1953125),"RPY literal pair");
    s.mobility=bdmc::HydrodynamicMobility::Oseen;
    check(near(bdmc::hydrodynamic_mobility(q,s)[3],.375),"Oseen literal pair");
    s.mobility=bdmc::HydrodynamicMobility::RPY;
    check(near(bdmc::hydrodynamic_mobility({0,0,0,1,0,0},s)[3],.8125),"overlap RPY literal");
    bdmc::HydrodynamicSimulation modern(s,q),replay(s,q);
    gsl_rng* rng=gsl_rng_alloc(gsl_rng_mt19937);gsl_rng_set(rng,810729UL);
    auto native_q=q;
    for(int step=0;step<7;++step) {
        const auto dq=HI_BD::increment(native_q,zero,zero,s.dt,s.thermal_energy,rng,{s.radius,s.viscosity,HI_BD::Model::RPY});
        for(std::size_t i=0;i<6;++i) native_q[i]+=dq[i];
        modern.advance(1);
    }
    gsl_rng_free(rng);
    check(modern.observe().positions==native_q,"exact frozen native HI path parity");
    replay.advance(7);check(replay.observe().positions==native_q,"same seed replay");
    const auto frames=modern.sample(5,2,true);
    check(frames.size()==4&&frames[0].step==7&&frames[1].step==9&&frames[3].step==12,"batch timing partial stride");
    check(near(frames.back().time,12*s.dt),"physical time");
    const auto before=modern.observe();rejects([&]{modern.sample(1,0);});
    rejects([&]{modern.sample(std::numeric_limits<std::uint64_t>::max(),1);});
    check(modern.observe().positions==before.positions&&modern.observe().step==before.step,"invalid batch transaction");
    auto one=s;one.n_particles=1;one.thermal_energy=0;one.external_forces={2,0,0};one.background_velocity={.2,0,0};one.shear_rate=.3;
    bdmc::HydrodynamicSimulation driven(one,{0,2,0});driven.advance(1);
    check(near(driven.observe().positions[0],one.dt*(2+.2+.6)),"Stokes drag and coherent shear");
    auto soft=s;soft.thermal_energy=0;soft.epsilon=6;soft.sigma=8;
    bdmc::HydrodynamicSimulation repulsive(soft,q);auto force=repulsive.observe();
    check(near(force.forces[0],-.5625)&&near(force.forces[3],.5625),"direct energy soft force scaling");
    check(near(force.potential_energy,1.25)&&near(force.pair_virial[0],2.25),"soft potential and pair virial");
    repulsive.advance(1);check(repulsive.observe().positions[0]<0&&repulsive.observe().positions[3]>4,"free-space repulsion direction");
    rejects([&]{bdmc::HydrodynamicSimulation bad(s,{0,0,0,0,0,0});});
    rejects([&]{auto bad=s;bad.n_particles=129;bad.validate();});
    rejects([&]{auto bad=s;bad.seed=std::numeric_limits<std::uint64_t>::max();bad.validate();});
    rejects([&]{auto bad=s;bad.sigma=0;bad.validate();});
    rejects([&]{bdmc::hydrodynamic_mobility({0,0,0},s);});
    rejects([&]{auto bad=s;bad.external_forces={1,2};bad.validate();});
    rejects([&]{auto bad=s;bad.radius=std::numeric_limits<double>::quiet_NaN();bad.validate();});
    // Candidate overlap is rejected after random draws. Repeated retries must
    // reproduce that rejected trial instead of resampling until success.
    auto trial=s;trial.mobility=bdmc::HydrodynamicMobility::Oseen;trial.seed=3519;trial.dt=.1;trial.thermal_energy=1;
    bdmc::HydrodynamicSimulation failure(trial,{0,0,0,2.000001,0,0});
    const auto initial=failure.observe();
    bool rejected=false;try{failure.advance(1);}catch(const std::domain_error&){rejected=true;}
    check(rejected,"fixed seed candidate overlap rejected");
    for(int i=0;i<16;++i) rejects([&]{failure.advance(1);});
    check(failure.observe().positions==initial.positions&&failure.observe().step==0,"failed candidate preserves RNG and state");

    bdmc::StokesianSpec sd;sd.radius=.7;sd.viscosity=1.3;
    const D force1{1,2,3},torque1{2,-1,.5},ambient{.2,.3,.4},omega{.1,-.2,.3},strain{1,.2,0,.2,-.3,0,0,0,-.7};
    const auto response=bdmc::solve_stokesian({0,0,0},force1,torque1,ambient,omega,strain,sd);
    for(std::size_t i=0;i<3;++i) {
        check(near(response.translation[i],ambient[i]+force1[i]/(6*pi*sd.viscosity*sd.radius)),"isolated SD Stokes drag");
        check(near(response.rotation[i],omega[i]+torque1[i]/(8*pi*sd.viscosity*std::pow(sd.radius,3))),"isolated SD rotational drag");
    }
    for(std::size_t i=0;i<9;++i) check(near(response.stresslet[i],(20*pi/3)*sd.viscosity*std::pow(sd.radius,3)*strain[i]),"isolated Einstein stresslet");
    const D pair{0,0,0,1.5,0,0},forces{1,0,0,-1,0,0},e0(9,0);
    const auto direct=HI_SD::solve(pair,forces,zero,zero,zero,e0,{sd.radius,sd.viscosity,true,3.99});
    const auto wrapped=bdmc::solve_stokesian(pair,forces,zero,zero,zero,e0,sd);
    check(direct.translation==wrapped.translation&&direct.rotation==wrapped.rotation&&direct.stresslet==wrapped.stresslet,"frozen native SD parity");
    check(wrapped.relative_residual<1e-9,"SD residual");
    rejects([&]{bdmc::solve_stokesian({0,0,0,1.40001,0,0},forces,zero,zero,zero,e0,sd);});
    rejects([&]{bdmc::solve_stokesian({0,0,0},force1,torque1,ambient,omega,D{0,1,0,0,0,0,0,0,0},sd);});
    std::cout<<"{\"status\":\"passed\",\"scope\":\"modern HI bridge: analytic limits, frozen-source parity, state/RNG transactions, free-space forces, copied observations checked separately\"}\n";
    return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
