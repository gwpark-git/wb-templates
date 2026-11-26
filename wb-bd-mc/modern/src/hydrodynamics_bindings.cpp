#include "bdmc/hydrodynamics.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>

namespace py=pybind11;
namespace {
using Array=py::array_t<double,py::array::c_style|py::array::forcecast>;
std::vector<double> vectors(py::object object,std::size_t n,const char* message,bool optional=false) {
    if(optional&&object.is_none()) return std::vector<double>(3*n,0);
    auto a=Array::ensure(object);
    if(!a||a.ndim()!=2||a.shape(0)!=static_cast<py::ssize_t>(n)||a.shape(1)!=3) throw std::invalid_argument(message);
    return {a.data(),a.data()+a.size()};
}
std::vector<double> positions(py::object object) {
    auto a=Array::ensure(object);
    if(!a||a.ndim()!=2||a.shape(1)!=3||a.shape(0)<1||a.shape(0)>128) throw std::invalid_argument("positions must have shape (N,3), 1<=N<=128");
    return {a.data(),a.data()+a.size()};
}
template<class Values> py::array_t<double> copied(const Values& v,std::vector<py::ssize_t> shape) {
    py::array_t<double> a(shape);std::copy(v.begin(),v.end(),a.mutable_data());return a;
}
py::dict spec_dict(const bdmc::HydrodynamicSpec& s) {
    py::dict d;d["n_particles"]=s.n_particles;d["radius"]=s.radius;d["viscosity"]=s.viscosity;
    d["dt"]=s.dt;d["thermal_energy"]=s.thermal_energy;d["mobility"]=bdmc::hydrodynamic_mobility_name(s.mobility);
    d["epsilon"]=s.epsilon;d["sigma"]=s.sigma;d["shear_rate"]=s.shear_rate;d["seed"]=s.seed;
    d["external_forces"]=s.external_forces;d["background_velocity"]=s.background_velocity;return d;
}
py::dict metadata(const bdmc::HydrodynamicSpec& s) {
    auto d=spec_dict(s);d["schema"]="bdmc.hydrodynamics.observation.v1";
    d["model"]="free_space_brownian_hi";d["boundary"]="unbounded free space; no periodic wrapping or minimum-image forces";
    d["execution_backend"]="serial_dense_gsl_reference";d["dense_particle_limit"]=128;
    d["integrator"]="old-state Euler-Maruyama: dq=(u+M*F)*dt+sqrt(2*kBT*dt)*L*z; L*L^T=M";
    d["flow"]="prescribed u_x+=shear_rate*y plus optional constant per-particle ambient velocity; no solved solvent";
    d["thermal_drift"]="div_q(M)=0 for this equal-radius free-space translational tensor only";
    d["units"]="one caller-selected coherent length/time/energy system; radius independent of soft-force cutoff sigma";
    d["potential"]="U=epsilon/3*(1-r/sigma)^2*(2+r/sigma) for r<sigma, zero outside; epsilon direct energy";
    d["force"]="F=-grad(U); radial epsilon/sigma*(1-(r/sigma)^2); external forces additive";
    d["pair_virial"]="sum unordered pairs r_ij outer F_on_j; unnormalised soft-force virial only";
    d["rng"]="GSL mt19937; Gaussian ziggurat; uint32 seed; seed zero selects GSL default 4357";
    d["rng_compatibility"]="matches frozen native HI stream for the same inputs; distinct from modern std::mt19937_64";
    d["array_ownership"]="independent owned C-contiguous float64 NumPy copies";
    d["failure_contract"]="failed individual step preserves state and RNG; prior successful batch steps remain completed; no retry or clipping";
    d["configuration_checks"]="dense positive-definiteness checked at construction and before each candidate commit; no jitter";
    d["observation_phase"]="current complete state; time=absolute_step*dt; sample includes final partial stride";
    d["restart_contract"]="not implemented; observation arrays are not restart checkpoints";
    d["scientific_scope"]="small equal-sphere free-space translational HI reference; no lubrication, stresslet constraint, walls, PBC, network MC or large-N accelerator";
    return d;
}
py::dict frame(const bdmc::HydrodynamicFrame& f,const bdmc::HydrodynamicSpec& s) {
    py::dict d;const std::vector<py::ssize_t> shape{static_cast<py::ssize_t>(s.n_particles),3};
    d["positions"]=copied(f.positions,shape);d["unwrapped_positions"]=copied(f.positions,shape);
    d["forces"]=copied(f.forces,shape);d["potential_energy"]=f.potential_energy;
    d["pair_virial"]=copied(f.pair_virial,{3,3});d["step"]=f.step;d["time"]=f.time;d["metadata"]=metadata(s);return d;
}
py::dict batch(const std::vector<bdmc::HydrodynamicFrame>& fs,const bdmc::HydrodynamicSpec& s) {
    const auto n=static_cast<py::ssize_t>(s.n_particles),count=static_cast<py::ssize_t>(fs.size());
    const std::vector<py::ssize_t> shape{count,n,3};
    py::array_t<double> q(shape),f(shape),virial({count,py::ssize_t{3},py::ssize_t{3}}),energy(count),times(count);
    py::array_t<std::uint64_t> steps(count);
    for(std::size_t i=0;i<fs.size();++i) {
        std::copy(fs[i].positions.begin(),fs[i].positions.end(),q.mutable_data()+3*s.n_particles*i);
        std::copy(fs[i].forces.begin(),fs[i].forces.end(),f.mutable_data()+3*s.n_particles*i);
        std::copy(fs[i].pair_virial.begin(),fs[i].pair_virial.end(),virial.mutable_data()+9*i);
        energy.mutable_data()[i]=fs[i].potential_energy;times.mutable_data()[i]=fs[i].time;steps.mutable_data()[i]=fs[i].step;
    }
    py::array_t<double> unwrapped(shape);
    if(q.size()) std::copy(q.data(),q.data()+q.size(),unwrapped.mutable_data());
    py::dict d;d["positions"]=q;d["unwrapped_positions"]=unwrapped;
    d["forces"]=f;d["pair_virial"]=virial;d["potential_energy"]=energy;d["times"]=times;d["steps"]=steps;d["metadata"]=metadata(s);return d;
}
}

PYBIND11_MODULE(_hydrodynamics,m) {
    m.doc()="Frozen native free-space Brownian HI and athermal Stokesian response with owned arrays";
    py::class_<bdmc::HydrodynamicSpec>(m,"HydrodynamicSpec")
      .def(py::init([](std::size_t n,double radius,double eta,double dt,double kt,const std::string& mobility,double epsilon,double sigma,double shear,std::uint64_t seed,py::object force,py::object ambient) {
        bdmc::HydrodynamicSpec s;s.n_particles=n;s.radius=radius;s.viscosity=eta;s.dt=dt;s.thermal_energy=kt;
        s.mobility=bdmc::parse_hydrodynamic_mobility(mobility);s.epsilon=epsilon;s.sigma=sigma;s.shear_rate=shear;s.seed=seed;
        s.validate();s.external_forces=vectors(force,n,"external_forces must have shape (N,3)",true);
        s.background_velocity=vectors(ambient,n,"background_velocity must have shape (N,3)",true);s.validate();return s;
      }),py::arg("n_particles")=1,py::arg("radius")=1.,py::arg("viscosity")=1.,py::arg("dt")=.001,py::arg("thermal_energy")=1.,py::arg("mobility")="rpy",py::arg("epsilon")=0.,py::arg("sigma")=1.,py::arg("shear_rate")=0.,py::arg("seed")=0,py::arg("external_forces")=py::none(),py::arg("background_velocity")=py::none())
      .def_readwrite("n_particles",&bdmc::HydrodynamicSpec::n_particles)
      .def_readwrite("radius",&bdmc::HydrodynamicSpec::radius).def_readwrite("viscosity",&bdmc::HydrodynamicSpec::viscosity)
      .def_readwrite("dt",&bdmc::HydrodynamicSpec::dt).def_readwrite("thermal_energy",&bdmc::HydrodynamicSpec::thermal_energy)
      .def_property("mobility",[](const bdmc::HydrodynamicSpec& s){return bdmc::hydrodynamic_mobility_name(s.mobility);},[](bdmc::HydrodynamicSpec& s,const std::string& v){s.mobility=bdmc::parse_hydrodynamic_mobility(v);})
      .def_readwrite("epsilon",&bdmc::HydrodynamicSpec::epsilon).def_readwrite("sigma",&bdmc::HydrodynamicSpec::sigma)
      .def_readwrite("shear_rate",&bdmc::HydrodynamicSpec::shear_rate).def_readwrite("seed",&bdmc::HydrodynamicSpec::seed)
      .def_readwrite("external_forces",&bdmc::HydrodynamicSpec::external_forces).def_readwrite("background_velocity",&bdmc::HydrodynamicSpec::background_velocity)
      .def("validate",&bdmc::HydrodynamicSpec::validate).def("to_dict",&spec_dict);
    py::class_<bdmc::HydrodynamicSimulation>(m,"HydrodynamicSimulation")
      .def(py::init([](bdmc::HydrodynamicSpec s,py::object q){s.validate();return std::make_unique<bdmc::HydrodynamicSimulation>(s,vectors(q,s.n_particles,"initial_positions must have shape (N,3)"));}),py::arg("spec"),py::arg("initial_positions"))
      .def_property_readonly("spec",&bdmc::HydrodynamicSimulation::spec)
      .def_property_readonly("metadata",[](const bdmc::HydrodynamicSimulation& s){return metadata(s.spec());})
      .def("advance",&bdmc::HydrodynamicSimulation::advance,py::arg("steps"),py::call_guard<py::gil_scoped_release>())
      .def("observe",[](const bdmc::HydrodynamicSimulation& s){auto f=[&](){py::gil_scoped_release release;return s.observe();}();return frame(f,s.spec());})
      .def("sample",[](bdmc::HydrodynamicSimulation& s,std::uint64_t steps,std::uint64_t stride,bool initial){auto fs=[&](){py::gil_scoped_release release;return s.sample(steps,stride,initial);}();return batch(fs,s.spec());},py::arg("steps"),py::arg("stride")=1,py::arg("include_initial")=true);
    m.def("mobility",[](py::object q,const bdmc::HydrodynamicSpec& s){auto pos=vectors(q,s.n_particles,"positions must have shape (N,3)");auto values=[&](){py::gil_scoped_release release;return bdmc::hydrodynamic_mobility(pos,s);}();const auto dim=static_cast<py::ssize_t>(3*s.n_particles);return copied(values,{dim,dim});},py::arg("positions"),py::arg("spec"));
    py::class_<bdmc::StokesianSpec>(m,"StokesianSpec")
      .def(py::init([](double radius,double eta,bool lub,double cutoff){bdmc::StokesianSpec s{radius,eta,lub,cutoff};s.validate();return s;}),py::arg("radius")=1.,py::arg("viscosity")=1.,py::arg("lubrication")=true,py::arg("lubrication_cutoff")=3.99)
      .def_readwrite("radius",&bdmc::StokesianSpec::radius).def_readwrite("viscosity",&bdmc::StokesianSpec::viscosity)
      .def_readwrite("lubrication",&bdmc::StokesianSpec::lubrication).def_readwrite("lubrication_cutoff",&bdmc::StokesianSpec::lubrication_cutoff)
      .def("validate",&bdmc::StokesianSpec::validate);
    m.def("solve_stokesian",[](py::object q,const bdmc::StokesianSpec& s,py::object f,py::object t,py::object u,py::object w,py::object e) {
        auto pos=positions(q);const auto n=pos.size()/3;
        if(n>64) throw std::invalid_argument("Dense Stokesian solve requires 1..64 particles");
        auto force=vectors(f,n,"forces must have shape (N,3)",true),torque=vectors(t,n,"torques must have shape (N,3)",true);
        auto ambient=vectors(u,n,"background_velocity must have shape (N,3)",true),angular=vectors(w,n,"background_angular_velocity must have shape (N,3)",true);
        std::vector<double> strain(9,0);
        if(!e.is_none()) { auto a=Array::ensure(e);if(!a||!((a.ndim()==2&&a.shape(0)==3&&a.shape(1)==3)||(a.ndim()==3&&a.shape(0)==static_cast<py::ssize_t>(n)&&a.shape(1)==3&&a.shape(2)==3))) throw std::invalid_argument("strain must have shape (3,3) or (N,3,3)");strain.assign(a.data(),a.data()+a.size()); }
        auto r=[&](){py::gil_scoped_release release;return bdmc::solve_stokesian(pos,force,torque,ambient,angular,strain,s);}();
        const auto count=static_cast<py::ssize_t>(n);py::dict d;
        d["translation"]=copied(r.translation,{count,3});d["rotation"]=copied(r.rotation,{count,3});d["stresslet"]=copied(r.stresslet,{count,3,3});
        d["relative_residual"]=r.relative_residual;d["minimum_gap_over_radius"]=r.minimum_gap_over_radius;
        py::dict info;info["schema"]="bdmc.stokesian.response.v1";info["boundary"]="unbounded free space; no periodic images or walls";
        info["thermal"]="athermal only; no Brownian increment or thermal drift";info["scope"]="instantaneous translation/rotation/stresslet; no time integration, network coupling or complete suspension stress";
        info["radius"]=s.radius;info["viscosity"]=s.viscosity;info["lubrication"]=s.lubrication;info["lubrication_cutoff"]=s.lubrication_cutoff;
        info["ambient"]="supplied centre translation/angular velocity and symmetric traceless strain; no solved fluid";
        info["array_ownership"]="independent owned C-contiguous float64 NumPy copies";info["dense_particle_limit"]=64;
        info["source"]="vendored dev_HI native invariant-tensor port; pinned JFSD source/table provenance in HI_SOURCE_MANIFEST.json";
        d["metadata"]=info;return d;
    },py::arg("positions"),py::arg("spec")=bdmc::StokesianSpec{},py::arg("forces")=py::none(),py::arg("torques")=py::none(),py::arg("background_velocity")=py::none(),py::arg("background_angular_velocity")=py::none(),py::arg("strain")=py::none());
}
