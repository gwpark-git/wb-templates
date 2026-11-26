#include "bdmc/reference.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <limits>
#include <memory>

namespace py = pybind11;
namespace {
py::ssize_t size(std::size_t n) {
    if(n>static_cast<std::size_t>(std::numeric_limits<py::ssize_t>::max())) throw std::overflow_error("array shape exceeds ssize_t");
    return static_cast<py::ssize_t>(n);
}
template<class Range> py::array_t<double> copied(const Range& values,std::vector<py::ssize_t> shape) {
    py::array_t<double> output(shape); std::copy(values.begin(),values.end(),output.mutable_data()); return output;
}
std::vector<double> initial(py::object object,std::size_t n,const char* message) {
    if(object.is_none()) return {};
    const auto a=py::array_t<double,py::array::c_style|py::array::forcecast>::ensure(object);
    if(!a || a.ndim()!=2 || a.shape(0)!=size(n) || a.shape(1)!=3) throw std::invalid_argument(message);
    return {a.data(),a.data()+a.size()};
}
py::dict spec_dict(const bdmc::RepulsiveSpec& s) {
    py::dict d;
    d["execution_backend"]=s.execution_backend;d["backend_threads"]=s.backend_threads;
    d["model"]="repulsive_brownian"; d["dimension"]=3; d["n_particles"]=s.n_particles;
    d["box_length"]=s.box_length; d["dt"]=s.dt; d["mobility"]=s.mobility; d["diffusion"]=s.diffusion;
    d["epsilon"]=s.epsilon; d["sigma"]=s.sigma; d["seed"]=s.seed;
    d["noise_policy"]=bdmc::noise_name(s.noise_policy); d["force_search"]=bdmc::spatial_backend_name(s.force_search);
    return d;
}
py::dict spec_dict(const bdmc::DumbbellSpec& s) {
    py::dict d;
    d["execution_backend"]=s.execution_backend;d["backend_threads"]=s.backend_threads;
    d["model"]="kinetic_dumbbell"; d["dimension"]=3; d["n_dumbbells"]=s.n_dumbbells;
    d["dt"]=s.dt; d["mobility"]=s.mobility; d["diffusion"]=s.diffusion; d["spring_constant"]=s.spring_constant;
    d["spring_law"]=bdmc::spring_law_name(s.spring_law); d["q0"]=s.q0;
    d["protocol"]=bdmc::flow_protocol_name(s.protocol); d["shear_rate"]=s.shear_rate;
    d["step_strain"]=s.step_strain; d["strain_amplitude"]=s.strain_amplitude;
    d["angular_frequency"]=s.angular_frequency; d["seed"]=s.seed; d["noise_policy"]=bdmc::noise_name(s.noise_policy);
    return d;
}
template<class Spec> py::dict common_metadata(const Spec& s) {
    // The exchange contract includes units, observation phase and numerical
    // policy so a learned dataset cannot infer physics from a class name alone.
    // Keep these declarations aligned when adding a model or execution backend.
    auto d=spec_dict(s);
    d["schema"]="bdmc.reference.observation.v1";
    d["execution_backend"]=s.execution_backend; d["backend_threads"]=s.backend_threads;
    d["arithmetic_precision"]=s.execution_backend=="metal" ? "float32 bulk kernels; float64 host storage/reductions" : "float64";
    d["backend"]=s.execution_backend=="serial" ? "modern_cpu_reference" : "modern_"+s.execution_backend;
    d["execution"]="explicit bulk kernel backend; host RNG and state/observable orchestration; MPI independent jobs are a separate layer; no HI";
    d["thermal_energy"]=s.diffusion/s.mobility;
    d["fluctuation_dissipation"]="D=mobility*kT; bead coefficients; kT=diffusion/mobility";
    d["time_unit"]="caller-selected physical/reduced time; legacy clock requires an explicit mapping";
    d["energy_unit"]="direct epsilon or spring_constant energy units; no implicit temperature scaling";
    d["rng"]="std::mt19937_64; top-53-bit uniform; Box-Muller Gaussian with cached spare";
    d["backend_reproducibility"]="fixed host RNG draw order; backend arithmetic may change paths; Metal float32 requires statistical/tolerance checks";
    d["rng_compatibility"]="modern seed stream; not legacy GSL trajectory identity; no cross-platform bitwise promise";
    d["initialization_rng"]="initialization consumes no RNG";
    d["noise_finite_step_law"]=s.noise_policy==bdmc::NoisePolicy::Gaussian ? "standard Gaussian z" : "standardized uniform z in [-sqrt(3),sqrt(3)); variance matched, not Gaussian finite-step law";
    d["array_ownership"]="independent owned C-contiguous float64 NumPy copies";
    d["observation_phase"]="current complete state; after completed steps; optional initial frame";
    d["observation_time"]="time=absolute_step*dt; stride relative to sample call; final partial stride included";
    d["restart_contract"]="unsupported for reference engines; observations are not restart checkpoints";
    d["failure_contract"]="failed individual step preserves state and RNG; earlier successful steps in a batch remain completed";
    d["scientific_scope"]="standard kinetic reference equations; focused checks are not historical regression or complete scientific validation";
    return d;
}
py::dict metadata(const bdmc::RepulsiveSpec& s,bool explicit_state) {
    auto d=common_metadata(s);
    d["boundary"]="3D cubic periodic; minimum image in [-L/2,L/2); sigma<=L/2";
    d["initialization"]=explicit_state ? "explicit supplied positions copied; unwrapped retained" : "deterministic regular grid; not an equilibrated ensemble";
    d["potential"]="U(r)=epsilon/3*(1-r/sigma)^2*(2+r/sigma) for r<sigma; zero otherwise; each unordered pair once";
    d["force_convention"]="F=-gradient U; outward radial magnitude epsilon/sigma*(1-(r/sigma)^2)";
    d["coincidence_policy"]="coincident distinct periodic particles rejected: radial potential has no unique Cartesian gradient";
    d["integrator"]="simultaneous old-state Euler: dx=mobility*F*dt+sqrt(2*diffusion*dt)*z";
    d["noise_increment_variance"]="2*diffusion*dt per particle coordinate";
    d["pair_virial_convention"]="sum unordered pairs r_ij outer F_on_j; positive repulsive tensor; no volume normalization, solvent or hydrodynamic stress";
    d["axis_order"]="sample: frame,particle,coordinate; tensor frame,row,column";
    d["legacy_mapping"]="soft-potential shape preserved; direct epsilon energy and explicit mobility/diffusion clock must be mapped before legacy comparison";
    return d;
}
py::dict metadata(const bdmc::DumbbellSpec& s,bool explicit_state) {
    auto d=common_metadata(s);
    d["boundary"]="unbounded connector space; independent kinetic dumbbells; no periodic bead box";
    d["initialization"]=explicit_state ? "explicit supplied connectors copied" : "zero connectors; not an equilibrated ensemble";
    d["step_shear_timing"]="if selected, qx+=step_strain*qy applied once during construction before time-zero observation";
    d["flow_convention"]="K_xy=dgamma/dt; x flow, y gradient; simple rate constant; oscillatory rate=strain_amplitude*angular_frequency*cos(angular_frequency*t_old)";
    d["integrator"]="simultaneous old-state Euler: dq=(K*q-2*mobility*F_spring)*dt+sqrt(4*diffusion*dt)*z";
    d["noise_increment_variance"]="4*diffusion*dt per connector coordinate";
    d["spring_force_convention"]="positive connector tension F_spring=H*q (Hookean) or H*q/(1-|q|^2/q0^2) (FENE); H=spring_constant";
    d["spring_energy"]="Hookean U=H*|q|^2/2; FENE U=-H*q0^2/2*log(1-|q|^2/q0^2); mean per dumbbell";
    d["polymer_stress_convention"]="positive-Cauchy polymer extra stress per dumbbell: mean(q outer F_spring)-kT*I; multiply by specified number density for stress per volume";
    d["viscosity_scope"]="no number density or solvent viscosity supplied; no absolute suspension viscosity claim";
    d["conformation_convention"]="mean(q outer q) per dumbbell";
    d["fene_domain_policy"]="strict |q|<q0 before/after every step; failure is transactional including RNG; no clamp, rejection resampling or discarded sample; Euler is not domain-safe";
    d["stability_scope"]="mobility*H*dt<1 checks the Hookean eigenvalue limit; finite-step bias and nonlinear FENE stability still require convergence checks";
    d["axis_order"]="sample: frame,dumbbell,coordinate; tensor frame,row,column";
    return d;
}
py::dict frame_dict(const bdmc::RepulsiveFrame& f,const bdmc::RepulsiveSpec& s,bool explicit_state) {
    py::dict d; const std::vector<py::ssize_t> shape{size(s.n_particles),3};
    d["positions"]=copied(f.positions,shape); d["unwrapped_positions"]=copied(f.unwrapped_positions,shape);
    d["forces"]=copied(f.forces,shape); d["potential_energy"]=f.potential_energy;
    d["pair_virial"]=copied(f.pair_virial,{3,3}); d["step"]=f.step; d["time"]=f.time; d["metadata"]=metadata(s,explicit_state); return d;
}
py::dict frame_dict(const bdmc::DumbbellFrame& f,const bdmc::DumbbellSpec& s,bool explicit_state) {
    py::dict d; const std::vector<py::ssize_t> shape{size(s.n_dumbbells),3};
    d["connectors"]=copied(f.connectors,shape); d["spring_forces"]=copied(f.spring_forces,shape);
    d["conformation"]=copied(f.conformation,{3,3}); d["polymer_stress"]=copied(f.polymer_stress,{3,3});
    d["mean_energy"]=f.mean_energy; d["shear_rate"]=f.shear_rate; d["step"]=f.step; d["time"]=f.time; d["metadata"]=metadata(s,explicit_state); return d;
}
py::dict batch_dict(const std::vector<bdmc::RepulsiveFrame>& frames,const bdmc::RepulsiveSpec& s,bool explicit_state) {
    const auto count=size(frames.size()); const std::vector<py::ssize_t> shape{count,size(s.n_particles),3};
    py::array_t<double> positions(shape),unwrapped(shape),forces(shape),virial({count,py::ssize_t{3},py::ssize_t{3}}),energies(count),times(count);
    py::array_t<std::uint64_t> steps(count);
    for(std::size_t i=0;i<frames.size();++i) {
        const auto& f=frames[i]; const auto offset=i*s.n_particles*3;
        std::copy(f.positions.begin(),f.positions.end(),positions.mutable_data()+offset);
        std::copy(f.unwrapped_positions.begin(),f.unwrapped_positions.end(),unwrapped.mutable_data()+offset);
        std::copy(f.forces.begin(),f.forces.end(),forces.mutable_data()+offset);
        std::copy(f.pair_virial.begin(),f.pair_virial.end(),virial.mutable_data()+i*9);
        energies.mutable_data()[i]=f.potential_energy; times.mutable_data()[i]=f.time; steps.mutable_data()[i]=f.step;
    }
    py::dict d; d["positions"]=positions; d["unwrapped_positions"]=unwrapped; d["forces"]=forces;
    d["pair_virial"]=virial; d["potential_energy"]=energies; d["times"]=times; d["steps"]=steps;
    d["metadata"]=metadata(s,explicit_state); return d;
}
py::dict batch_dict(const std::vector<bdmc::DumbbellFrame>& frames,const bdmc::DumbbellSpec& s,bool explicit_state) {
    const auto count=size(frames.size()); const std::vector<py::ssize_t> shape{count,size(s.n_dumbbells),3};
    const std::vector<py::ssize_t> tensors{count,3,3};
    py::array_t<double> connectors(shape),forces(shape),conformation(tensors),stress(tensors),energies(count),times(count),rates(count);
    py::array_t<std::uint64_t> steps(count);
    for(std::size_t i=0;i<frames.size();++i) {
        const auto& f=frames[i]; const auto offset=i*s.n_dumbbells*3;
        std::copy(f.connectors.begin(),f.connectors.end(),connectors.mutable_data()+offset);
        std::copy(f.spring_forces.begin(),f.spring_forces.end(),forces.mutable_data()+offset);
        std::copy(f.conformation.begin(),f.conformation.end(),conformation.mutable_data()+i*9);
        std::copy(f.polymer_stress.begin(),f.polymer_stress.end(),stress.mutable_data()+i*9);
        energies.mutable_data()[i]=f.mean_energy; times.mutable_data()[i]=f.time; steps.mutable_data()[i]=f.step; rates.mutable_data()[i]=f.shear_rate;
    }
    py::dict d; d["connectors"]=connectors; d["spring_forces"]=forces; d["conformation"]=conformation;
    d["polymer_stress"]=stress; d["mean_energy"]=energies; d["times"]=times; d["steps"]=steps; d["shear_rates"]=rates;
    d["metadata"]=metadata(s,explicit_state); return d;
}
template<class Sim,class Spec> void simulation_methods(py::class_<Sim>& c) {
    c.def_property_readonly("spec",&Sim::spec)
     .def_property_readonly("metadata",[](const Sim& s){return metadata(s.spec(),s.explicit_initial_state());})
     .def("advance",&Sim::advance,py::arg("steps"),py::call_guard<py::gil_scoped_release>())
     .def("observe",[](const Sim& s){auto f=[&](){py::gil_scoped_release release;return s.observe();}(); return frame_dict(f,s.spec(),s.explicit_initial_state());})
     .def("sample",[](Sim& s,std::uint64_t steps,std::uint64_t stride,bool initial){auto frames=[&](){py::gil_scoped_release release;return s.sample(steps,stride,initial);}();return batch_dict(frames,s.spec(),s.explicit_initial_state());},py::arg("steps"),py::arg("stride")=1,py::arg("include_initial")=true);
}
} // namespace

PYBIND11_MODULE(_reference,m) {
    m.doc()="Additive repulsive-particle and connector-space kinetic dumbbell references";
    py::class_<bdmc::RepulsiveSpec> r(m,"RepulsiveSpec");
    r.def(py::init([](std::size_t n,double box,double dt,double mu,double D,double epsilon,double sigma,std::uint64_t seed,const std::string& noise,const std::string& search,const std::string& backend,int threads){
        bdmc::RepulsiveSpec s; s.n_particles=n;s.box_length=box;s.dt=dt;s.mobility=mu;s.diffusion=D;s.epsilon=epsilon;s.sigma=sigma;s.seed=seed;s.noise_policy=bdmc::parse_noise_policy(noise);s.force_search=bdmc::parse_spatial_backend(search);s.execution_backend=backend;s.backend_threads=threads;s.validate();return s;
    }),py::arg("n_particles")=64,py::arg("box_length")=10.,py::arg("dt")=0.001,py::arg("mobility")=1.,py::arg("diffusion")=1.,py::arg("epsilon")=1.,py::arg("sigma")=1.,py::arg("seed")=0,py::arg("noise_policy")="gaussian",py::arg("force_search")="cell_list",py::arg("execution_backend")="serial",py::arg("backend_threads")=1)
     .def_readwrite("n_particles",&bdmc::RepulsiveSpec::n_particles).def_readwrite("box_length",&bdmc::RepulsiveSpec::box_length)
     .def_readwrite("dt",&bdmc::RepulsiveSpec::dt).def_readwrite("mobility",&bdmc::RepulsiveSpec::mobility).def_readwrite("diffusion",&bdmc::RepulsiveSpec::diffusion)
     .def_readwrite("epsilon",&bdmc::RepulsiveSpec::epsilon).def_readwrite("sigma",&bdmc::RepulsiveSpec::sigma).def_readwrite("seed",&bdmc::RepulsiveSpec::seed)
     .def_property("noise_policy",[](const bdmc::RepulsiveSpec& s){return bdmc::noise_name(s.noise_policy);},[](bdmc::RepulsiveSpec& s,const std::string& v){s.noise_policy=bdmc::parse_noise_policy(v);})
     .def_property("force_search",[](const bdmc::RepulsiveSpec& s){return bdmc::spatial_backend_name(s.force_search);},[](bdmc::RepulsiveSpec& s,const std::string& v){s.force_search=bdmc::parse_spatial_backend(v);})
     .def_readwrite("execution_backend",&bdmc::RepulsiveSpec::execution_backend).def_readwrite("backend_threads",&bdmc::RepulsiveSpec::backend_threads)
     .def_property_readonly("dimension",[](const bdmc::RepulsiveSpec&){return 3;})
     .def("validate",&bdmc::RepulsiveSpec::validate).def("to_dict",[](const bdmc::RepulsiveSpec& s){return spec_dict(s);});
    py::class_<bdmc::RepulsiveSimulation> rs(m,"RepulsiveSimulation");
    rs.def(py::init([](bdmc::RepulsiveSpec s,py::object a){auto values=initial(a,s.n_particles,"initial_positions must have shape (n_particles,3)");py::gil_scoped_release release;return std::make_unique<bdmc::RepulsiveSimulation>(s,std::move(values));}),py::arg("spec"),py::arg("initial_positions")=py::none());
    simulation_methods<bdmc::RepulsiveSimulation,bdmc::RepulsiveSpec>(rs);

    py::class_<bdmc::DumbbellSpec> d(m,"DumbbellSpec");
    d.def(py::init([](std::size_t n,double dt,double mu,double D,double H,const std::string& law,double q0,const std::string& protocol,double rate,double strain,double amp,double omega,std::uint64_t seed,const std::string& noise,const std::string& backend,int threads){
        bdmc::DumbbellSpec s;s.n_dumbbells=n;s.dt=dt;s.mobility=mu;s.diffusion=D;s.spring_constant=H;s.spring_law=bdmc::parse_spring_law(law);s.q0=q0;s.protocol=bdmc::parse_flow_protocol(protocol);s.shear_rate=rate;s.step_strain=strain;s.strain_amplitude=amp;s.angular_frequency=omega;s.seed=seed;s.noise_policy=bdmc::parse_noise_policy(noise);s.execution_backend=backend;s.backend_threads=threads;s.validate();return s;
    }),py::arg("n_dumbbells")=64,py::arg("dt")=0.001,py::arg("mobility")=1.,py::arg("diffusion")=1.,py::arg("spring_constant")=1.,py::arg("spring_law")="hookean",py::arg("q0")=10.,py::arg("protocol")="no_flow",py::arg("shear_rate")=0.,py::arg("step_strain")=0.,py::arg("strain_amplitude")=0.,py::arg("angular_frequency")=1.,py::arg("seed")=0,py::arg("noise_policy")="gaussian",py::arg("execution_backend")="serial",py::arg("backend_threads")=1)
     .def_readwrite("n_dumbbells",&bdmc::DumbbellSpec::n_dumbbells).def_readwrite("dt",&bdmc::DumbbellSpec::dt)
     .def_readwrite("mobility",&bdmc::DumbbellSpec::mobility).def_readwrite("diffusion",&bdmc::DumbbellSpec::diffusion)
     .def_readwrite("spring_constant",&bdmc::DumbbellSpec::spring_constant).def_readwrite("q0",&bdmc::DumbbellSpec::q0)
     .def_readwrite("shear_rate",&bdmc::DumbbellSpec::shear_rate).def_readwrite("step_strain",&bdmc::DumbbellSpec::step_strain)
     .def_readwrite("strain_amplitude",&bdmc::DumbbellSpec::strain_amplitude).def_readwrite("angular_frequency",&bdmc::DumbbellSpec::angular_frequency)
     .def_readwrite("seed",&bdmc::DumbbellSpec::seed)
     .def_property("spring_law",[](const bdmc::DumbbellSpec& s){return bdmc::spring_law_name(s.spring_law);},[](bdmc::DumbbellSpec& s,const std::string& v){s.spring_law=bdmc::parse_spring_law(v);})
     .def_property("protocol",[](const bdmc::DumbbellSpec& s){return bdmc::flow_protocol_name(s.protocol);},[](bdmc::DumbbellSpec& s,const std::string& v){s.protocol=bdmc::parse_flow_protocol(v);})
     .def_property("noise_policy",[](const bdmc::DumbbellSpec& s){return bdmc::noise_name(s.noise_policy);},[](bdmc::DumbbellSpec& s,const std::string& v){s.noise_policy=bdmc::parse_noise_policy(v);})
     .def_readwrite("execution_backend",&bdmc::DumbbellSpec::execution_backend).def_readwrite("backend_threads",&bdmc::DumbbellSpec::backend_threads)
     .def_property_readonly("dimension",[](const bdmc::DumbbellSpec&){return 3;})
     .def("validate",&bdmc::DumbbellSpec::validate).def("to_dict",[](const bdmc::DumbbellSpec& s){return spec_dict(s);});
    py::class_<bdmc::DumbbellSimulation> ds(m,"DumbbellSimulation");
    ds.def(py::init([](bdmc::DumbbellSpec s,py::object a){auto values=initial(a,s.n_dumbbells,"initial_connectors must have shape (n_dumbbells,3)");py::gil_scoped_release release;return std::make_unique<bdmc::DumbbellSimulation>(s,std::move(values));}),py::arg("spec"),py::arg("initial_connectors")=py::none());
    simulation_methods<bdmc::DumbbellSimulation,bdmc::DumbbellSpec>(ds);
    m.def("execution_backends",[]{py::list result; for(const auto& c:bdmc::execution_backends()){py::dict d;d["name"]=c.name;d["precision"]=c.precision;d["compiled"]=c.compiled;d["available"]=c.available;d["reason"]=c.reason;result.append(d);} return result;});
    m.def("repulsive_pair_energy",&bdmc::repulsive_pair_energy,py::arg("r"),py::arg("epsilon"),py::arg("sigma"));
    m.def("repulsive_radial_force",&bdmc::repulsive_radial_force,py::arg("r"),py::arg("epsilon"),py::arg("sigma"));
    m.def("dumbbell_spring_energy",&bdmc::dumbbell_spring_energy,py::arg("q"),py::arg("spec"));
    m.def("dumbbell_spring_force",&bdmc::dumbbell_spring_force,py::arg("q"),py::arg("spec"));
}
