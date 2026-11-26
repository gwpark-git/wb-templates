#include "bdmc/structure.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <limits>
#include <optional>

namespace py=pybind11;
namespace {
py::ssize_t ssize(std::size_t n) {
    if(n>static_cast<std::size_t>(std::numeric_limits<py::ssize_t>::max()))
        throw std::overflow_error("analysis array exceeds ssize_t");
    return static_cast<py::ssize_t>(n);
}
template<class T> py::array_t<T> copied(const std::vector<T>& values,std::vector<py::ssize_t> shape={}) {
    if(shape.empty()) shape={ssize(values.size())};
    py::array_t<T> out(shape);
    std::copy(values.begin(),values.end(),out.mutable_data());
    out.attr("setflags")(false);
    return out;
}
struct Input {
    std::vector<double> values;
    std::size_t frames,particles,dimension;
};
Input positions(py::object object) {
    const auto a=py::array_t<double,py::array::c_style|py::array::forcecast>::ensure(object);
    if(!a || (a.ndim()!=2 && a.ndim()!=3))
        throw std::invalid_argument("positions must have shape (N,d) or (F,N,d)");
    Input out;
    out.dimension=static_cast<std::size_t>(a.shape(a.ndim()-1));
    out.particles=static_cast<std::size_t>(a.shape(a.ndim()-2));
    out.frames=a.ndim()==3?static_cast<std::size_t>(a.shape(0)):1;
    out.values.assign(a.data(),a.data()+a.size());
    return out;
}
bdmc::PairCorrelationSpec spec(const Input& input,const std::vector<double>& lengths,
    const std::vector<double>& radial,const std::string& normalization,const std::string& backend,
    double shear,std::size_t flow,std::size_t gradient) {
    if(input.dimension!=lengths.size()) throw std::invalid_argument("position dimension and box lengths differ");
    bdmc::PairCorrelationSpec out;
    out.box_lengths=lengths; out.radial_edges=radial;
    if(normalization!="finite_N" && normalization!="bulk")
        throw std::invalid_argument("normalization must be 'finite_N' or 'bulk'");
    out.finite_population=normalization=="finite_N";
    out.backend=bdmc::parse_spatial_backend(backend=="auto"?(shear==0.?"cell_list":"brute_force"):backend);
    out.shear_displacement=shear; out.flow_axis=flow; out.gradient_axis=gradient;
    return out;
}
py::dict output(const bdmc::PairCorrelationResult& r,const bdmc::PairCorrelationSpec& s,
                std::vector<py::ssize_t> shape,const char* kind) {
    py::dict out;
    out["schema"]="bdmc.pair_correlation.v1"; out["kind"]=kind;
    out["counts"]=copied(r.counts,shape); out["g"]=copied(r.g,shape);
    out["bin_volumes"]=copied(r.bin_volumes,shape); out["expected_counts"]=copied(r.expected_counts,shape);
    out["n_frames"]=r.n_frames; out["n_particles"]=r.n_particles;
    out["dimension"]=s.box_lengths.size(); out["volume"]=r.volume;
    out["number_density"]=r.number_density; out["normalization_density"]=r.normalization_density;
    out["normalization"]=s.finite_population?"finite_N":"bulk";
    out["backend"]=bdmc::spatial_backend_name(s.backend); out["candidate_pairs"]=r.candidate_pairs;
    out["box_lengths"]=copied(s.box_lengths); out["shear_displacement"]=s.shear_displacement;
    out["flow_axis"]=s.flow_axis; out["gradient_axis"]=s.gradient_axis;
    out["pair_convention"]="ordered distinct-particle pairs; each unordered pair contributes +/- displacement";
    out["bin_convention"]="left-closed/right-open; final upper edge included";
    out["ownership"]="independent C-contiguous readonly owned NumPy copies";
    out["boundary"]="periodic; gradient lattice shift then flow/other-axis wrapping";
    return out;
}
} // namespace

PYBIND11_MODULE(_analysis,m) {
    m.doc()="Structural observers with explicit shell/solid-angle/3D-projection normalization";
    m.def("pair_correlation",[](py::object p,const std::vector<double>& L,const std::vector<double>& edges,
        const std::string& normalization,const std::string& backend,double shear,std::size_t flow,std::size_t gradient) {
        const auto input=positions(p); const auto s=spec(input,L,edges,normalization,backend,shear,flow,gradient);
        bdmc::PairCorrelationResult r;
        {py::gil_scoped_release release; r=bdmc::pair_correlation(input.values,input.frames,input.particles,s);}
        auto out=output(r,s,{ssize(r.g.size())},"radial");
        out["radial_edges"]=copied(edges); out["cumulative_coordination"]=copied(r.cumulative_coordination);
        out["coordination_lower_radius"]=edges.front(); return out;
    },py::arg("positions"),py::arg("box_lengths"),py::arg("radial_edges"),py::kw_only(),
       py::arg("normalization")="finite_N",py::arg("backend")="auto",py::arg("shear_displacement")=0.,
       py::arg("flow_axis")=0,py::arg("gradient_axis")=1);
    m.def("directional_pair_correlation",[](py::object p,const std::vector<double>& L,const std::vector<double>& edges,
        const std::vector<double>& mu,const std::vector<double>& phi,const std::array<double,3>& axis,
        const std::array<double,3>& reference,const std::string& normalization,const std::string& backend,
        double shear,std::size_t flow,std::size_t gradient) {
        const auto input=positions(p);
        bdmc::DirectionalCorrelationSpec s;
        s.radial=spec(input,L,edges,normalization,backend,shear,flow,gradient);
        s.mu_edges=mu; s.phi_edges=phi; s.polar_axis=axis; s.azimuth_reference=reference;
        bdmc::DirectionalCorrelationResult r;
        {py::gil_scoped_release release; r=bdmc::directional_pair_correlation(input.values,input.frames,input.particles,s);}
        auto out=output(r.histogram,s.radial,{ssize(r.shape[0]),ssize(r.shape[1]),ssize(r.shape[2])},"directional_3d");
        out["radial_edges"]=copied(edges); out["mu_edges"]=copied(mu); out["phi_edges"]=copied(phi);
        out["polar_axis"]=axis; out["azimuth_reference"]=reference;
        out["zero_distance_pairs"]=r.zero_distance_pairs;
        out["angular_convention"]="mu=cos(theta) about normalized polar_axis; phi=0 along projected azimuth_reference; right-handed";
        out["pole_convention"]="azimuth at a polar-axis displacement is assigned atan2(0,0)=0; coincident distinct particles are excluded and reported";
        return out;
    },py::arg("positions"),py::arg("box_lengths"),py::arg("radial_edges"),
       py::arg("mu_edges")=std::vector<double>{-1.,1.},
       py::arg("phi_edges")=std::vector<double>{0.,6.283185307179586476925286766559},py::kw_only(),
       py::arg("polar_axis")=std::array<double,3>{{0.,0.,1.}},py::arg("azimuth_reference")=std::array<double,3>{{1.,0.,0.}},
       py::arg("normalization")="finite_N",py::arg("backend")="auto",py::arg("shear_displacement")=0.,
       py::arg("flow_axis")=0,py::arg("gradient_axis")=1);
    m.def("planar_pair_correlation",[](py::object p,const std::vector<double>& L,const std::vector<double>& x,
        const std::vector<double>& y,double slab,const std::string& normalization,double shear,
        std::size_t flow,std::size_t gradient) {
        const auto input=positions(p);
        bdmc::PlanarCorrelationSpec s;
        // Projected displacements can exceed a radial half-box cutoff, so the
        // complete brute-force route is intentional and reported.
        s.geometry=spec(input,L,{},normalization,"brute_force",shear,flow,gradient);
        s.flow_edges=x; s.gradient_edges=y; s.slab_half_width=slab;
        bdmc::PlanarCorrelationResult r;
        {py::gil_scoped_release release; r=bdmc::planar_pair_correlation(input.values,input.frames,input.particles,s);}
        auto out=output(r.histogram,s.geometry,{ssize(r.shape[0]),ssize(r.shape[1])},"planar_3d_projection");
        out["flow_edges"]=copied(x); out["gradient_edges"]=copied(y);
        out["slab_width"]=r.slab_width; out["slab_half_width"]=slab;
        out["projection"]="3D pair density integrated over the unplotted displacement coordinate; normalized by actual slab volume";
        return out;
    },py::arg("positions"),py::arg("box_lengths"),py::arg("flow_edges"),py::arg("gradient_edges"),py::kw_only(),
       py::arg("slab_half_width")=-1.,py::arg("normalization")="finite_N",py::arg("shear_displacement")=0.,
       py::arg("flow_axis")=0,py::arg("gradient_axis")=1);
    m.def("autocorrelation",[](py::object object,std::optional<py::ssize_t> max_lag,bool center,const std::string& window) {
        const auto a=py::array_t<double,py::array::c_style|py::array::forcecast>::ensure(object);
        if(!a || (a.ndim()!=1 && a.ndim()!=2) || a.shape(0)==0 || (a.ndim()==2 && a.shape(1)==0))
            throw std::invalid_argument("autocorrelation requires shape (time,) or (time,components)");
        const auto n=static_cast<std::size_t>(a.shape(0));
        const auto c=static_cast<std::size_t>(a.ndim()==1?1:a.shape(1));
        if(max_lag && *max_lag<0) throw std::invalid_argument("max_lag must be nonnegative");
        const auto lag=max_lag?static_cast<std::size_t>(*max_lag):(n/2?n/2-1:0);
        if(window!="all_origins" && window!="legacy_fixed_window")
            throw std::invalid_argument("window must be all_origins or legacy_fixed_window");
        const auto mode=window=="all_origins"?bdmc::CorrelationWindow::AllOrigins:bdmc::CorrelationWindow::LegacyFixedWindow;
        const std::vector<double> values(a.data(),a.data()+a.size());
        bdmc::AutocorrelationResult r;
        {py::gil_scoped_release release; r=bdmc::autocorrelation(values,n,c,lag,center,mode);}
        py::dict out;
        out["values"]=copied(r.values,a.ndim()==1?std::vector<py::ssize_t>{ssize(lag+1)}:
                                                   std::vector<py::ssize_t>{ssize(lag+1),ssize(c)});
        out["origin_counts"]=copied(r.origin_counts); out["window"]=window;
        out["centered"]=center; out["schema"]="bdmc.autocorrelation.v1";
        return out;
    },py::arg("values"),py::arg("max_lag")=py::none(),py::kw_only(),py::arg("center")=false,py::arg("window")="all_origins");
}
