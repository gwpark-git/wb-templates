#include "bdmc/structure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace bdmc {
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
std::size_t product(std::size_t a, std::size_t b) {
    if (a && b > std::numeric_limits<std::size_t>::max()/a)
        throw std::overflow_error("analysis array shape overflows size_t");
    return a*b;
}
double represented(long double value, const char* message) {
    const auto out = static_cast<double>(value);
    if (!std::isfinite(out) || (value != 0.0L && out == 0.0)) throw std::overflow_error(message);
    return out;
}
void edges(const std::vector<double>& values, const char* message) {
    require(values.size() >= 2, message);
    for (std::size_t i=0; i<values.size(); ++i)
        require(std::isfinite(values[i]) && (i==0 || values[i] > values[i-1]), message);
}
std::size_t bin(const std::vector<double>& e, double value) {
    if (value < e.front() || value > e.back()) return e.size();
    if (value == e.back()) return e.size()-2;
    return static_cast<std::size_t>(std::upper_bound(e.begin(), e.end(), value)-e.begin()-1);
}
void increment(std::uint64_t& count, std::uint64_t delta=1) {
    if (delta > std::numeric_limits<std::uint64_t>::max()-count)
        throw std::overflow_error("analysis pair count overflows uint64");
    count += delta;
}
void shape(const std::vector<double>& positions, std::size_t f, std::size_t n, std::size_t d) {
    require(f>0 && n>=2, "pair correlation needs at least one frame and two particles");
    require(positions.size()==product(product(f,n),d), "positions shape differs from frames*particles*dimension");
    for (double x : positions) require(std::isfinite(x), "positions must be finite");
}
PairCorrelationResult initialized(std::size_t f, std::size_t n, const PairCorrelationSpec& s,
                                  const std::vector<double>& volumes) {
    PairCorrelationResult r;
    r.n_frames=f; r.n_particles=n; r.bin_volumes=volumes;
    long double volume=1;
    for(double L : s.box_lengths) volume*=L;
    r.volume=represented(volume,"box volume is outside floating-point range");
    r.number_density=represented(static_cast<long double>(n)/volume,"number density is outside floating-point range");
    // A particle has N-1 distinct possible neighbors in a finite box. Keep that
    // ideal-gas normalization distinct from N/V and preserve raw ordered counts
    // so coordination can be recovered without mixing density conventions.
    const auto density=static_cast<long double>(s.finite_population ? n-1 : n)/volume;
    r.normalization_density=represented(density,"normalization density is outside floating-point range");
    r.counts.resize(volumes.size()); r.g.resize(volumes.size()); r.expected_counts.reserve(volumes.size());
    for(double v:volumes) {
        require(std::isfinite(v) && v>0,"bin volumes must be finite and positive");
        r.expected_counts.push_back(represented(static_cast<long double>(f)*n*density*v,
                                               "expected pair count is outside floating-point range"));
    }
    return r;
}
void finish(PairCorrelationResult& r, bool cumulative=false) {
    std::uint64_t total=0;
    for(std::size_t k=0;k<r.counts.size();++k) {
        r.g[k]=represented(static_cast<long double>(r.counts[k])/r.expected_counts[k],"pair correlation overflows");
        if(cumulative) {
            increment(total,r.counts[k]);
            r.cumulative_coordination.push_back(static_cast<double>(static_cast<long double>(total)/r.n_frames/r.n_particles));
        }
    }
}
double shell(double a, double b, std::size_t d) {
    const auto lo=static_cast<long double>(a), hi=static_cast<long double>(b), delta=hi-lo;
    const auto value=d==1 ? 2*delta : d==2 ? pi*delta*(hi+lo) : (4*pi/3)*delta*(hi*hi+hi*lo+lo*lo);
    return represented(value,"shell volume is outside floating-point range");
}
std::array<double,3> image(std::array<double,3> dr,const PairCorrelationSpec& s,const OrthorhombicBox& box) {
    if(s.shear_displacement!=0.) {
        const double Ly=s.box_lengths[s.gradient_axis];
        const double crossings=std::floor(dr[s.gradient_axis]/Ly+0.5);
        dr[s.gradient_axis]-=crossings*Ly;
        dr[s.flow_axis]-=crossings*std::remainder(s.shear_displacement,s.box_lengths[s.flow_axis]);
    }
    for(std::size_t c=0;c<s.box_lengths.size();++c) dr[c]=box.minimum_image(dr[c],c);
    return dr;
}
template<class F> void pairs(const std::vector<double>& positions,std::size_t frames,std::size_t n,
                             const PairCorrelationSpec& s,bool within_radius,
                             std::uint64_t& candidates,F&& accept) {
    const auto d=s.box_lengths.size();
    const OrthorhombicBox box(s.box_lengths);
    for(std::size_t frame=0;frame<frames;++frame) {
        const auto begin=positions.begin()+static_cast<std::ptrdiff_t>(frame*n*d);
        std::vector<double> p(begin,begin+static_cast<std::ptrdiff_t>(n*d));
        if(within_radius && s.shear_displacement==0.) {
            const auto found=SpatialSearch({s.box_lengths,s.radial_edges.back(),s.backend})(p);
            increment(candidates,found.stats.candidate_pairs);
            for(const auto& q:found.pairs) accept(q.displacement,q.distance);
        } else {
            if(s.shear_displacement!=0.) {
                const long double shift=std::remainder(s.shear_displacement,s.box_lengths[s.flow_axis]);
                for(std::size_t i=0;i<n;++i) {
                    const long double crossings=std::floor(static_cast<long double>(p[i*d+s.gradient_axis])/
                                                          s.box_lengths[s.gradient_axis]);
                    const long double shifted=static_cast<long double>(p[i*d+s.flow_axis])-crossings*shift;
                    p[i*d+s.flow_axis]=static_cast<double>(std::fmod(shifted,s.box_lengths[s.flow_axis]));
                }
            }
            for(std::size_t k=0;k<p.size();++k) p[k]=box.wrap(p[k],k%d);
            for(std::size_t i=0;i<n;++i) for(std::size_t j=i+1;j<n;++j) {
                increment(candidates);
                std::array<double,3> dr{};
                for(std::size_t c=0;c<d;++c) dr[c]=p[j*d+c]-p[i*d+c];
                dr=image(dr,s,box);
                double r=0;
                for(std::size_t c=0;c<d;++c) r=std::hypot(r,dr[c]);
                if(!within_radius || r<=s.radial_edges.back()) accept(dr,r);
            }
        }
    }
}
double norm(const std::array<double,3>& a) {return std::hypot(a[0],a[1],a[2]);}
double dot(const std::array<double,3>& a,const std::array<double,3>& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
std::array<std::array<double,3>,3> basis(const DirectionalCorrelationSpec& s) {
    auto u=s.polar_axis, e=s.azimuth_reference;
    for(double v:u) require(std::isfinite(v),"polar_axis must be finite");
    for(double v:e) require(std::isfinite(v),"azimuth_reference must be finite");
    const double nu=norm(u), ne=norm(e);
    require(std::isfinite(nu) && nu>0 && std::isfinite(ne) && ne>0,"direction axes must be nonzero and finite");
    for(auto& v:u) v/=nu;
    for(auto& v:e) v/=ne;
    const auto projection=dot(e,u);
    for(std::size_t c=0;c<3;++c) e[c]-=projection*u[c];
    const auto length=norm(e);
    require(length>64*std::numeric_limits<double>::epsilon(),"azimuth_reference must not be parallel to polar_axis");
    for(auto& v:e) v/=length;
    return {{u,e,{{u[1]*e[2]-u[2]*e[1],u[2]*e[0]-u[0]*e[2],u[0]*e[1]-u[1]*e[0]}}}};
}
void validate_geometry(const PairCorrelationSpec& s) {
    (void)OrthorhombicBox(s.box_lengths);
    (void)spatial_backend_name(s.backend);
    require(std::isfinite(s.shear_displacement),"shear displacement must be finite");
    if(s.shear_displacement!=0.) {
        require(s.box_lengths.size()>=2 && s.flow_axis<s.box_lengths.size() &&
                s.gradient_axis<s.box_lengths.size() && s.flow_axis!=s.gradient_axis,"invalid shear axes");
        require(s.backend==SpatialBackend::BruteForce,"sheared analysis currently requires brute_force backend");
    }
}
} // namespace

void PairCorrelationSpec::validate() const {
    validate_geometry(*this); edges(radial_edges,"radial edges must be finite and strictly increasing");
    require(radial_edges.front()>=0,"radial edges must be nonnegative");
    require(radial_edges.back()<=0.5**std::min_element(box_lengths.begin(),box_lengths.end()),
            "radial cutoff must not exceed half the shortest periodic box");
}
PairCorrelationResult pair_correlation(const std::vector<double>& positions,std::size_t f,std::size_t n,
                                      const PairCorrelationSpec& s) {
    s.validate(); shape(positions,f,n,s.box_lengths.size());
    std::vector<double> volumes;
    for(std::size_t k=1;k<s.radial_edges.size();++k) volumes.push_back(shell(s.radial_edges[k-1],s.radial_edges[k],s.box_lengths.size()));
    auto out=initialized(f,n,s,volumes);
    pairs(positions,f,n,s,true,out.candidate_pairs,[&](const auto&,double r) {
        const auto k=bin(s.radial_edges,r);
        if(k<out.counts.size()) increment(out.counts[k],2);
    });
    finish(out,true); return out;
}

void DirectionalCorrelationSpec::validate() const {
    radial.validate(); require(radial.box_lengths.size()==3,"directional solid-angle correlation requires 3D positions");
    edges(mu_edges,"mu edges must be finite and strictly increasing");
    edges(phi_edges,"phi edges must be finite and strictly increasing");
    require(mu_edges.front()>=-1 && mu_edges.back()<=1,"mu edges must lie in [-1,1]");
    require(phi_edges.front()>=0 && phi_edges.back()<=2*pi,"phi edges must lie in [0,2*pi]");
    (void)basis(*this);
}
DirectionalCorrelationResult directional_pair_correlation(const std::vector<double>& positions,std::size_t f,
    std::size_t n,const DirectionalCorrelationSpec& s) {
    s.validate(); shape(positions,f,n,3);
    DirectionalCorrelationResult out;
    out.shape={{s.radial.radial_edges.size()-1,s.mu_edges.size()-1,s.phi_edges.size()-1}};
    const auto size=product(product(out.shape[0],out.shape[1]),out.shape[2]);
    std::vector<double> volumes; volumes.reserve(size);
    for(std::size_t r=0;r<out.shape[0];++r) for(std::size_t m=0;m<out.shape[1];++m) for(std::size_t p=0;p<out.shape[2];++p) {
        volumes.push_back(represented(static_cast<long double>(shell(s.radial.radial_edges[r],s.radial.radial_edges[r+1],3))*
            (s.mu_edges[m+1]-s.mu_edges[m])*(s.phi_edges[p+1]-s.phi_edges[p])/(4*pi),"sector volume overflows"));
    }
    out.histogram=initialized(f,n,s.radial,volumes);
    const auto axes=basis(s);
    pairs(positions,f,n,s.radial,true,out.histogram.candidate_pairs,[&](auto dr,double r) {
        const auto ir=bin(s.radial.radial_edges,r);
        if(ir>=out.shape[0]) return;
        if(r==0.) { increment(out.zero_distance_pairs,2); return; }
        for(double sign : {1.,-1.}) {
            std::array<double,3> u{};
            for(std::size_t c=0;c<3;++c) u[c]=sign*dr[c]/r;
            const double mu=std::max(-1.,std::min(1.,dot(u,axes[0])));
            double phi=std::atan2(dot(u,axes[2]),dot(u,axes[1]));
            if(phi<0.) phi+=2*pi;
            if(phi==2*pi) phi=0.;
            const auto im=bin(s.mu_edges,mu), ip=bin(s.phi_edges,phi);
            if(im<out.shape[1] && ip<out.shape[2]) increment(out.histogram.counts[(ir*out.shape[1]+im)*out.shape[2]+ip]);
        }
    });
    finish(out.histogram); return out;
}

void PlanarCorrelationSpec::validate() const {
    validate_geometry(geometry);
    require(geometry.box_lengths.size()==3,"planar projection requires 3D positions");
    require(geometry.flow_axis<3 && geometry.gradient_axis<3 && geometry.flow_axis!=geometry.gradient_axis,"invalid planar axes");
    edges(flow_edges,"flow edges must be finite and strictly increasing");
    edges(gradient_edges,"gradient edges must be finite and strictly increasing");
    const double Lx=geometry.box_lengths[geometry.flow_axis],Ly=geometry.box_lengths[geometry.gradient_axis];
    require(flow_edges.front()>=-Lx/2 && flow_edges.back()<=Lx/2 &&
            gradient_edges.front()>=-Ly/2 && gradient_edges.back()<=Ly/2,"planar edges must lie within half-box displacement ranges");
    const std::size_t normal=3-geometry.flow_axis-geometry.gradient_axis;
    require(std::isfinite(slab_half_width) && (slab_half_width==-1 || (slab_half_width>0 && slab_half_width<=geometry.box_lengths[normal]/2)),
            "slab_half_width must be -1 (full box) or positive and at most half the normal box");
}
PlanarCorrelationResult planar_pair_correlation(const std::vector<double>& positions,std::size_t f,std::size_t n,
                                               const PlanarCorrelationSpec& s) {
    s.validate(); shape(positions,f,n,3);
    const auto& g=s.geometry;
    const std::size_t normal=3-g.flow_axis-g.gradient_axis;
    PlanarCorrelationResult out;
    out.shape={{s.flow_edges.size()-1,s.gradient_edges.size()-1}};
    out.slab_width=s.slab_half_width==-1?g.box_lengths[normal]:2*s.slab_half_width;
    std::vector<double> volumes; volumes.reserve(product(out.shape[0],out.shape[1]));
    for(std::size_t i=0;i<out.shape[0];++i) for(std::size_t j=0;j<out.shape[1];++j)
        volumes.push_back(represented(static_cast<long double>(s.flow_edges[i+1]-s.flow_edges[i])*
            (s.gradient_edges[j+1]-s.gradient_edges[j])*out.slab_width,"planar bin volume overflows"));
    out.histogram=initialized(f,n,g,volumes);
    pairs(positions,f,n,g,false,out.histogram.candidate_pairs,[&](const auto& dr,double) {
        if(s.slab_half_width!=-1 && std::abs(dr[normal])>s.slab_half_width) return;
        for(double sign:{1.,-1.}) {
            const auto i=bin(s.flow_edges,sign*dr[g.flow_axis]),j=bin(s.gradient_edges,sign*dr[g.gradient_axis]);
            if(i<out.shape[0] && j<out.shape[1]) increment(out.histogram.counts[i*out.shape[1]+j]);
        }
    });
    finish(out.histogram); return out;
}

AutocorrelationResult autocorrelation(const std::vector<double>& values,std::size_t n,std::size_t c,
    std::size_t lag,bool center,CorrelationWindow window) {
    require(n>0 && c>0 && values.size()==product(n,c),"autocorrelation needs a nonempty time/component array");
    require(lag<n,"max_lag must be smaller than sample count");
    require(window==CorrelationWindow::AllOrigins || window==CorrelationWindow::LegacyFixedWindow,"unknown correlation window");
    if(window==CorrelationWindow::LegacyFixedWindow)
        require(n>=4 && lag<n/2,"legacy_fixed_window requires Nt>=4 and max_lag<floor(Nt/2)");
    std::vector<long double> mean(c,0.);
    for(std::size_t i=0;i<n;++i) for(std::size_t k=0;k<c;++k) {
        const auto value=values[i*c+k]; require(std::isfinite(value),"autocorrelation values must be finite");
        if(center) mean[k]+=static_cast<long double>(value)/n;
    }
    AutocorrelationResult out;
    out.n_components=c; out.values.resize(product(lag+1,c)); out.origin_counts.resize(lag+1);
    for(std::size_t t=0;t<=lag;++t) {
        const auto origins=window==CorrelationWindow::AllOrigins?n-t:n/2-1;
        out.origin_counts[t]=origins;
        for(std::size_t k=0;k<c;++k) {
            long double sum=0.;
            for(std::size_t i=0;i<origins;++i) sum+=(values[i*c+k]-mean[k])*(values[(i+t)*c+k]-mean[k]);
            out.values[t*c+k]=represented(sum/origins,"autocorrelation result overflows");
        }
    }
    return out;
}
} // namespace bdmc
