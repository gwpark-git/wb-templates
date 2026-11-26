#include "bdmc/structure.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double pi=3.14159265358979323846;
void check(bool condition,const char* message) {if(!condition) throw std::runtime_error(message);}
void close(double a,double b,const char* message) {check(std::abs(a-b)<=1e-12*std::max(1.,std::abs(b)),message);}
template<class F> void rejects(F&& f) {
    bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"invalid analysis contract must be rejected");
}
void radial_pair_oracles() {
    bdmc::PairCorrelationSpec s;
    s.box_lengths={10,10,10};s.radial_edges={0,1.5,2};
    const auto r=bdmc::pair_correlation({0,0,0,1,0,0},1,2,s);
    check(r.counts==std::vector<std::uint64_t>({2,0}),"ordered pair shell count");
    close(r.bin_volumes[0],4*pi/3*1.5*1.5*1.5,"exact 3D shell volume");
    close(r.g[0],1000/r.bin_volumes[0],"finite-population RDF");
    close(r.cumulative_coordination[0],1,"one exact neighbor per particle");
    s.finite_population=false;
    const auto bulk=bdmc::pair_correlation({0,0,0,1,0,0},1,2,s);
    close(bulk.g[0],r.g[0]/2,"bulk finite-N factor");
    close(bulk.cumulative_coordination[0],1,"direct coordination independent of normalization");
    for(std::size_t d=1;d<=3;++d) {
        s.box_lengths=std::vector<double>(d,10.);s.radial_edges={0,1,2};
        std::vector<double> p(2*d,0.);p[d]=1.;
        const auto q=bdmc::pair_correlation(p,1,2,s);
        check(q.counts==std::vector<std::uint64_t>({0,2}),"internal edge belongs to upper bin");
        const double expected=d==1?2.:d==2?pi:4*pi/3;
        close(q.bin_volumes[0],expected,"dimension-dependent shell Jacobian");
        p[d]=2.;const auto edge=bdmc::pair_correlation(p,1,2,s);
        check(edge.counts.back()==2,"final outer edge included");
    }
}
void directional_solid_angle_oracle() {
    bdmc::DirectionalCorrelationSpec s;
    s.radial.box_lengths={10,10,10};s.radial.radial_edges={0,2};
    s.mu_edges={-1,0,1};s.phi_edges={0,pi,2*pi};
    const auto r=bdmc::directional_pair_correlation({0,0,0,1,1,1},1,2,s);
    check(r.histogram.counts==std::vector<std::uint64_t>({0,1,1,0}),"opposite solid-angle sectors");
    for(double volume:r.histogram.bin_volumes) close(volume,8*pi/3,"sector has quarter sphere volume");
    const auto radial=bdmc::pair_correlation({0,0,0,1,1,1},1,2,s.radial);
    double integrated=0.;
    for(std::size_t k=0;k<r.histogram.g.size();++k) integrated+=r.histogram.g[k]*r.histogram.bin_volumes[k];
    close(integrated,radial.g[0]*radial.bin_volumes[0],"angular integration recovers radial density");
    const auto overlap=bdmc::directional_pair_correlation({0,0,0,0,0,0},1,2,s);
    check(overlap.zero_distance_pairs==2,"undefined coincident directions explicitly counted");
    for(auto value:overlap.histogram.counts) check(value==0,"overlap omitted from angular sectors");
}
void planar_projection_is_3d() {
    bdmc::PlanarCorrelationSpec s;
    s.geometry.box_lengths={10,10,10};s.flow_edges={-2,0,2};s.gradient_edges={-1,0,1};
    const auto full=bdmc::planar_pair_correlation({0,0,0,1,.5,.2},1,2,s);
    check(full.histogram.counts==std::vector<std::uint64_t>({1,0,0,1}),"opposite flow-gradient bins");
    close(full.histogram.bin_volumes[0],20,"full 3D projection integrates normal box length");
    s.slab_half_width=.1;
    const auto omitted=bdmc::planar_pair_correlation({0,0,0,1,.5,.2},1,2,s);
    check(omitted.histogram.counts==std::vector<std::uint64_t>({0,0,0,0}),"finite slab removes distant-normal pairs");
    s.slab_half_width=.2;
    const auto slab=bdmc::planar_pair_correlation({0,0,0,1,.5,.2},1,2,s);
    check(slab.histogram.counts==full.histogram.counts,"slab boundary included");
    close(slab.histogram.bin_volumes[0],.8,"actual finite slab volume");
    close(slab.histogram.g[0],25*full.histogram.g[0],"3D slab normalization");
}
void lees_edwards_lattice_oracle() {
    bdmc::PairCorrelationSpec s;
    s.box_lengths={40,8,6};s.radial_edges={0,.3,1.,3.};
    s.shear_displacement=127.;s.backend=bdmc::SpatialBackend::BruteForce;
    const std::vector<double> base{1,0,1,8,7.8,1};
    auto translated=base;
    translated[3]+=5*127.-3*40.;translated[4]+=5*8.;translated[5]-=4*6.;
    auto reference=[&](const std::vector<double>& p) {
        double nearest=std::numeric_limits<double>::infinity();
        for(int iy=-12;iy<=12;++iy) for(int ix=-64;ix<=64;++ix) for(int iz=-8;iz<=8;++iz) {
            const double dx=p[3]-p[0]-ix*40.-iy*127.;
            const double dy=p[4]-p[1]-iy*8.;
            const double dz=p[5]-p[2]-iz*6.;
            nearest=std::min(nearest,std::hypot(dx,dy,dz));
        }
        return nearest;
    };
    close(reference(base),.2,"independent elongated-cell image enumeration");
    close(reference(translated),.2,"independent multiple lattice crossings");
    const auto a=bdmc::pair_correlation(base,1,2,s),b=bdmc::pair_correlation(translated,1,2,s);
    check(a.counts==std::vector<std::uint64_t>({2,0,0}) && a.counts==b.counts,"LE lattice invariant radial bins");
    bdmc::DirectionalCorrelationSpec direction;direction.radial=s;
    direction.mu_edges={-1,1};direction.phi_edges={0,pi,2*pi};
    const auto da=bdmc::directional_pair_correlation(base,1,2,direction);
    const auto db=bdmc::directional_pair_correlation(translated,1,2,direction);
    check(da.histogram.counts==db.histogram.counts,"LE lattice invariant directional bins");
    bdmc::PlanarCorrelationSpec planar;planar.geometry=s;
    planar.flow_edges={-1,.1,1};planar.gradient_edges={-1,0,1};
    check(bdmc::planar_pair_correlation(base,1,2,planar).histogram.counts==
          bdmc::planar_pair_correlation(translated,1,2,planar).histogram.counts,"LE lattice invariant planar bins");
}
void batch_and_backend_oracles() {
    bdmc::PairCorrelationSpec s;s.box_lengths={5,5,5};s.radial_edges={0,.5,1,2.5};
    const std::vector<double> p{0,0,0,.2,.3,.4,4.9,4.9,4.9,2,2,2};
    const auto cell=bdmc::pair_correlation(p,1,4,s);
    s.backend=bdmc::SpatialBackend::BruteForce;
    const auto brute=bdmc::pair_correlation(p,1,4,s);
    check(cell.counts==brute.counts,"cell-list and all-pairs exact count identity");
    auto batch=p;batch.insert(batch.end(),p.begin(),p.end());
    const auto twice=bdmc::pair_correlation(batch,2,4,s);
    for(std::size_t k=0;k<brute.counts.size();++k) {
        check(twice.counts[k]==2*brute.counts[k],"frame count accumulation");
        close(twice.g[k],brute.g[k],"batch RDF normalization");
        close(twice.cumulative_coordination[k],brute.cumulative_coordination[k],"batch coordination normalization");
    }
}
void temporal_estimators() {
    const auto r=bdmc::autocorrelation({1,2,3,4},4,1,3);
    const std::vector<double> expected{7.5,20./3,5.5,4};
    check(r.origin_counts==std::vector<std::uint64_t>({4,3,2,1}),"all origins count");
    for(std::size_t k=0;k<4;++k) close(r.values[k],expected[k],"independent all-origins product sum");
    const auto old=bdmc::autocorrelation({1,2,3,4},4,1,1,false,bdmc::CorrelationWindow::LegacyFixedWindow);
    check(old.values==std::vector<double>({1,2}),"explicit old fixed window");
    const auto centered=bdmc::autocorrelation({3,3,3,3},4,1,3,true);
    for(auto v:centered.values) close(v,0,"centered constant series");
    rejects([]{bdmc::autocorrelation({1,2,3},3,1,1,false,bdmc::CorrelationWindow::LegacyFixedWindow);});
}
void invalid_contracts() {
    bdmc::PairCorrelationSpec s;s.box_lengths={10,10,10};s.radial_edges={0,1};
    rejects([&]{bdmc::pair_correlation({0,0,0},1,1,s);});
    rejects([&]{bdmc::pair_correlation({0,0,0,1,0,0},0,2,s);});
    s.radial_edges={0,1,1};rejects([&]{s.validate();});
    s.radial_edges={0,6};rejects([&]{s.validate();});
    s.radial_edges={0,1};s.shear_displacement=1;rejects([&]{s.validate();});
    bdmc::DirectionalCorrelationSpec d;d.radial.box_lengths={10,10,10};d.radial.radial_edges={0,1};
    d.azimuth_reference={0,0,1};rejects([&]{d.validate();});
    rejects([]{bdmc::autocorrelation({1,std::numeric_limits<double>::quiet_NaN()},2,1,0);});
}
} // namespace
int main() {
    const std::vector<std::pair<std::string,std::function<void()>>> tests{
        {"radial_pair_oracles",radial_pair_oracles},{"directional_solid_angle_oracle",directional_solid_angle_oracle},
        {"planar_projection_is_3d",planar_projection_is_3d},{"lees_edwards_lattice_oracle",lees_edwards_lattice_oracle},
        {"batch_and_backend_oracles",batch_and_backend_oracles},{"temporal_estimators",temporal_estimators},
        {"invalid_contracts",invalid_contracts}};
    int failed=0;
    for(const auto& test:tests) try{test.second();std::cout<<"PASS "<<test.first<<'\n';}
        catch(const std::exception& e){++failed;std::cerr<<"FAIL "<<test.first<<": "<<e.what()<<'\n';}
    return failed?1:0;
}
