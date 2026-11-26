// Independent bounded oracles for the legacy simulator. No assert() is used:
// release builds with NDEBUG retain every check. Repaired issue invariants
// are ordinary passing requirements; no expected-failure allowance remains.
#include "../lib/geometry.h"
#include "../lib/handle_association.h"
#include "../lib/time_evolution.h"
#include "../src/stochastic_HEUR_flowers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Vec = std::array<double, 3>;
using Pair = std::pair<long, long>;

void require(bool result, const std::string& message) {
    if (!result) throw std::runtime_error(message);
}

void near(double actual, double expected, const std::string& message,
          double tolerance = 2e-11) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance * (1. + std::abs(expected))) {
        std::ostringstream out;
        out << std::setprecision(17) << message << ": actual=" << actual
            << ", expected=" << expected;
        throw std::runtime_error(out.str());
    }
}

// Populate the public legacy condition table in memory. These tests exercise
// geometry/topology APIs, not the independent .inp parser and CLI tests.
struct Conditions : COND {
    Conditions() {
        const std::pair<const char*, const char*> values[] = {
            {"N_dimension", "3"}, {"Np", "4"}, {"box_dimension", "10"},
            {"LBx", "10"}, {"LBy", "10"}, {"LBz", "10"},
            {"Method", "NAPLE"}, {"Integrator", "Euler"},
            {"dt/tauR", "0.001"}, {"dt/tauB", "0.001"},
            {"basic_random_seed", "17"}, {"basic_random_seed_SS", "29"},
            {"CONTINUATION_TRAJ", "FALSE"}, {"reset_initial_time", "FALSE"},
            {"cell_list", "FALSE"}, {"cutoff_connection", "1"},
            {"MECHANICAL_PERTURBATION", "FALSE"},
            {"shear_axis", "0"}, {"shear_grad_axis", "1"},
            {"tracking_individual_chain", "TRUE"},
            {"CONTINUATION_CHAIN", "FALSE"}, {"N_chains_per_particle", "4"},
            {"N_THREADS_BD","1"}, {"N_THREADS_SS","1"}, {"Nt","5"},
            {"N_skip_ener","1"}, {"N_skip_file","1"}, {"N_steps_block","1"},
            {"Rt","1"}, {"MC_LOG","FALSE"}, {"MC_ASSOCIATION_MAP","FALSE"},
            {"MC_renewal","FALSE"}, {"record_RDIST","FALSE"},
            {"connector","Gaussian"}, {"association_probability","DEFAULT"},
            {"transition_probability","UNIFORM"}, {"chain_selection","UNIFORM"},
            {"repulsion_type","SOFT_REPULSION_P2"}, {"repulsion_coupling","GEOMETRIC_MEAN"},
            {"repulsion_coefficient","4"}, {"repulsion_coefficient_base",".04"},
            {"effective_distance","1"}, {"scale_factor_chain","1"}, {"l_cap",".12"},
            {"friction_junction","DEFAULT"}, {"probability_area","FALSE"},
            {"Step","NO_EQUILIBRATION"}
        };
        N_arg = sizeof(values) / sizeof(values[0]);
        arg = new std::string*[N_arg];
        for (long i = 0; i < N_arg; ++i) {
            arg[i] = new std::string[2];
            arg[i][0] = values[i].first;
            arg[i][1] = values[i].second;
        }
        ERR = "ERR";
    }
    void set(const std::string& key, const std::string& value) {
        (*this)(key) = value;
    }
};

double wrapped(double x, double length) {
    return x - length * std::floor(x / length);
}

Vec position(TRAJECTORY& t, long particle) {
    return {{t(0, particle, 0), t(0, particle, 1), t(0, particle, 2)}};
}

void position(TRAJECTORY& t, long particle, const Vec& p) {
    for (long k = 0; k < 3; ++k) t(0, particle, k) = p[k];
}

struct Image { Vec vector; double distance; };

// Independent integer-lattice enumeration, wider than the legacy {-1,0,1}
// search, initialized to infinity rather than a box-length sentinel.
Image lattice(const Vec& from, const Vec& to, const Vec& box, double shift) {
    Image best = {{{0., 0., 0.}}, std::numeric_limits<double>::infinity()};
    for (int nx = -4; nx <= 4; ++nx)
        for (int ny = -4; ny <= 4; ++ny)
            for (int nz = -4; nz <= 4; ++nz) {
                Vec v = {{to[0] - from[0] + nx * box[0] + ny * shift,
                          to[1] - from[1] + ny * box[1],
                          to[2] - from[2] + nz * box[2]}};
                const double distance = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
                if (distance < best.distance) best = {v, distance};
            }
    return best;
}

void compare_vector(MATRIX& actual, const Vec& expected, const std::string& what) {
    for (long k = 0; k < 3; ++k) near(actual(k), expected[k], what);
}

void periodic_geometry() {
    Conditions c;
    TRAJECTORY t(c, 2);
    const Vec box = {{10., 13., 16.}};
    for (long k = 0; k < 3; ++k) t.box_dimension[k] = box[k];
    const Vec raw[] = {{{-20.25, 39., 80.2}}, {{30., -26.5, -48.}},
                       {{9.8, 12.7, 15.6}}, {{0.3, 0.4, 0.7}}};
    for (long i = 0; i < 4; ++i) position(t, i, raw[i]);
    GEOMETRY::minimum_image_convention_loop(t, 0);
    for (long i = 0; i < 4; ++i)
        for (long k = 0; k < 3; ++k)
            near(t(0, i, k), wrapped(raw[i][k], box[k]), "periodic many-box wrap");
    MATRIX v(3, 1), reverse(3, 1);
    for (long i = 0; i < 4; ++i)
        for (long j = i + 1; j < 4; ++j) {
            const Image expected = lattice(position(t,i), position(t,j), box, 0.);
            near(GEOMETRY::get_minimum_distance(t,0,i,j,v), expected.distance,
                 "ordinary periodic lattice norm");
            compare_vector(v, expected.vector, "ordinary target-minus-subject vector");
            GEOMETRY::get_minimum_distance(t,0,j,i,reverse);
            for (long k=0; k<3; ++k) near(v(k), -reverse(k), "pair antisymmetry");
        }
}

void shear_geometry() {
    Conditions c;
    c.set("MECHANICAL_PERTURBATION", "SIMPLE_SHEAR");
    TRAJECTORY t(c, 2);
    RDIST r(c);
    const Vec box = {{10., 10., 10.}};
    MATRIX v(3,1), reverse(3,1);
    // Non-tied, deterministic interior/boundary samples; no stochastic oracle.
    for (double shift : {-5., -4.7, -2.3, 0., 1.9, 4.8, 5.}) {
        double reported_shift=0.;
        GEOMETRY::set_box_shift_factor(reported_shift,r,shift,box[0]);
        for (int n=0; n<24; ++n) {
            const Vec a = {{std::fmod(.17+1.37*n,10.), std::fmod(.09+2.19*n,10.),
                            std::fmod(.31+3.07*n,10.)}};
            const Vec b = {{std::fmod(9.81+2.43*n,10.), std::fmod(9.73+1.23*n,10.),
                            std::fmod(9.49+.89*n,10.)}};
            position(t,0,a); position(t,1,b);
            const Image expected=lattice(a,b,box,shift);
            near(r.measure_minimum_distance(r,t,0,0,1,nullptr), expected.distance,
                 "selected shear distance policy vs independent lattice");
            compare_vector(r.Rvec[0][1],expected.vector,"selected shear vector policy");
            near(GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,v,shift),
                 expected.distance,"shear lattice norm");
            compare_vector(v,expected.vector,"shear lattice image vector");
            GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,1,0,reverse,shift);
            for(long k=0;k<3;++k) near(v(k),-reverse(k),"shear antisymmetry");
        }
    }
}

void step_shear_and_shift() {
    Conditions c;
    TRAJECTORY t(c,2);
    RDIST r(c);
    const Vec original[] = {{{9.7,.4,1.}}, {{.2,9.6,2.}},
                            {{5.1,4.2,3.}}, {{1.3,8.8,4.}}};
    for(double gamma : {-2.7,-.3,.4,2.8}) {
        for(long i=0;i<4;++i) position(t,i,original[i]);
        GEOMETRY::apply_step_shear(t,0,0,1,gamma,10.);
        for(long i=0;i<4;++i) {
            near(t(0,i,0),wrapped(original[i][0]+gamma*original[i][1],10.),
                 "step shear affine x");
            near(t(0,i,1),original[i][1],"step shear unchanged gradient");
            near(t(0,i,2),original[i][2],"step shear unchanged vorticity coordinate");
        }
    }
    for(double displacement : {-32.7,-15.1,-5.,0.,5.,15.1,32.7}) {
        double shift=0.,plain=0.;
        GEOMETRY::set_box_shift_factor(shift,r,displacement,10.);
        GEOMETRY::set_box_shift_factor(plain,displacement,10.);
        // Shift is equivalent modulo L; central representation lies in [-L/2,L/2].
        near(wrapped(shift,10.),wrapped(displacement,10.),"box shift modulo L");
        near(plain,shift,"box shift overloads agree");
        require(std::abs(r.map_to_central_box_image)<=5.+1e-12,"central shift domain");
        near(wrapped(r.map_to_central_box_image,10.),wrapped(displacement,10.),
             "central image lattice equivalence");
    }
}

void single_gradient_crossings() {
    Conditions c;
    TRAJECTORY t(c,2);
    const Vec raw[] = {{{1.2,10.3,.4}}, {{8.7,-.2,1.}},
                       {{4.4,10.,2.}}, {{5.1,0.,3.}}};
    for(double shift : {-4.3,0.,3.7}) {
        for(long i=0;i<4;++i) position(t,i,raw[i]);
        GEOMETRY::apply_shear_boundary_condition(t,0,0,1,shift);
        GEOMETRY::minimum_image_convention_loop(t,0);
        for(long i=0;i<4;++i) {
            const double winding=std::floor(raw[i][1]/10.);
            near(t(0,i,0),wrapped(raw[i][0]-winding*shift,10.),"single-crossing flow shift");
            near(t(0,i,1),wrapped(raw[i][1],10.),"single-crossing gradient wrap");
            near(t(0,i,2),raw[i][2],"single-crossing z");
        }
    }
}

void connector_laws() {
    POTENTIAL_VARIABLE p;
    p.constructor_POTENTIAL_VARIABLE();
    p.N_dimension=3; p.scale_factor_chain=1.4; p.ratio_RM_R0=2.8;
    using Law=double (*)(double,POTENTIAL_VARIABLE&);
    const Law forces[]={FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_force};
    const Law energies[]={FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_potential,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_potential};
    const Law probabilities[]={FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_Boltzmann,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_Boltzmann,
        FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_Boltzmann};
    for(int law=0;law<3;++law)
        for(double x : {.13,.6,1.2,2.1}) {
            const double h=1e-5;
            const double derivative=(energies[law](x+h,p)-energies[law](x-h,p))/(2*h);
            near(forces[law](x,p),derivative,"connector dU/dr = positive tension",2e-8);
            near(probabilities[law](x,p),std::exp(-energies[law](x,p)),"Boltzmann exp(-U)");
            const double stiffness=3.*(law==0 ? 1. : 1.4*1.4);
            const double expected=stiffness*x/(law==2 ? 1.-x*x/(2.8*2.8) : 1.);
            near(forces[law](x,p),expected,"independent connector stiffness/domain law");
        }
    // Valid-domain checks stop before Rmax; out-of-domain behavior is not a
    // singular FENE law and cannot be validated by the derivative oracle above.
    near(FORCE::NAPLE::MC_ASSOCIATION::MAP_minimum_R0_modified_Gaussian_Boltzmann(.4,p),
         std::exp(-.5*3.*1.4*1.4),"supported modified Gaussian minimum-R0 clamp");
}

void kinetic_probabilities() {
    POTENTIAL_VARIABLE p;
    p.constructor_POTENTIAL_VARIABLE();
    p.N_dimension=3; p.scale_factor_chain=1.2; p.l_cap=.4; p.energy_barrier=2.;
    auto force=FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force;
    for(double dt : {0.,.001,.1,2.})
        for(double x : {.2,1.,1.7}) {
            p.delta_t0=dt;
            const double tension=3.*1.2*1.2*x;
            const double rate=std::exp(tension*.4);
            near(KINETICS::dissociation_probability(x,force,p),std::min(1.,dt*rate),
                 "linear-rate dissociation with saturation");
            near(KINETICS::FIRST_ORDER::dissociation_probability(x,force,p),
                 -std::expm1(-dt*rate),"first-order integrated hazard");
            near(KINETICS::METROPOLIS::transition_probability(x,force,p),
                 std::min(1.,std::exp(tension*.4-2.)),"Metropolis barrier");
            near(KINETICS::UNIFORM::transition_probability(x,force,p),1.,"uniform acceptance");
        }
}

std::map<Pair,long> graph(ASSOCIATION& a) {
    std::map<Pair,long> result;
    long ends=0;
    for(long i=0;i<a.Np;++i) {
        require(a.TOKEN[i]>=1 && a.HASH[i](0)==i,"permanent self entry");
        std::set<long> targets;
        long particle_ends=0;
        for(long j=0;j<a.TOKEN[i];++j) {
            const double raw_id=a.HASH[i](j), raw_weight=a.weight[i](j);
            const long id=static_cast<long>(raw_id), w=static_cast<long>(raw_weight);
            require(raw_id==id && id>=0 && id<a.Np,"valid integer adjacency ID");
            require(targets.insert(id).second,"unique adjacency targets");
            require(raw_weight==w && w>=0 && (j==0 || w>0),"integer positive bridge weight");
            if(j==0) require(w%2==0,"even self end count");
            else near(a.return_multiplicity(id,i),w,"symmetric multiplicity");
            particle_ends+=w;
            if(i<=id && w>0) result[{i,id}]=(i==id ? w/2 : w);
        }
        near(a.N_ASSOCIATION_PARTICLE(i),particle_ends,"maintained per-particle end count");
        near(a.N_CONNECTED_ENDS(i),particle_ends,"reported per-particle ends");
        ends+=particle_ends;
    }
    require(ends==2*a.Nc*a.Np,"total physical chain-end conservation");
    near(a.N_TOTAL_CONNECTED_ENDS(),ends,"reported total chain ends");
    long bridges=0;
    for(const auto& entry:result) if(entry.first.first!=entry.first.second) bridges+=entry.second;
    near(a.N_TOTAL_ASSOCIATION(),bridges,"bridge multiplicity total");
    return result;
}

void compare_tracking(ASSOCIATION& a,CHAIN_HANDLE& chains) {
    const auto actual=graph(a);
    std::map<Pair,long> reconstructed;
    std::vector<long> ends(a.Np,0), seen(2*chains.N_chains,0);
    for(long k=0;k<chains.N_chains;++k) {
        const long i=chains.HEAD(k),j=chains.TAIL(k);
        require(i>=0 && i<a.Np && j>=0 && j<a.Np,"labelled chain endpoint IDs");
        ++reconstructed[{std::min(i,j),std::max(i,j)}];
        ++ends[i]; ++ends[j];
    }
    require(actual==reconstructed,"independently reconstructed chain graph");
    for(long i=0;i<a.Np;++i) {
        require(chains.P_TOKEN[i]==ends[i],"reverse endpoint-token count");
        for(long h=0;h<chains.P_TOKEN[i];++h) {
            const long endpoint=chains.PARTICLE[i][h];
            require(endpoint>=0 && endpoint<2*chains.N_chains,"reverse endpoint ID");
            require(chains.CE_ATTACHED_REF(endpoint)==i,"reverse endpoint ownership");
            ++seen[endpoint];
        }
    }
    for(long count:seen) require(count==1,"each labelled endpoint occurs exactly once");
}

void topology_actions() {
    Conditions c;
    ASSOCIATION a(4,4,6,true);
    CHAIN_HANDLE chains(c,a);
    POTENTIAL_SET p;
    p.w_function=KINETICS::UNIFORM::detachment_weight;
    p.f_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force;
    p.force_variables.N_dimension=3;
    MATRIX distances[4];
    for(auto& row:distances) row.initial(4,1,.7);
    compare_tracking(a,chains);
    // Each row is (anchor, old endpoint, new endpoint, expected action).
    // Repeated bridges test multiplicity; deletion tests hash compaction.
    const long transitions[][4]={{0,0,1,INDEX_MC::ADD},{0,0,1,INDEX_MC::ADD},
        {0,0,2,INDEX_MC::ADD},{0,1,3,INDEX_MC::MOV},
        {0,2,2,INDEX_MC::CANCEL},{0,1,0,INDEX_MC::OPP_DEL},
        {0,2,0,INDEX_MC::OPP_DEL},{0,3,0,INDEX_MC::OPP_DEL}};
    for(const auto& row:transitions) {
        INDEX_MC idx;
        idx.beads[a.flag_itself]=row[0]; idx.beads[a.flag_other]=row[1];
        idx.beads[a.flag_new]=row[2];
        idx.beads[a.flag_hash_other]=(row[0]==row[1] ? 0 : a.FIND_HASH_INDEX(row[0],row[1]));
        const long action=ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(a,idx);
        require(action==row[3],"serial action classification");
        auto expected=graph(a);
        if(action!=INDEX_MC::CANCEL) {
            const Pair old={std::min(row[0],row[1]),std::max(row[0],row[1])};
            const Pair next={std::min(row[0],row[2]),std::max(row[0],row[2])};
            require(expected[old]>0,"oracle transition has an old chain");
            if(--expected[old]==0) expected.erase(old);
            ++expected[next];
        }
        // Match production order: aggregate action, then endpoint tracking.
        ACTION::ACT(p,a,idx,distances,action);
        chains.TRACKING_ACTION(a,action,idx);
        require(graph(a)==expected,"independent one-chain transition oracle");
        compare_tracking(a,chains);
    }
    require(a.N_TOTAL_ASSOCIATION()==0,"add/move/delete round trip returns to loops");

    // Endpoint-count constraint and simple-pair restriction are separate rules.
    ASSOCIATION fixed(4,2,0,true);
    INDEX_MC idx;
    idx.beads[fixed.flag_itself]=0; idx.beads[fixed.flag_other]=0;
    idx.beads[fixed.flag_new]=1;
    require(ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(fixed,idx)==INDEX_MC::CANCEL,
            "fixed end-count rejects an add");
    ASSOCIATION single(4,2,3,false);
    single.add_association(0,1);
    require(ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(single,idx)==INDEX_MC::CANCEL,
            "single-bridge policy rejects duplicate pair");
    graph(single);
}

void suggestion_distribution() {
    ASSOCIATION a(4,4,6,true);
    a.add_association(0,1); a.add_association(0,1); a.add_association(0,2);
    const double cases[]={.7,1.3,2.1};
    double sum=0.;
    for(long j=0;j<a.TOKEN[0];++j) {
        a.CASE_SUGGESTION[0](j)=cases[j];
        sum+=a.weight[0](j)*cases[j];
    }
    near(a.update_Z_SUGGESTION_particle(0),sum,"suggestion partition sum");
    a.update_dPDF_SUGGESTION_particle(0);
    a.update_dCDF_SUGGESTION_particle(0);
    double lower=0.;
    for(long j=0;j<a.TOKEN[0];++j) {
        const double probability=a.weight[0](j)*cases[j]/sum;
        near(a.dPDF_SUGGESTION[0](j),probability,"multiplicity weighted discrete PDF");
        near(a.dCDF_SUGGESTION[0](j),lower+probability,"suggestion cumulative PDF");
        require(a.GET_INDEX_HASH_FROM_ROLL(0,lower+probability/2.)==j,"CDF interior selection");
        if(j+1<a.TOKEN[0])
            require(a.GET_INDEX_HASH_FROM_ROLL(0,a.dCDF_SUGGESTION[0](j))==j+1,
                    "strict upper-bound CDF selection convention");
        lower+=probability;
    }
    near(lower,1.,"suggestion normalisation");
}

void network_force_and_virial() {
    Conditions c; c.set("Np","2");
    c.set("MECHANICAL_PERTURBATION","SIMPLE_SHEAR");
    TRAJECTORY t(c,2); RDIST r(c);
    position(t,0,{{.2,.3,.5}}); position(t,1,{{9.3,9.6,.9}});
    const Vec box={{10.,10.,10.}};
    const double shift=2.7;
    const Image expected=lattice(position(t,0),position(t,1),box,shift);
    double shift_out=0.; GEOMETRY::set_box_shift_factor(shift_out,r,shift,10.);
    for(long i=0;i<2;++i)
        for(long j=0;j<2;++j)
            r.Rsca[i](j)=r.measure_minimum_distance(r,t,0,i,j,nullptr);
    ASSOCIATION a(2,3,4,true);
    a.add_association(0,1); a.add_association(0,1);
    POTENTIAL_SET p;
    p.force_variables.N_dimension=3; p.force_variables.scale_factor_chain=1.3;
    p.f_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force;
    p.e_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential;
    const double stiffness=2.*3.*1.3*1.3;
    MATRIX force(3,1);
    for(long i=0;i<2;++i) {
        double xx=0.,yy=0.,zz=0.,xy=0.,xz=0.,yz=0.,energy=0.;
        INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost_with_RF(
            p,a,force,i,r.Rvec,r.Rsca,xx,yy,zz,xy,xz,yz,energy);
        for(long k=0;k<3;++k)
            near(force(k),(i==0 ? 1.:-1.)*stiffness*expected.vector[k],
                 "network force direction and bridge multiplicity");
        near(energy,.5*stiffness*expected.distance*expected.distance,
             "directed network connector energy");
        const double tensors[]={xx,yy,zz,xy,xz,yz};
        const long axes[][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
        for(int k=0;k<6;++k)
            near(tensors[k],stiffness*expected.vector[axes[k][0]]*expected.vector[axes[k][1]],
                 "directed force-phase connector tensor");
    }
}

void known_multibox_shear() {
    Conditions c; TRAJECTORY t(c,2);
    position(t,0,{{1.2,20.3,.4}});
    GEOMETRY::apply_shear_boundary_condition(t,0,0,1,2.3);
    GEOMETRY::minimum_image_convention_loop(t,0);
    near(t(0,0,0),wrapped(1.2-2.*2.3,10.),
                     "F09 two gradient crossings require two flow shifts");
}

void known_elongated_shear() {
    Conditions c; TRAJECTORY t(c,2);
    const Vec box={{10.,40.,40.}};
    for(long k=0;k<3;++k) t.box_dimension[k]=box[k];
    position(t,0,{{0.,0.,0.}}); position(t,1,{{0.,19.,19.}});
    MATRIX v(3,1);
    const double actual=GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,v,2.);
    near(actual,lattice(position(t,0),position(t,1),box,2.).distance,
         "F09 elongated shear minimum must not be capped at Lx");
}

void known_minimum_R0_gaussian() {
    POTENTIAL_VARIABLE p; p.constructor_POTENTIAL_VARIABLE(); p.N_dimension=3;
    near(FORCE::NAPLE::MC_ASSOCIATION::MAP_minimum_R0_Gaussian_Boltzmann(.4,p),
         std::exp(-1.5),"F04 Gaussian proposal minimum-R0 clamp");
}

void known_minimum_R0_fene() {
    POTENTIAL_VARIABLE p; p.constructor_POTENTIAL_VARIABLE();
    p.N_dimension=3; p.scale_factor_chain=1.2; p.ratio_RM_R0=2.8;
    const double expected=std::pow(1.-1./(2.8*2.8),1.5*1.2*1.2*2.8*2.8);
    near(FORCE::NAPLE::MC_ASSOCIATION::MAP_minimum_R0_FENE_Boltzmann(.4,p),
         expected,"F04 FENE proposal minimum-R0 clamp");
}

void known_stale_bond_geometry() {
    Conditions c; c.set("Np","2"); c.set("cell_list","TRUE");
    TRAJECTORY t(c,2); RDIST r(c);
    position(t,0,{{.2,.2,.2}}); position(t,1,{{.7,.2,.2}});
    long scratch[3];
    r.allocate_cells_from_positions(t,0,scratch);
    GEOMETRY::compute_RDIST_particle(r,0,t,0);
    near(r.Rsca[0](1),.5,"F03 initial covered bond distance");
    ASSOCIATION a(2,2,3,true); a.add_association(0,1);
    position(t,1,{{4.7,.2,.2}});
    r.allocate_cells_from_positions(t,0,scratch);
    GEOMETRY::compute_RDIST_particle(r,0,t,0);
    GEOMETRY::compute_bonded_RDIST_particle(r,a,0,t,0);
    POTENTIAL_SET p; p.force_variables.N_dimension=3;
    p.f_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force;
    MATRIX force(3,1);
    INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(p,a,force,0,r.Rvec,r.Rsca);
    near(force(0),3.*4.5,
                     "F03 bonded force must use current non-neighbor geometry");
}


void coupled_energy_and_force() {
    Conditions c; c.set("Np","2");
    TRAJECTORY t(c,2); RDIST r(c);
    position(t,0,{{.1,.2,.3}}); position(t,1,{{.5,.2,.3}});
    long scratch[3]; r.allocate_cells_from_positions(t,0,scratch);
    GEOMETRY::compute_RDIST_particle(r,0,t,0);
    MATRIX populations(2,1); populations(0)=2.; populations(1)=6.;
    POTENTIAL_SET p;
    p.force_variables.N_dimension=3; p.force_variables.effective_distance=1.;
    p.force_variables.repulsion_coefficient_base=.2;
    p.force_variables.N_ASSOCIATION_PARTICLE=&populations;
    p.f_repulsion=FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force;
    p.e_repulsion=FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential;
    p.f_repulsion_coupling=FORCE::NAPLE::SIMPLE_REPULSION::geometrical_mean_repulsive_coefficient;
    MATRIX force(3,1), plain(3,1);
    double xx=0.,yy=0.,zz=0.,xy=0.,xz=0.,yz=0.,energy=0.;
    INTEGRATOR::EULER::cal_repulsion_force_R_boost_with_RF(p,force,0,r,xx,yy,zz,xy,xz,yz,energy);
    INTEGRATOR::EULER::cal_repulsion_force_R_boost(p,plain,0,r);
    const double coupling=.2*2.*6., radius=.4, h=1e-6;
    const auto potential=[&](double x) { return coupling*(1.-x)*(1.-x)*(2.+x)/3.; };
    near(energy,potential(radius),"F05 directed coupled potential energy");
    near(force(0),(potential(radius+h)-potential(radius-h))/(2*h),"F05 force is dU/dr",1e-9);
    near(force(0),plain(0),"RF and plain coupling conventions agree");
    near(xx,force(0)*radius,"coupled repulsive virial");
    GEOMETRY::compute_RDIST_particle(r,1,t,0);
    near(ANALYSIS::cal_potential_energy_R_boost(p,r),energy,
         "standalone energy matches coupled unique-pair physical total");
}

void fractional_cutoff_and_move_weight() {
    Conditions c; ASSOCIATION a(4,4,6,true); RDIST r(c);
    a.add_association(0,1); r.Rsca[0](1)=1.5; r.Rsca[1](0)=1.5;
    long i=-1,h=-1,j=-1;
    require(!HEUR::get_index_bridge_distance_larger_than_cutoff_range(i,h,j,a,r,1.8),
            "F06 fractional cutoff must not truncate");
    require(HEUR::get_index_bridge_distance_larger_than_cutoff_range(i,h,j,a,r,1.4),
            "F06 cutoff selects truly long bridge");
    POTENTIAL_SET p; p.force_variables.N_dimension=3;
    p.f_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force;
    p.w_function=+[](double radius,double(*)(double,POTENTIAL_VARIABLE&),POTENTIAL_VARIABLE&) { return 1.+radius; };
    r.Rsca[0](2)=.7; r.Rsca[1](2)=2.3;
    INDEX_MC idx; idx.beads[a.flag_itself]=0;idx.beads[a.flag_other]=1;
    idx.beads[a.flag_new]=2;idx.beads[a.flag_hash_other]=a.FIND_HASH_INDEX(0,1);
    ACTION::ACT(p,a,idx,r.Rsca,INDEX_MC::MOV);
    near(a.CASE_SUGGESTION[0](a.FIND_HASH_INDEX(0,2)),1.7,"move uses new anchor-to-target distance");
}

void renewal_and_saturated_tracking() {
    Conditions c; ASSOCIATION a(4,4,6,true); CHAIN_HANDLE chains(c,a);
    a.add_association(0,1); a.add_association(0,1); a.add_association(2,3);
    chains.allocate_existing_bridges(a); compare_tracking(a,chains);
    a.set_initial_condition(); chains.allocate_existing_bridges(a); compare_tracking(a,chains);
    require(a.N_TOTAL_ASSOCIATION()==0,"F07 renewal returns both representations to loops");
    // All chains at one particle exercise full endpoint hashes and 2*N scratch.
    CHAIN_HANDLE concentrated(8,1);
    for(int n=0;n<100;++n)
        require(concentrated.get_index_degeneracy(0,0)>=0,"loop degeneracy scratch capacity");
    concentrated.mov_attachment(0,0);
    require(concentrated.P_TOKEN[0]==16,"full reverse hash removal and insertion conserve ends");
    std::set<long> seen;
    for(long h=0; h<16; ++h) seen.insert(concentrated.PARTICLE[0][h]);
    require(seen.size()==16,"full reverse hash has unique endpoints");
}

double zero_random(MATRIX& random, POTENTIAL_VARIABLE&) { random.set_value(0.); return 0.; }
double degree_friction(long degree,POTENTIAL_VARIABLE&) { return 1.+degree; }

void heterogeneous_particle_friction() {
    Conditions c; c.set("Method","NAPLE_ASSOCIATION");
    TRAJECTORY t(c,2); RDIST r(c); ASSOCIATION a(4,4,6,true);
    a.add_association(0,1);a.add_association(0,2);
    const Vec positions[]={{{.1,.2,.3}},{{.5,.6,.7}},{{.9,1.,1.1}},{{2.,2.,2.}}};
    for(long i=0;i<4;++i) position(t,i,positions[i]);
    long scratch[3];r.allocate_cells_from_positions(t,0,scratch);
    for(long i=0;i<4;++i) GEOMETRY::compute_RDIST_particle(r,i,t,0);
    POTENTIAL_SET p; p.force_variables.N_dimension=3;p.force_variables.repulsion_coefficient=1.;
    p.f_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force;
    p.e_connector=FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_potential;
    p.f_repulsion=FORCE::DEFAULT::EMPTY_force_contribution;
    p.e_repulsion=FORCE::DEFAULT::EMPTY_force_contribution;
    p.f_repulsion_coupling=FORCE::NAPLE::SIMPLE_REPULSION::pre_averaged_repulsive_coefficient;
    p.zeta_particle=degree_friction;p.scale_random=zero_random;
    HEUR::TEMPORAL_VARIABLE_HEUR v(c,2); RNG_BOOST rng(c);
    HEUR::OMP_time_evolution_Euler(t,0,1,a,p,v.force_repulsion,v.force_random,v.force_spring,
                                  r,v.vec_boost_Nd_parallel,rng,1,c,v);
    for(long i=0;i<4;++i) {
        Vec expected={{0.,0.,0.}};long degree=0;
        for(long h=1;h<a.TOKEN[i];++h) {
            const long j=a.HASH[i](h);const long w=a.weight[i](h);degree+=w;
            for(long k=0;k<3;++k) expected[k]+=3.*w*(positions[j][k]-positions[i][k]);
        }
        for(long k=0;k<3;++k)
            near(t(1,i,k),positions[i][k]+t.dt*expected[k]/(1.+degree),
                 "F02 all coordinates use particle degree friction");
    }
}

void serialized_topology_and_renewal() {
    Conditions one,many; one.set("Method","NAPLE_ASSOCIATION"); many.set("Method","NAPLE_ASSOCIATION");
    many.set("N_THREADS_SS","4");
    TRAJECTORY t(one,2); RDIST r(one);
    for(long i=0;i<4;++i) position(t,i,{{.1+.2*i,.2,.3}});
    long scratch[3];r.allocate_cells_from_positions(t,0,scratch);
    for(long i=0;i<4;++i) GEOMETRY::compute_RDIST_particle(r,i,t,0);
    ASSOCIATION a(4,4,6,true),b(4,4,6,true);
    CHAIN_HANDLE ca(one,a), cb(many,b);
    POTENTIAL_SET pa,pb;
    require(FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(pa,one)==0,"one-thread law setup");
    require(FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(pb,many)==0,"many-thread law setup");
    pa.force_variables.N_ASSOCIATION_PARTICLE=&a.N_ASSOCIATION_PARTICLE;
    pb.force_variables.N_ASSOCIATION_PARTICLE=&b.N_ASSOCIATION_PARTICLE;
    RNG_BOOST ra(one),rb(many);HEUR::TEMPORAL_VARIABLE_HEUR va(one,2),vb(many,2);
    RECORD_DATA da(RECORD_DATA::NoOutput{}),db(RECORD_DATA::NoOutput{});
    require(!da.traj.is_open() && !da.ener.is_open() && !da.MC_LOG.is_open(),
            "explicit no-output topology fixture opens no files");
    INDEX_MC ia[1],ib[4];LOCK la(4),lb(4);
    for(long step=0;step<80;++step) {
        HEUR::OMP_SS_topological_time_evolution(step,a,ca,pa,r,ra,ia,da,one,la,va);
        HEUR::OMP_SS_topological_time_evolution(step,b,cb,pb,r,rb,ib,db,many,lb,vb);
        compare_tracking(a,ca);compare_tracking(b,cb);
        require(graph(a)==graph(b),"F08 topology is invariant to requested SS thread count");
        for(long k=0;k<ca.N_chains;++k)
            require(ca.HEAD(k)==cb.HEAD(k) && ca.TAIL(k)==cb.TAIL(k),"F08 labelled RNG schedule agrees");
        for(long k=0;k<6;++k) require(va.cnt_arr[k]==vb.cnt_arr[k],"F08 action counts agree");
    }
    require(va.cnt_arr[INDEX_MC::ADD]>0 && va.cnt_arr[INDEX_MC::MOV]>0,
            "F08 test performs nontrivial graph transitions");
    va.MC_renewal=true;
    HEUR::OMP_SS_update_STATISTICS(a,ca,pa,r,va,da);
    compare_tracking(a,ca);require(a.N_TOTAL_ASSOCIATION()==0,"F07 production renewal path");
}

void shear_domain_and_raw_distance() {
    for(const auto& item:std::vector<std::pair<std::string,std::string>>{
            {"cell_list","TRUE"},{"shear_axis","2"},{"shear_grad_axis","0"},{"N_dimension","2"}}) {
        Conditions c;c.set("MECHANICAL_PERTURBATION","SIMPLE_SHEAR");c.set(item.first,item.second);
        bool rejected=false;try { RDIST r(c); }catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"F09 unsupported shear configuration must fail explicitly");
    }
    Conditions c;TRAJECTORY t(c,2);position(t,0,{{.1,.2,.3}});position(t,1,{{3.1,4.2,.3}});
    near(GEOMETRY::get_simple_distance(t,0,0,1),5.,"unwrapped direct Euclidean distance");
    t.box_dimension[0]=80.;t.box_dimension[1]=1.;t.box_dimension[2]=4.;
    position(t,0,{{0.,0.,0.}});position(t,1,{{20.,.1,.3}});MATRIX vector(3,1);
    double expected=std::numeric_limits<double>::infinity();
    for(long nx=-30;nx<=30;++nx) for(long ny=-30;ny<=30;++ny) {
        const double x=20.+nx*80.+ny*7.,y=.1+ny;
        expected=std::min(expected,std::sqrt(x*x+y*y+.3*.3));
    }
    near(GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,vector,7.),expected,
         "F09 shortest lattice may require gradient images outside minus/plus one");
}

void extreme_geometry_domain() {
    Conditions c;TRAJECTORY t(c,2);MATRIX vector(3,1);
    const auto rejected=[&](const std::function<void()>& operation,const char* message) {
        bool caught=false;try { operation(); }catch(const std::domain_error&) { caught=true; }
        require(caught,message);
    };
    // Independent reviewer counterexamples: the former unchecked casts reached
    // +/-4.9e19, or requested 9.8e9 iterations, for this finite displacement.
    for(double ly:{1e-20,1e-10}) {
        t.box_dimension[0]=1.;t.box_dimension[1]=ly;t.box_dimension[2]=1.;
        position(t,0,{{0.,0.,0.}});position(t,1,{{.49,0.,0.}});
        rejected([&] { GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,vector,0.); },
                 "extreme shear aspect ratio must reject before cast or scan");
        Conditions configuration;configuration.set("MECHANICAL_PERTURBATION","SIMPLE_SHEAR");
        configuration.set("box_dimension",ly==1e-20 ? "1e-20" : "1e-10");
        configuration.set("LBx","1");configuration.set("LBy",ly==1e-20 ? "1e-20" : "1e-10");
        configuration.set("LBz","1");
        rejected([&] { RDIST distances(configuration); },
                 "unsupported shear box rejects before parallel distance allocation");
    }
    for(long k=0;k<3;++k)t.box_dimension[k]=1.;
    for(double y:{1e20,-1e20,std::numeric_limits<double>::max()}) {
        position(t,0,{{.1,y,.3}});
        rejected([&] { GEOMETRY::apply_shear_boundary_condition(t,0,0,1,.2); },
                 "huge finite shear winding must reject explicitly");
        rejected([&] { GEOMETRY::minimum_image_convention_loop(t,0); },
                 "huge finite periodic winding must reject explicitly");
    }
    position(t,0,{{.1,.2,.3}});position(t,1,{{.4,.5,.6}});
    rejected([&] { GEOMETRY::apply_shear_boundary_condition(t,0,0,1,
                       std::numeric_limits<double>::infinity()); },"nonfinite shear shift rejected");
    t.box_dimension[1]=std::numeric_limits<double>::quiet_NaN();
    rejected([&] { GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,vector,.2); },
             "nonfinite shear length rejected");
    t.box_dimension[1]=1.;position(t,1,{{std::numeric_limits<double>::infinity(),.2,.3}});
    rejected([&] { GEOMETRY::get_minimum_distance_simple_shear_fixed_axis(t,0,0,1,vector,.2); },
             "nonfinite shear pair displacement rejected");
}

void fractional_box_cell_geometry() {
    Conditions c;c.set("cell_list","TRUE");c.set("box_dimension",".8");
    c.set("LBx",".8");c.set("LBy",".8");c.set("LBz",".8");c.set("cutoff_connection",".2");
    TRAJECTORY trajectory(c,2);RDIST distances(c);
    require(distances.CELL_LIST_BOOST && distances.N_cells==64 && distances.N_neighbor_cells==27,
            "fractional box must retain 4x4x4 cell subdivision");
    for(long axis=0;axis<3;++axis) {
        require(distances.N_div[axis]==4,"division precedes integer conversion for L=.8/cutoff=.2");
        near(distances.cell_length[axis],.2,"fractional cell length");
    }
    const Vec samples[]={{{.01,.02,.03}},{{.79,.02,.03}},{{.01,.79,.03}},{{.01,.02,.79}}};
    for(long i=0;i<4;++i) position(trajectory,i,samples[i]);
    long scratch[3];distances.allocate_cells_from_positions(trajectory,0,scratch);
    for(long i=0;i<4;++i) GEOMETRY::compute_RDIST_particle(distances,i,trajectory,0);
    for(long i=0;i<4;++i) for(long j=0;j<4;++j) {
        double norm_squared=0.;
        for(long axis=0;axis<3;++axis) {
            const double displacement=samples[j][axis]-samples[i][axis];
            const double expected=displacement-.8*std::floor(displacement/.8+.5);
            near(distances.Rvec[i][j](axis),expected,"fractional cell/PBC vector vs direct oracle");
            norm_squared+=expected*expected;
        }
        near(distances.Rsca[i](j),std::sqrt(norm_squared),"fractional cell/PBC distance vs direct oracle");
    }
    // Invalid subdivision and product domains reject before owned allocation.
    for(const char* cutoff:{"1e-30","1e-9"}) {
        Conditions invalid;invalid.set("cell_list","TRUE");invalid.set("cutoff_connection",cutoff);
        bool rejected=false;try { CLIST cells(invalid); }catch(const std::domain_error&) { rejected=true; }
        require(rejected,"unrepresentable cell counts/products must reject before allocation");
    }
}

void report_duration_oracle() {
    Conditions c;TRAJECTORY trajectory(c,2);MATRIX energy(1,31,0.);
    BROWNIAN::BROWNIAN_VARIABLE brownian(c,2);
    REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE repulsive(c,2);
    HEUR::TEMPORAL_VARIABLE_HEUR network(c,2);
    network.time_SS=1e12;network.time_LV_init=1e12;
    network.time_LV_force=1e12;network.time_LV_update=1e12;
    network.N_associations=0;network.N_diff_associations=0;
    const std::vector<double> initial(trajectory.data,trajectory.data+trajectory.size);
    const auto check=[&](auto& phases,const auto& report) {
        // These synthetic cumulative times are deliberately unrelated to the
        // actual call duration. An independent wall-time bracket detects a
        // returned cumulative phase sum, without a fragile speed threshold.
        phases.time_LV=1e12;phases.time_AN=2e12;
        phases.time_file=3e12;phases.time_DIST=4e12;
        const double begin=omp_get_wtime();
        const double elapsed=report();
        const double end=omp_get_wtime();
        require(std::isfinite(elapsed) && elapsed>=0. && elapsed<=end-begin+1e-9,
                "report return must be elapsed call time, not cumulative phase totals");
        require(phases.time_LV==1e12 && phases.time_AN==2e12 &&
                phases.time_file==3e12 && phases.time_DIST==4e12,
                "report preserves printed cumulative phase totals");
        require(std::equal(initial.begin(),initial.end(),trajectory.data),
                "diagnostic report does not modify physical trajectory state");
        for(long i=0;i<energy.size;++i) near(energy(i),0.,"diagnostic report preserves energy state");
    };
    check(brownian,[&] { return BROWNIAN::report_simulation_info(trajectory,energy,brownian); });
    check(repulsive,[&] { return REPULSIVE_BROWNIAN::report_simulation_info(trajectory,energy,repulsive); });
    check(network,[&] { return HEUR::report_simulation_info(trajectory,energy,network); });
}

void trajectory_deep_copy() {
    Conditions c;TRAJECTORY source(c,2);source.box_dimension[1]=17.;source.c_t=11;
    source(0)=1.2;position(source,0,{{1.,2.,3.}});source.Method="TEST";
    TRAJECTORY copy(source), assigned(c,2);assigned=source;assigned=assigned;
    for(TRAJECTORY* target:{&copy,&assigned}) {
        require(target->N_dimension==3 && target->rows==2 && target->Np==4 && target->c_t==11,
                "F11 trajectory metadata survives copy");
        require(target->Method=="TEST","F11 copy model identity");
        near(target->box_dimension[1],17.,"F11 anisotropic box copied");
        near((*target)(0,0,2),3.,"F11 trajectory samples copied");
        target->box_dimension[1]=23.;(*target)(0,0,2)=7.;
    }
    near(source.box_dimension[1],17.,"F11 copied box has independent ownership");
    near(source(0,0,2),3.,"F11 copied samples have independent ownership");
}

struct Test { const char* name; std::function<void()> run; };
} // namespace

int main() {
    const Test tests[]={
        {"coupled_energy_force_derivative",coupled_energy_and_force},
        {"fractional_cutoff_and_move_weight",fractional_cutoff_and_move_weight},
        {"renewal_and_saturated_tracking",renewal_and_saturated_tracking},
        {"heterogeneous_particle_friction",heterogeneous_particle_friction},
        {"serialized_topology_and_renewal",serialized_topology_and_renewal},
        {"shear_supported_domain_and_raw_distance",shear_domain_and_raw_distance},
        {"trajectory_copy_ownership",trajectory_deep_copy},
        {"extreme_geometry_domain",extreme_geometry_domain},
        {"report_elapsed_duration",report_duration_oracle},
        {"fractional_box_cell_geometry",fractional_box_cell_geometry},
        {"ordinary_periodic_geometry",periodic_geometry},
        {"cubic_xy_shear_lattice",shear_geometry},
        {"step_affine_and_box_shift",step_shear_and_shift},
        {"single_gradient_crossings",single_gradient_crossings},
        {"connector_energy_derivatives",connector_laws},
        {"kinetic_acceptance_laws",kinetic_probabilities},
        {"serial_topology_and_labelled_chains",topology_actions},
        {"multiplicity_suggestion_cdf",suggestion_distribution},
        {"network_shear_force_energy_virial",network_force_and_virial},
        {"multi_box_gradient_crossing",known_multibox_shear},
        {"elongated_shear_box",known_elongated_shear},
        {"gaussian_minimum_R0",known_minimum_R0_gaussian},
        {"fene_minimum_R0",known_minimum_R0_fene},
        {"non_neighbor_bond_cache",known_stale_bond_geometry}
    };
    int passed=0,failed=0;
    for(const auto& test:tests) {
        try {
            test.run(); ++passed; std::cout << "PASS " << test.name << '\n';
        } catch(const std::exception& error) {
            ++failed;
            std::cerr << "FAIL " << test.name << ": " << error.what() << '\n';
        }
    }
    std::cout << "Legacy physics: " << passed << " passed, " << failed << " failed\n";
    return failed ? 1 : 0;
}
