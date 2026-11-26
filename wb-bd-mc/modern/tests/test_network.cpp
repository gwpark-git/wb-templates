#include "bdmc/network.hpp"

#ifdef BDMC_NETWORK_LEGACY_ORACLE
#include "../../lib/potential.h"
#endif

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace {
using namespace bdmc;
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
void close(double actual, double expected, const char* message, double tolerance=2e-14) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance)
        throw std::runtime_error(std::string(message)+" actual="+std::to_string(actual)+" expected="+std::to_string(expected));
}
template<class F> void rejects(F fn, const char* message) {
    bool rejected=false; try {fn();} catch (const std::invalid_argument&) {rejected=true;}
    check(rejected,message);
}
void independent_counts(const NetworkTopology& topology) {
    const auto n=topology.n_particles();
    std::vector<std::size_t> ends(n), dense(n*n);
    std::size_t bridges=0, pairs=0;
    for (const auto& chain : topology.endpoints()) {
        ++ends[chain[0]]; ++ends[chain[1]];
        ++dense[chain[0]*n+chain[1]]; ++dense[chain[1]*n+chain[0]];
        bridges += chain[0]!=chain[1];
    }
    for (std::size_t i=0;i<n;++i) {
        std::size_t row_sum=0;
        check(topology.adjacency()[i].front().particle==i,"permanent self comes first");
        check(topology.weight(i,i)%2==0,"even loop-end weight");
        for (std::size_t j=0;j<n;++j) {
            check(topology.weight(i,j)==dense[i*n+j],"independent endpoint-to-dense reconstruction");
            row_sum+=dense[i*n+j];
            if (i<j && dense[i*n+j]) ++pairs;
        }
        check(row_sum==ends[i] && ends[i]==topology.attached_ends()[i],"local attached-end conservation");
    }
    check(std::accumulate(ends.begin(),ends.end(),std::size_t{0})==2*topology.n_chains(),"global 2Nc conservation");
    check(topology.n_bridges()==bridges && topology.n_distinct_pairs()==pairs,"multiplicity differs from distinct pairs");
    topology.validate();
}
void owned_topology_and_transactions() {
    auto t=NetworkTopology::loops(3,2);
    check(t.n_chains()==6 && t.weight(0,0)==4,"loop initialization owns all chains");
    check(t.particle_of_end(0)==0 && t.opposite_end(0)==6 && t.opposite_end(6)==0,"legacy stable end IDs");
    t.transact({{0,0,1},{1,0,1}});
    check(t.weight(0,1)==2 && t.n_bridges()==2 && t.n_distinct_pairs()==1,"duplicate physical bridges counted");
    independent_counts(t);
    const auto before=t.endpoints();const auto epoch=t.epoch();
    rejects([&]{t.transact({{0,1,2},{1,2,0}});},"late expectation failure rejects whole batch");
    check(t.endpoints()==before && t.epoch()==epoch,"failed batch preserves endpoints and epoch");
    rejects([&]{t.transact({{0,1,2},{99,0,1}});},"invalid end after valid event rejects whole batch");
    check(t.endpoints()==before && t.epoch()==epoch,"invalid-end rollback preserves caches");
    t.transact({{0,1,2},{0,2,1}});
    check(t.epoch()==epoch,"net-zero transaction does not invalidate topology");
    t.transact({{0,1,0},{1,1,0}});independent_counts(t);
    check(t.n_bridges()==0 && t.weight(0,0)==4,"deletion returns exactly to loops");
    rejects([]{NetworkTopology invalid(0,{});},"zero particles rejected");
    rejects([]{NetworkTopology invalid(2,{{0,2}});},"bad particle endpoint rejected");
    NetworkTopology empty(2,{}); independent_counts(empty);
    rejects([&]{empty.particle_of_end(0);},"empty network rejects end without modulo zero");
}
void action_policy_and_attempts() {
    NetworkTopology t(3,{{0,0},{0,0},{1,1}});
    NetworkConstraints c;
    NetworkKinetics k;
    auto attempt=[&](std::size_t anchor,std::size_t other,std::size_t dest,std::size_t end) {
        NetworkProposal p;p.anchor=anchor;p.other=other;p.destination=dest;p.end=end;p.topology_epoch=t.epoch();
        return attempt_network_transition(t,p,0,.2,c,k);
    };
    check(attempt(0,0,1,0).action==NetworkAction::Add,"loop to bridge action");
    check(attempt(0,0,1,1).action==NetworkAction::Add,"second bridge action");
    check(attempt(0,1,2,0).action==NetworkAction::Move,"move opposite bridge end");
    check(t.endpoints()[0]==ChainEndpoints{2,0},"anchor retained, opposite end moved");
    check(attempt(0,1,1,1).rejection==NetworkRejection::SameDestination,"unchanged destination cancels");
    check(attempt(0,1,0,1).action==NetworkAction::OppositeDelete,"bridge to loop action");
    check(attempt(0,2,0,0).action==NetworkAction::OppositeDelete,"last bridge to loop");
    independent_counts(t);
    c.min_attached_ends=4;c.max_attached_ends=4;
    check(classify_network_action(t,0,0,1,c).rejection==NetworkRejection::CountOrMultiplicity,"count constraint boundary cancels");
    c=NetworkConstraints{}; c.multiple_connections=false;
    check(attempt(0,0,1,0).action==NetworkAction::Add,"single policy first bridge");
    check(attempt(0,0,1,1).rejection==NetworkRejection::CountOrMultiplicity,"single policy rejects duplicate without retry");
    c.min_attached_ends=100;c.max_attached_ends=100;
    check(attempt(0,1,0,0).action==NetworkAction::OppositeDelete,"single policy ignores endpoint count bounds, as legacy");
    c=NetworkConstraints{};
    NetworkProposal p;p.anchor=0;p.other=0;p.destination=2;p.end=0;p.topology_epoch=t.epoch();
    const auto before=t.endpoints();k.law=NetworkTransitionLaw::Dissociation;k.delta_t0=.25;
    check(attempt_network_transition(t,p,0,.25,c,k).rejection==NetworkRejection::Acceptance,"strict roll<p boundary");
    check(t.endpoints()==before,"kinetic rejection never mutates state");
    check(attempt_network_transition(t,p,0,.249,c,k).action==NetworkAction::Add,"below acceptance boundary succeeds");
    rejects([&]{attempt_network_transition(t,p,0,.1,c,k);},"stale proposal rejected");
    p.topology_epoch=t.epoch();
    rejects([&]{attempt_network_transition(t,p,0,1,c,k);},"roll one invalid");
    rejects([&]{attempt_network_transition(t,p,0,-.1,c,k);},"negative roll invalid");
}
void pure_laws_and_schedule() {
    NetworkKinetics k;k.l_cap=.4;k.energy_barrier=2;
    for (double tension : {-2.,0.,.7,9.}) for (double dt : {0.,.001,.1,2.}) {
        k.delta_t0=dt;const double rate=std::exp(tension*k.l_cap);
        k.law=NetworkTransitionLaw::Uniform;close(network_transition_probability(tension,k),1,"uniform law");
        k.law=NetworkTransitionLaw::Weighted;close(network_transition_probability(tension,k),std::min(1.,rate),"weighted law");
        k.law=NetworkTransitionLaw::Metropolis;close(network_transition_probability(tension,k),std::min(1.,std::exp(tension*k.l_cap-2)),"Metropolis old-tension barrier law");
        k.law=NetworkTransitionLaw::Dissociation;close(network_transition_probability(tension,k),std::min(1.,dt*rate),"saturated linear-rate law");
        k.law=NetworkTransitionLaw::FirstOrder;close(network_transition_probability(tension,k),1-std::exp(-dt*rate),"first-order integrated hazard");
        close(network_detachment_weight(tension,k.l_cap,NetworkDetachmentLaw::Weighted),rate,"weighted detachment is not clipped");
    }
    k.law=NetworkTransitionLaw::FirstOrder;k.delta_t0=1e-20;k.l_cap=0;
    close(network_transition_probability(1,k),1e-20,"tiny positive first-order probability retained",1e-35);
    k.law=NetworkTransitionLaw::Dissociation;k.l_cap=1;k.delta_t0=0;
    close(network_transition_probability(1000,k),0,"zero duration with huge finite tension gives zero");
    k.delta_t0=.1;rejects([&]{network_transition_probability(1000,k);},"overflowing kinetic intermediate rejected explicitly");
    for(auto law:{NetworkTransitionLaw::Dissociation,NetworkTransitionLaw::FirstOrder}) {
        k.law=law;k.delta_t0=1e-320;
        rejects([&]{network_transition_probability(720,k);},"overflowing rate with tiny duration cannot silently return one");
        k.delta_t0=1e308;
        rejects([&]{network_transition_probability(-750,k);},"underflowing rate with huge duration cannot silently return zero");
    }
    rejects([&]{network_detachment_weight(1000,1,NetworkDetachmentLaw::Weighted);},"nonfinite suggestion weight fails explicitly");
    NetworkSchedule s;s.dt_over_tauR=.003;s.steps_per_block=7;s.repulsion_coefficient=2;s.Rt=5;
    close(s.delta_t0(),.0021,"legacy block-to-rate conversion");
    check(s.attempts_per_block(NetworkTopology::loops(3,2))==6,"legacy ordinary Nchains attempt budget");
    s.steps_per_block=0;rejects([&]{s.delta_t0();},"zero topology block interval rejected");
}
void finite_distributions() {
    const auto p=network_discrete_distribution(3,{{0,1},{2,3},{1,0}});
    close(p[0].probability,.25,"finite candidate measure has no radial Jacobian");close(p[1].probability,.75,"finite candidate normalized");
    rejects([]{network_discrete_distribution(2,{{0,0},{1,0}});},"allzero candidates fail instead of resampling");
    rejects([]{network_discrete_distribution(2,{{0,1},{0,2}});},"duplicate candidate ID rejected");
    rejects([]{network_discrete_distribution(2,{{0,-1},{1,2}});},"negative candidate rejected");
    NetworkTopology t(3,{{0,0},{0,1},{1,0},{0,2}});
    const auto s=network_suggestion_distribution(t,0,{{0,2},{1,3},{2,5}});
    // Independent count: two loop ends*2, two bridge ends*3, one*5.
    close(s[0].probability,4./15,"loop end degeneracy in suggestion");
    close(s[1].probability,6./15,"bridge multiplicity in suggestion");
    close(s[2].probability,5./15,"weighted suggestion denominator");
    rejects([&]{network_suggestion_distribution(t,0,{{0,1},{1,1}});},"missing existing old neighbor fails");
}
std::vector<ChainEndpoints> decoded(std::size_t state,std::size_t n,std::size_t chains) {
    std::vector<ChainEndpoints> e(chains);
    for (auto& row:e) for (auto& particle:row) {particle=state%n;state/=n;}
    return e;
}
std::size_t encoded(const std::vector<ChainEndpoints>& endpoints,std::size_t n) {
    std::size_t value=0,place=1;
    for (const auto& row:endpoints) for (const auto particle:row) {value+=particle*place;place*=n;}
    return value;
}
void exhaustive_finite_proposal_oracle() {
    for (std::size_t n:{2,3}) for (std::size_t chains:{1,2}) {
        std::size_t states=1;for(std::size_t i=0;i<2*chains;++i)states*=n;
        for (std::size_t state=0;state<states;++state) {
            NetworkTopology t(n,decoded(state,n,chains));independent_counts(t);
            std::vector<std::vector<NetworkWeight>> a(n),d(n),f(n);
            for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
                a[i].push_back({j,static_cast<double>(j+1)});d[i].push_back({j,static_cast<double>(i+j+1)});f[i].push_back({j,.2*static_cast<double>(i+j)});
            }
            NetworkConstraints c;NetworkKinetics k;k.law=NetworkTransitionLaw::Dissociation;k.delta_t0=.3;k.l_cap=.4;
            const auto proposals=enumerate_network_proposals(t,a,d,f,c,k);
            std::map<std::array<std::size_t,4>,double> actual,expected;
            double total=0;
            for(const auto& p:proposals) {
                total+=p.proposal_probability;actual[{p.anchor,p.other,p.destination,p.end}]+=p.proposal_probability;
                auto next=t;const auto r=attempt_network_transition(next,p,.2*static_cast<double>(p.anchor+p.other),0,c,k);
                auto expected_ends=t.endpoints();
                if(p.other!=p.destination)expected_ends[p.end%chains][p.end/chains]=p.destination;
                check(next.endpoints()==expected_ends,"all finite labeled transitions move exactly one opposite end");
                independent_counts(next);
                close(p.mutation_probability,p.action==NetworkAction::Cancel?0:p.proposal_probability*p.acceptance_probability,"accepted mass retains cancellation");
                check(r.action==p.action,"enumerated classification agrees with applied proposal");
            }
            // Oracle enumerates individual anchor-attached ends directly,
            // without traversing the implementation's aggregate adjacency.
            const auto& endpoint=t.endpoints();std::size_t occupied=0;
            for(std::size_t i=0;i<n;++i)occupied+=t.attached_ends()[i]!=0;
            const double destination_Z=static_cast<double>(n*(n+1))/2;
            for(std::size_t anchor=0;anchor<n;++anchor)if(t.attached_ends()[anchor]) {
                double old_Z=0;
                for(std::size_t e=0;e<2*chains;++e) {
                    const auto at=endpoint[e%chains][e/chains];
                    if(at==anchor)old_Z+=static_cast<double>(anchor+endpoint[e%chains][1-e/chains]+1);
                }
                for(std::size_t fixed=0;fixed<2*chains;++fixed)if(endpoint[fixed%chains][fixed/chains]==anchor) {
                    const auto moving=fixed<chains?fixed+chains:fixed-chains;
                    const auto other=endpoint[moving%chains][moving/chains];
                    for(std::size_t dest=0;dest<n;++dest)
                        expected[{anchor,other,dest,moving}]+=static_cast<double>(anchor+other+1)/old_Z/
                            static_cast<double>(occupied)*static_cast<double>(dest+1)/destination_Z;
                }
            }
            close(total,1,"exhaustive proposal mass sums to one");check(actual.size()==expected.size(),"exhaustive event possibility count");
            for(const auto& item:expected)close(actual[item.first],item.second,"independent labeled-end event probability");
        }
    }
}
void one_chain_uniform_detailed_balance() {
    // This deliberately tiny UNIFORM/no-count-restriction subset has the
    // uniform distribution on labeled endpoint states. No general DB claim.
    for(std::size_t n:{2,3}) {
        const auto states=n*n;std::vector<double> matrix(states*states);
        std::vector<std::vector<NetworkWeight>> a(n),d(n),f(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {a[i].push_back({j,1});d[i].push_back({j,1});f[i].push_back({j,0});}
        for(std::size_t state=0;state<states;++state) {
            NetworkTopology t(n,decoded(state,n,1));
            for(const auto& p:enumerate_network_proposals(t,a,d,f,NetworkConstraints{},NetworkKinetics{})) {
                auto next=t;attempt_network_transition(next,p,0,0,NetworkConstraints{},NetworkKinetics{});
                matrix[state*states+encoded(next.endpoints(),n)]+=p.proposal_probability;
            }
        }
        for(std::size_t i=0;i<states;++i) {
            double total=0;
            for(std::size_t j=0;j<states;++j) {total+=matrix[i*states+j];close(matrix[i*states+j],matrix[j*states+i],"one-chain uniform detailed balance");}
            close(total,1,"one-chain transition row normalized");
        }
    }
}
void independent_bond_geometry_and_energy() {
    NetworkTopology t(4,{{0,3},{3,0},{1,1},{1,2}});
    std::vector<double> x{.1,.2,.3, 1,2,3, 1.3,2.4,3.2, 6.1,.2,.3};
    NetworkConnectorSpec s;s.spring_constant=2.4;s.scale=.25;s.q0=8;
    for(auto law:{NetworkConnectorLaw::Hookean,NetworkConnectorLaw::FENE}) {
        s.law=law;const auto f=evaluate_network_connectors(t,x,{20,20,20},s);
        check(std::abs(f.forces[0])>1,"bond beyond any short neighbor cutoff still contributes");
        for(std::size_t c=0;c<3;++c) {
            double sum=0;for(std::size_t i=0;i<t.n_particles();++i)sum+=f.forces[i*3+c];close(sum,0,"Newton force conservation");
        }
        constexpr double h=1e-5;
        for(std::size_t c=0;c<x.size();++c) {
            auto p=x,m=x;p[c]+=h;m[c]-=h;
            const auto derivative=(evaluate_network_connectors(t,p,{20,20,20},s).potential_energy-
                evaluate_network_connectors(t,m,{20,20,20},s).potential_energy)/(2*h);
            close(f.forces[c],-derivative,"dynamic multiplicity force is negative energy gradient",2e-9);
        }
        double expected_energy=0;std::array<double,9> virial{};
        // Independent physical-chain traversal, counting both duplicate chains.
        for(const auto& ends:t.endpoints())if(ends[0]!=ends[1]) {
            std::array<double,3> q{};double r2=0;
            for(std::size_t c=0;c<3;++c) {q[c]=x[ends[1]*3+c]-x[ends[0]*3+c];r2+=q[c]*q[c];}
            const double factor=s.spring_constant*s.scale/(law==NetworkConnectorLaw::FENE?1-r2/(s.q0*s.q0):1);
            expected_energy+=law==NetworkConnectorLaw::FENE?-.5*s.spring_constant*s.scale*s.q0*s.q0*std::log(1-r2/(s.q0*s.q0)):.5*s.spring_constant*s.scale*r2;
            for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b)virial[a*3+b]+=q[a]*factor*q[b];
        }
        close(f.potential_energy,expected_energy,"physical-chain energy oracle");
        for(std::size_t c=0;c<9;++c)close(f.connector_virial[c],virial[c],"physical-chain connector virial oracle");
    }
    NetworkTopology pair(2,{{0,1}});s.law=NetworkConnectorLaw::Hookean;s.spring_constant=3;s.scale=1;
    const auto periodic=evaluate_network_connectors(pair,{9.8,0,0,.3,0,0},{10,10,10},s);
    close(periodic.forces[0],1.5,"periodic bonded geometry orientation");close(periodic.forces[3],-1.5,"periodic attractive partner");
    s.law=NetworkConnectorLaw::FENE;s.q0=.5;
    rejects([&]{evaluate_network_connectors(pair,{0,0,0,.5,0,0},{10,10,10},s);},"FENE exact boundary fails explicitly");
    check(pair.endpoints()==std::vector<ChainEndpoints>{{0,1}},"force-domain failure leaves topology untouched");
}
#ifdef BDMC_NETWORK_LEGACY_ORACLE
void actual_legacy_header_laws() {
    POTENTIAL_VARIABLE p;p.constructor_POTENTIAL_VARIABLE();p.N_dimension=3;p.scale_factor_chain=1.2;
    p.ratio_RM_R0=4;p.l_cap=.4;p.energy_barrier=2;
    auto force=FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force;
    using Legacy=double(*)(double,double(*)(double,POTENTIAL_VARIABLE&),POTENTIAL_VARIABLE&);
    const std::pair<NetworkTransitionLaw,Legacy> laws[]={{NetworkTransitionLaw::Uniform,KINETICS::UNIFORM::transition_probability},
        {NetworkTransitionLaw::Metropolis,KINETICS::METROPOLIS::transition_probability},
        {NetworkTransitionLaw::Weighted,KINETICS::WEIGHTED::transition_probability},
        {NetworkTransitionLaw::Dissociation,KINETICS::dissociation_probability},
        {NetworkTransitionLaw::FirstOrder,KINETICS::FIRST_ORDER::dissociation_probability}};
    NetworkKinetics k;k.l_cap=p.l_cap;k.energy_barrier=p.energy_barrier;
    for(double dt:{0.,.001,.1,2.})for(double r:{0.,.2,1.,1.7})for(const auto& law:laws) {
        p.delta_t0=dt;k.delta_t0=dt;k.law=law.first;
        close(network_transition_probability(force(r,p),k),law.second(r,force,p),"actual legacy inline transition law");
        close(network_detachment_weight(force(r,p),p.l_cap,NetworkDetachmentLaw::Weighted),KINETICS::WEIGHTED::detachment_weight(r,force,p),"actual legacy detachment law");
    }
    NetworkTopology pair(2,{{0,1}});NetworkConnectorSpec s;s.spring_constant=3*1.2*1.2;s.q0=4;s.scale=.25;
    for(double r:{.2,1.,1.7,3.7})for(auto law:{NetworkConnectorLaw::Hookean,NetworkConnectorLaw::FENE}) {
        s.law=law;const auto frame=evaluate_network_connectors(pair,{0,0,0,r,0,0},{20,20,20},s);
        const double tension=law==NetworkConnectorLaw::FENE?FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_force(r,p):force(r,p);
        const double energy=law==NetworkConnectorLaw::FENE?FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_potential(r,p):FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential(r,p);
        close(frame.forces[0],s.scale*tension,"actual legacy connector radial force");
        close(frame.potential_energy,s.scale*energy,"actual legacy connector energy");
    }
}
#endif
} // namespace
int main() {
    const std::pair<const char*,std::function<void()>> tests[]={{"owned_topology_and_transactions",owned_topology_and_transactions},
        {"action_policy_and_attempts",action_policy_and_attempts},{"pure_laws_and_schedule",pure_laws_and_schedule},
        {"finite_distributions",finite_distributions},{"exhaustive_finite_proposal_oracle",exhaustive_finite_proposal_oracle},
        {"one_chain_uniform_detailed_balance",one_chain_uniform_detailed_balance},{"independent_bond_geometry_and_energy",independent_bond_geometry_and_energy},
#ifdef BDMC_NETWORK_LEGACY_ORACLE
        {"actual_legacy_header_laws",actual_legacy_header_laws},
#endif
    };
    for(const auto& test:tests)try {test.second();std::cout<<"PASS "<<test.first<<'\n';}
        catch(const std::exception& error){std::cerr<<"FAIL "<<test.first<<": "<<error.what()<<'\n';return 1;}
    std::cout<<"PASS "<<std::size(tests)<<" network checks\n";return 0;
}
