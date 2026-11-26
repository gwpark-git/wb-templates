#include "bdmc/execution.hpp"
#include "bdmc/reference.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool x,const char* message){if(!x)throw std::runtime_error(message);}
void close(double x,double y,double tol,const char* message){check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol,message);}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::invalid_argument&){failed=true;}check(failed,"invalid selection/input accepted");}
void same(const std::vector<double>& a,const std::vector<double>& b,double tolerance){check(a.size()==b.size(),"array shape");for(std::size_t i=0;i<a.size();++i)close(a[i],b[i],tolerance,"array mismatch");}
void run(const std::string& backend,int threads){
 auto kernel=bdmc::make_execution_kernels(backend,threads);
 const std::vector<double> x{9.8,.5,.3,.3,.7,.4,1.1,1.,.9};
 bdmc::RepulsiveSpec s;s.n_particles=3;s.sigma=1.4;s.epsilon=2.3;s.force_search=bdmc::SpatialBackend::BruteForce;
 bdmc::RepulsiveSimulation serial(s,x);const auto expected=serial.observe();
 const auto result=kernel->repulsive(x,10,2.3,1.4);same(result.forces,expected.forces,5e-14);
 double energy=0;for(double u:result.energies)energy+=u;close(energy,expected.potential_energy,1e-14,"pair counting");
 for(std::size_t k=0;k<9;++k){double virial=0;for(std::size_t i=0;i<3;++i)virial+=result.pair_virials[9*i+k];close(virial,expected.pair_virial[k],2e-14,"pair virial");}
 const std::vector<double> q{.2,-.4,.1,.3,.2,-.1};
 for(bool fene:{false,true}){
  const auto spring=kernel->springs(q,1.7,3,fene);
  for(std::size_t i=0;i<2;++i){const double r2=q[3*i]*q[3*i]+q[3*i+1]*q[3*i+1]+q[3*i+2]*q[3*i+2];const double scale=fene?1.7/(1-r2/9):1.7;
   for(std::size_t c=0;c<3;++c)close(spring.forces[3*i+c],scale*q[3*i+c],1e-14,"spring factor");
   close(spring.energies[i],fene?-.5*1.7*9*std::log1p(-r2/9):.5*1.7*r2,1e-14,"spring energy");
  }
 }
 const auto step=kernel->linear_step(q,q,q,-.1,.3);for(std::size_t k=0;k<q.size();++k)close(step[k],1.2*q[k],1e-15,"linear update");
 rejects([&]{kernel->repulsive({0,0,0,0,0,0},10,1,1);});
 rejects([&]{kernel->springs({3,0,0},1,3,true);});
 rejects([&]{kernel->linear_step(q,{0,0,0},q,1,1);});
 if(backend!="serial"){
  auto parallel=s;parallel.execution_backend=backend;parallel.backend_threads=threads;
  bdmc::RepulsiveSimulation other(parallel,x);serial.advance(20);other.advance(20);same(serial.observe().unwrapped_positions,other.observe().unwrapped_positions,3e-12);
  for(const auto protocol:{bdmc::FlowProtocol::NoFlow,bdmc::FlowProtocol::SimpleShear,bdmc::FlowProtocol::OscillatoryShear}){
   bdmc::DumbbellSpec d;d.n_dumbbells=2;d.protocol=protocol;d.shear_rate=-.8;d.strain_amplitude=.3;d.angular_frequency=2.;d.seed=11517;
   auto p=d;p.execution_backend=backend;p.backend_threads=threads;
   bdmc::DumbbellSimulation a(d,q),b(p,q);a.advance(100);b.advance(100);same(a.observe().connectors,b.observe().connectors,1e-12);
  }
  bdmc::DumbbellSpec d;d.n_dumbbells=1;d.spring_law=bdmc::SpringLaw::FENE;d.q0=.7;d.dt=.05;d.seed=935003;d.execution_backend=backend;d.backend_threads=threads;
  bdmc::DumbbellSimulation failure(d);rejects([&]{failure.advance(1);});rejects([&]{failure.advance(1);});check(failure.observe().step==0,"failed step consumed clock/RNG");
 }
}
}
int main(){try{
 run("serial",1);
 for(const auto& c:bdmc::execution_backends()){
  if(c.name=="openmp"&&c.available){run("openmp",1);run("openmp",2);std::cout<<"PASS openmp parity/transactions threads1/2\n";}
  if(!c.available)rejects([&]{bdmc::make_execution_kernels(c.name);});
 }
 rejects([]{bdmc::make_execution_kernels("gpu_auto");});rejects([]{bdmc::make_execution_kernels("serial",2);});
 std::cout<<"PASS bulk force/energy/virial, integration, domains and explicit backend availability\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
