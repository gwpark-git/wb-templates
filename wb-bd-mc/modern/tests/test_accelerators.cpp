#include "bdmc/execution.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool condition,const std::string& message) {
    if(!condition) throw std::runtime_error(message);
}
void near(const std::vector<double>& actual,const std::vector<double>& expected,bool single,
          const std::string& label) {
    check(actual.size()==expected.size(),label+" shape");
    const double absolute=single ? 2e-6 : 2e-13;
    const double relative=single ? 3e-5 : 3e-13;
    for(std::size_t i=0;i<actual.size();++i)
        check(std::isfinite(actual[i]) && std::abs(actual[i]-expected[i])<=absolute+relative*std::abs(expected[i]),
              label+" differs at index "+std::to_string(i));
}
void near(const bdmc::KernelResult& actual,const bdmc::KernelResult& expected,bool single,const std::string& label) {
    near(actual.forces,expected.forces,single,label+" forces");
    near(actual.energies,expected.energies,single,label+" energies");
    near(actual.pair_virials,expected.pair_virials,single,label+" virials");
}
void rejects(const std::function<void()>& operation,const std::string& label) {
    try { operation(); }
    catch(const std::invalid_argument&) { return; }
    throw std::runtime_error(label+" did not reject invalid physics/input");
}
void exercise(bdmc::ExecutionKernels& kernels) {
    auto serial=bdmc::make_execution_kernels("serial");
    const bool single=kernels.precision()=="float32";
    const std::vector<double> pair{9.75,0,0,0.25,0,0,5,4,3};
    auto result=kernels.repulsive(pair,10,2,1);
    near(result,serial->repulsive(pair,10,2,1),single,"PBC pair");
    near(result.forces,{-1.5,0,0,1.5,0,0,0,0,0},single,"analytic PBC forces");
    near(result.energies,{5.0/24,5.0/24,0},single,"half-directed energy");
    std::vector<double> virial(27,0); virial[0]=virial[9]=0.375;
    near(result.pair_virials,virial,single,"half-directed virial");

    // Non-axis-aligned pairs catch transposition/sign errors in tensor output;
    // changing counts also exercises reusable device-buffer growth and reuse.
    for(std::size_t n:{1u,73u,9u,587u}) {
        std::vector<double> positions(3*n);
        for(std::size_t i=0;i<n;++i) {
            positions[3*i]=0.13+0.79*(i%9);
            positions[3*i+1]=0.21+0.83*((i/9)%9);
            positions[3*i+2]=0.17+0.73*(i/81);
        }
        near(kernels.repulsive(positions,10,1.3,1.3),serial->repulsive(positions,10,1.3,1.3),single,"allpairs count "+std::to_string(n));
    }
    const double h=1.7,q0=2;
    for(bool fene:{false,true}) {
        const std::vector<double> q{0.3,-0.4,0.2,0,0,0,-0.7,0.1,0.4};
        const auto spring=kernels.springs(q,h,q0,fene);
        near(spring,serial->springs(q,h,q0,fene),single,fene ? "FENE" : "Hookean");
        std::vector<double> force(q.size()),energy(q.size()/3);
        for(std::size_t i=0;i<energy.size();++i) {
            double norm2=0; for(int a=0;a<3;++a) norm2+=q[3*i+a]*q[3*i+a];
            const double ratio=fene ? norm2/(q0*q0):0;
            for(int a=0;a<3;++a) force[3*i+a]=h*q[3*i+a]/(1-ratio);
            energy[i]=fene ? -0.5*h*q0*q0*std::log1p(-ratio):0.5*h*norm2;
        }
        near(spring.forces,force,single,"analytic spring force");
        near(spring.energies,energy,single,"analytic spring energy");
    }
    std::vector<double> q(3*587),force(q.size()),noise(q.size());
    for(std::size_t k=0;k<q.size();++k) {
        q[k]=0.04*static_cast<double>(static_cast<int>(k%19)-9);
        force[k]=0.02*static_cast<double>(static_cast<int>(k%13)-6);
        noise[k]=0.13*static_cast<double>(static_cast<int>(k%7)-3);
    }
    near(kernels.springs(q,h,q0,true),serial->springs(q,h,q0,true),single,"FENE multiple threadgroups");
    const auto step=kernels.linear_step(q,force,noise,-0.07,0.23);
    std::vector<double> expected(q.size());
    for(std::size_t k=0;k<q.size();++k) expected[k]=q[k]-0.07*force[k]+0.23*noise[k];
    near(step,expected,single,"host-supplied noise Euler update");
    near(kernels.linear_step(q,force,noise,-0.07,0.23),step,single,"repeatability without internal RNG");

    rejects([&] { kernels.repulsive({0,0,0,0,0,0},10,1,1); },"coincidence");
    rejects([&] { kernels.repulsive({-0.1,0,0},10,1,1); },"unwrapped coordinate");
    rejects([&] { kernels.repulsive({0,0,0},10,1,6); },"cutoff");
    rejects([&] { kernels.springs({2,0,0},h,2,true); },"FENE boundary");
    rejects([&] { kernels.springs({2.1,0,0},h,2,true); },"FENE exterior");
    rejects([&] { kernels.springs({0,0,0},-1,2,false); },"spring constant");
    rejects([&] { kernels.springs({0,0},h,2,false); },"malformed connector");
    rejects([&] { kernels.springs({},h,2,false); },"empty connector");
    rejects([&] { kernels.springs({std::numeric_limits<double>::quiet_NaN(),0,0},h,2,false); },"NaN connector");
    rejects([&] { kernels.linear_step({0,0,0},{0,0,0},{0,0},1,1); },"noise shape");
    rejects([&] { kernels.linear_step({0,0,0},{0,0,0},{0,0,0},std::numeric_limits<double>::infinity(),1); },"infinite scale");
    if(single) {
        rejects([&] { kernels.springs({1e40,0,0},h,2,false); },"float32 range");
        rejects([&] { kernels.springs({1e-40,0,0},h,2,false); },"float32 subnormal input");
        rejects([&] { kernels.springs({2*(1-1e-9),0,0},h,2,true); },"float32 FENE rounded boundary");
        rejects([&] { kernels.repulsive({10-1e-9,0,0},10,1,1); },"float32 rounded wrapping boundary");
    }
    // Device errors must not poison a later valid call or change host inputs.
    near(kernels.repulsive(pair,10,2,1),result,single,"recovery after rejected inputs");
}
} // namespace

int main(int argc,char** argv) {
    if(argc!=2) { std::cerr<<"usage: bdmc_accelerator_tests cuda|metal|accelerate\n"; return 2; }
    const std::string name=argv[1];
    try {
        const auto capabilities=bdmc::execution_backends();
        const auto found=std::find_if(capabilities.begin(),capabilities.end(),[&](const auto& c) { return c.name==name; });
        check(found!=capabilities.end(),"backend is not registered: "+name);
        if(!found->available) {
            std::cout<<"SKIP "<<name<<": compiled="<<found->compiled<<"; "<<found->reason<<"\n";
            return 77;
        }
        auto kernels=bdmc::make_execution_kernels(name);
        check(kernels->name()==name,"factory changed requested backend");
        check(kernels->precision()==found->precision,"factory precision differs from capability");
        exercise(*kernels);
        std::cout<<"PASS "<<name<<" "<<kernels->precision()<<": analytic force/energy/virial, CPU parity, Euler, domain and buffer reuse\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"FAIL "<<name<<": "<<error.what()<<"\n";
        return 1;
    }
}
