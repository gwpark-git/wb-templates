#include "../lib/analysis.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void check(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
template<class F> void rejects(F&& f) {
    bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"invalid legacy autocorrelation contract was accepted");
}
}
int main() {
    try {
        using namespace POST_PROCESSING;
        MATRIX x(4,1,0.),serial(4,2,19.),parallel(4,2,-77.);
        for(long i=0;i<4;++i) x(i,0)=i+1.;
        compute_autocorrelation_brute_force(x,0,serial,1,CorrelationWindow::AllOrigins);
        compute_autocorrelation_brute_force_OMP(x,0,2,parallel,1,CorrelationWindow::AllOrigins);
        const double expected[]={7.5,20./3,5.5,4.};
        for(long i=0;i<4;++i) {
            check(std::abs(serial(i,1)-expected[i])<1e-13,"independent all-origins oracle");
            check(serial(i,1)==parallel(i,1),"serial/OpenMP exact identity");
            check(serial(i,0)==19.,"unselected output column preserved");
        }
        compute_autocorrelation_brute_force(x,0,serial,1,CorrelationWindow::AllOrigins);
        for(long i=0;i<4;++i) check(serial(i,1)==parallel(i,1),"repeat call must overwrite instead of accumulate");
        MATRIX old(2,1,123.);
        compute_autocorrelation_brute_force(x,0,old,0);
        check(old(0,0)==1. && old(1,0)==2.,"legacy fixed-window compatibility");
        compute_autocorrelation_brute_force(x,0,serial,1,CorrelationWindow::AllOrigins,true);
        check(serial(0,1)==1.25,"full-series centering oracle");
        rejects([&]{compute_autocorrelation_brute_force(x,1,serial,1);});
        rejects([&]{compute_autocorrelation_brute_force(x,0,serial,2);});
        rejects([&]{compute_autocorrelation_brute_force_OMP(x,0,0,serial,1);});
        rejects([&]{compute_autocorrelation_brute_force(x,0,x,0);});
        MATRIX short_x(3,1,0.);rejects([&]{compute_autocorrelation_brute_force(short_x,0,old,0);});
        x(0,0)=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{compute_autocorrelation_brute_force(x,0,old,0);});
        std::cout<<"PASS legacy autocorrelation formulas, overwrite, OMP and validation\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
