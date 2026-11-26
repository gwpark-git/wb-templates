#include "../lib/hi_sd.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <string>
#include <gsl/gsl_eigen.h>

typedef std::vector<double> D;
void require(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
double error(const D&a,const D&b){require(a.size()==b.size(),"test shape");double e=0;for(size_t k=0;k<a.size();++k)e=std::max(e,std::abs(a[k]-b[k]));return e;}
double inner(const D&a,const D&b){double v=0;for(size_t k=0;k<a.size();++k)v+=a[k]*b[k];return v;}
D rotate(const D& x,const D&r,double factor=1){D y(x.size(),0.);for(size_t i=0;i<x.size()/3;++i)for(int a=0;a<3;++a)for(int b=0;b<3;++b)y[3*i+a]+=factor*r[3*a+b]*x[3*i+b];return y;}
D rotate_tensor(const D& x,const D&r){D y(x.size(),0.);for(size_t i=0;i<x.size()/9;++i)for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int c=0;c<3;++c)for(int d=0;d<3;++d)y[9*i+3*a+b]+=r[3*a+c]*x[9*i+3*c+d]*r[3*b+d];return y;}
D reordered(const D& x,size_t stride,const std::vector<size_t>&perm){D y(x.size());for(size_t i=0;i<perm.size();++i)for(size_t k=0;k<stride;++k)y[stride*i+k]=x[stride*perm[i]+k];return y;}
int main(){try{
 const double pi=acos(-1.);int groups=0;double worst_residual=0.;
 HI_SD::Parameters p;p.radius=.7;p.viscosity=1.3;
 D q{0,0,0},F{.1,-.2,.3},T{.4,.2,-.1},zero(3,0.),E{.03,.02,-.01,.02,-.01,.04,-.01,.04,-.02};
 auto one=HI_SD::solve(q,F,T,zero,zero,E,p);D eu=F,ew=T,es=E;
 for(double&v:eu)v/=6*pi*p.viscosity*p.radius;
 for(double&v:ew)v/=8*pi*p.viscosity*pow(p.radius,3);
 for(double&v:es)v*=20*pi*p.viscosity*pow(p.radius,3)/3;
 require(error(one.translation,eu)<1e-13,"isolated Stokes drag");require(error(one.rotation,ew)<1e-13,"isolated rotational drag");require(error(one.stresslet,es)<1e-13,"isolated Einstein stresslet");++groups;
 // A compact non-collinear near-field configuration exercises every tensor mode.
 p.radius=1;p.viscosity=1./(6*pi);
 q=D{0,0,0,2.1,.4,.1,-.4,2.7,.6};F=D{.1,-.2,.3,-.2,.1,.03,.1,.1,-.2};T=D{.04,.02,-.01,-.02,.03,.01,.02,-.01,.04};zero=D(9,0.);
 D ee=E;ee.insert(ee.end(),E.begin(),E.end());ee.insert(ee.end(),E.begin(),E.end());
 D bg{.01,.02,-.03,.02,.01,.04,-.01,.03,.02},om{.04,.03,.02,.01,.02,.03,.03,.02,.01};
 auto base=HI_SD::solve(q,F,T,bg,om,ee,p);worst_residual=std::max(worst_residual,base.relative_residual);
 // Proper 90-degree permutation rotation, plus a non-axis-aligned orthogonal rotation.
 D R{.36,-.48,.8,.8,.6,0.,-.48,.64,.6};
 auto rot=HI_SD::solve(rotate(q,R),rotate(F,R),rotate(T,R),rotate(bg,R),rotate(om,R),rotate_tensor(ee,R),p);
 double covariance_error=std::max({error(rot.translation,rotate(base.translation,R)),error(rot.rotation,rotate(base.rotation,R)),error(rot.stresslet,rotate_tensor(base.stresslet,R))});
 require(covariance_error<2e-11,"proper rotation covariance");++groups;
 D reflection{-1,0,0,0,1,0,0,0,1};
 auto reflected=HI_SD::solve(rotate(q,reflection),rotate(F,reflection),rotate(T,reflection,-1),rotate(bg,reflection),rotate(om,reflection,-1),rotate_tensor(ee,reflection),p);
 double reflection_error=std::max({error(reflected.translation,rotate(base.translation,reflection)),error(reflected.rotation,rotate(base.rotation,reflection,-1)),error(reflected.stresslet,rotate_tensor(base.stresslet,reflection))});
 require(reflection_error<2e-11,"improper reflection with axial torque/rotation");++groups;
 std::vector<size_t> perm{2,0,1};
 auto permuted=HI_SD::solve(reordered(q,3,perm),reordered(F,3,perm),reordered(T,3,perm),reordered(bg,3,perm),reordered(om,3,perm),reordered(ee,9,perm),p);
 double permutation_error=std::max({error(permuted.translation,reordered(base.translation,3,perm)),error(permuted.rotation,reordered(base.rotation,3,perm)),error(permuted.stresslet,reordered(base.stresslet,9,perm))});
 require(permutation_error<2e-11,"particle label permutation invariance");++groups;
 // Lorentz reciprocal work: the strain vector uses full S:E tensor contraction.
 D e0(9,0.);auto force=HI_SD::solve(q,F,T,zero,zero,e0,p);auto strain=HI_SD::solve(q,zero,zero,zero,zero,ee,p);
 double reciprocal_work=inner(F,strain.translation)+inner(T,strain.rotation)+inner(force.stresslet,ee);
 require(std::abs(reciprocal_work)<2e-11,"force/torque--strain reciprocal work");++groups;
 // Effective translation/rotation mobility, with zero strain, must be symmetric positive.
 const size_t m=18;D mobility(m*m),column(m),ff(9),tt(9);
 for(size_t k=0;k<m;++k){std::fill(ff.begin(),ff.end(),0);std::fill(tt.begin(),tt.end(),0);if(k%6<3)ff[3*(k/6)+k%6]=1;else tt[3*(k/6)+k%6-3]=1;
  auto x=HI_SD::solve(q,ff,tt,zero,zero,e0,p);for(size_t i=0;i<3;++i)for(int d=0;d<3;++d){mobility[(6*i+d)*m+k]=x.translation[3*i+d];mobility[(6*i+3+d)*m+k]=x.rotation[3*i+d];}}
 double symmetry_error=0;for(size_t i=0;i<m;++i)for(size_t j=0;j<m;++j)symmetry_error=std::max(symmetry_error,std::abs(mobility[i*m+j]-mobility[j*m+i]));
 require(symmetry_error<2e-11,"Lorentz mobility symmetry");
 gsl_matrix_view mm=gsl_matrix_view_array(mobility.data(),m,m);gsl_vector* eval=gsl_vector_alloc(m);gsl_eigen_symm_workspace*w=gsl_eigen_symm_alloc(m);gsl_eigen_symm(&mm.matrix,eval,w);double mineig=gsl_vector_get(eval,0);for(size_t i=1;i<m;++i)mineig=std::min(mineig,gsl_vector_get(eval,i));gsl_eigen_symm_free(w);gsl_vector_free(eval);require(mineig>0,"positive effective mobility");++groups;
 // Squeeze mode approaches U1/(2h F)=1 as h/a ->0 with a=6*pi*eta=1.
 D pairF{1,0,0,-1,0,0},z6(6,0.);double squeeze_ratio=0.;
 for(double h:D{.001,.0002}){auto x=HI_SD::solve(D{0,0,0,2+h,0,0},pairF,z6,z6,z6,e0,p);squeeze_ratio=x.translation[0]/(2*h);require(std::abs(squeeze_ratio-1)<.025,"normal squeeze lubrication asymptotic");require(std::abs(x.translation[0]+x.translation[3])<1e-12,"opposite squeezing velocities");}++groups;
 // A perfectly symmetric triplet must not acquire an origin drift or unequal outer response.
 D triplet{0,-2.1,0,0,0,0,0,2.1,0},z9(9,0.),shear{0,.5,0,.5,0,0,0,0,0};
 auto tr=HI_SD::solve(triplet,z9,z9,z9,z9,shear,p);require(std::abs(tr.translation[0]+tr.translation[6])<1e-12&&std::abs(tr.translation[3])<1e-12,"triplet reflection velocity symmetry");require(std::abs(tr.rotation[2]-tr.rotation[8])<1e-12,"triplet outer rotation symmetry");require(std::abs(tr.stresslet[1]-tr.stresslet[19])<1e-12,"triplet outer stresslet symmetry");++groups;
 // Literal values from the separate unmodified JFSD source-operator probe,
 // third_party/hi_sd/reference_cases.json; parity corrections vanish in this case.
 auto axial=HI_SD::solve(D{0,0,0,2.1,0,0},pairF,z6,z6,z6,e0,p);
 require(std::abs(axial.translation[0]-.1349660045323884)<2e-13,"matched JFSD axial-pair velocity");
 require(std::abs(axial.stresslet[0]+.5126561534102703)<2e-13&&std::abs(axial.stresslet[4]-.25632807670513513)<2e-13,"matched JFSD axial-pair stresslet");++groups;
 for(size_t i=0;i<3;++i){require(std::abs(base.stresslet[9*i]+base.stresslet[9*i+4]+base.stresslet[9*i+8])<1e-12,"traceless stresslet");for(int a=0;a<3;++a)for(int b=0;b<3;++b)require(std::abs(base.stresslet[9*i+3*a+b]-base.stresslet[9*i+3*b+a])<1e-12,"symmetric stresslet");}++groups;
 int rejections=0;
 try{HI_SD::solve(D{0,0,0,1.99,0,0},pairF,z6,z6,z6,e0,p);}catch(const std::domain_error&){++rejections;}
 try{HI_SD::solve(D{0,0,0,2.00001,0,0},pairF,z6,z6,z6,e0,p);}catch(const std::domain_error&){++rejections;}
 try{D bad=E;bad[1]+=.1;HI_SD::solve(D{0,0,0},D(3),D(3),D(3),D(3),bad,p);}catch(const std::invalid_argument&){++rejections;}
 try{D bad=E;bad[0]+=.1;HI_SD::solve(D{0,0,0},D(3),D(3),D(3),D(3),bad,p);}catch(const std::invalid_argument&){++rejections;}
 try{HI_SD::Parameters bad=p;bad.viscosity=-1;HI_SD::solve(D{0,0,0},D(3),D(3),D(3),D(3),E,bad);}catch(const std::invalid_argument&){++rejections;}
 try{HI_SD::solve(D{0,0,0},D(2),D(3),D(3),D(3),E,p);}catch(const std::invalid_argument&){++rejections;}
 require(rejections==6,"invalid-input/overlap/unresolved-gap rejection");++groups;
 std::cout<<std::setprecision(17)<<"{\"status\":\"passed\",\"check_groups\":"<<groups<<",\"rotation_covariance_error\":"<<covariance_error<<",\"reflection_error\":"<<reflection_error<<",\"permutation_error\":"<<permutation_error<<",\"reciprocal_work_error\":"<<std::abs(reciprocal_work)<<",\"mobility_symmetry_error\":"<<symmetry_error<<",\"minimum_effective_mobility_eigenvalue\":"<<mineig<<",\"squeeze_ratio_at_gap_0_0002\":"<<squeeze_ratio<<",\"axial_pair_translation\":"<<axial.translation[0]<<",\"symmetric_triplet_outer_Ux\":"<<tr.translation[0]<<",\"symmetric_triplet_outer_Omegaz\":"<<tr.rotation[2]<<",\"symmetric_triplet_outer_Sxy\":"<<tr.stresslet[1]<<",\"worst_saddle_residual\":"<<worst_residual<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"SD test failure: "<<e.what()<<"\n";return 1;}}
