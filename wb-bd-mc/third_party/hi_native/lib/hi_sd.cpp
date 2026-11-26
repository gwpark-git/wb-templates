// SPDX-License-Identifier: Apache-2.0
// Selected formulas derived from JFSD, William Torre and contributors,
// commit5234c959c0fdb0ba7ba32d6015d0478014a512c7. See third_party/hi_sd/NOTICE.
// Native invariant-tensor/dense-GSL implementation; no JAX runtime.
#include "hi_sd.h"
#include "../third_party/hi_sd/resistance_table.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace HI_SD {
namespace {
typedef std::array<double,3> V;
typedef std::array<double,9> T;
typedef std::vector<double> D;
V zero_v() { return V{{0.,0.,0.}}; }
T zero_t() { return T{{0.,0.,0.,0.,0.,0.,0.,0.,0.}}; }
double dot(const V& a,const V& b) { double x=0;for(int k=0;k<3;++k)x+=a[k]*b[k];return x; }
V scale(V a,double b) { for(int k=0;k<3;++k)a[k]*=b;return a; }
V add(V a,const V& b,double c=1.) { for(int k=0;k<3;++k)a[k]+=c*b[k];return a; }
T add(T a,const T& b,double c=1.) { for(int k=0;k<9;++k)a[k]+=c*b[k];return a; }
T scale(T a,double b) { for(int k=0;k<9;++k)a[k]*=b;return a; }
V cross(const V& a,const V& b) { return V{{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}}; }
V matvec(const T& a,const V& b) { V c=zero_v();for(int k=0;k<3;++k)for(int l=0;l<3;++l)c[k]+=a[3*k+l]*b[l];return c; }
V getv(const D& x,size_t i,size_t stride=3,size_t offset=0) { return V{{x[stride*i+offset],x[stride*i+offset+1],x[stride*i+offset+2]}}; }
void addv(D& x,size_t i,const V& v,size_t stride=3,size_t offset=0,double s=1.) { for(int k=0;k<3;++k)x[stride*i+offset+k]+=s*v[k]; }
T stress_from5(const D& x,size_t offset) {
  return T{{x[offset],x[offset+1],x[offset+2],x[offset+1],x[offset+4],x[offset+3],x[offset+2],x[offset+3],-x[offset]-x[offset+4]}};
}
void add_dual(D& x,size_t offset,const T& e) {
  x[offset]+=2*e[0]+e[4];x[offset+1]+=2*e[1];x[offset+2]+=2*e[2];
  x[offset+3]+=2*e[5];x[offset+4]+=e[0]+2*e[4];
}
V axial(const V& n,const V& v,double xa,double ya) { return add(scale(v,ya),n,(xa-ya)*dot(n,v)); }
T G(const V& n,const V& v,double xg,double yg) {
  const double nv=dot(n,v);T out=zero_t();
  for(int a=0;a<3;++a)for(int b=0;b<3;++b)out[3*a+b]=xg*(n[a]*n[b]-(a==b?1./3.:0.))*nv+yg*(v[a]*n[b]+n[a]*v[b]-2*n[a]*n[b]*nv);
  return out;
}
T H(const V& n,const V& v,double yh) {
  const V e=cross(v,n);T out=zero_t();
  for(int a=0;a<3;++a)for(int b=0;b<3;++b)out[3*a+b]=yh*(n[a]*e[b]+e[a]*n[b]);
  return out;
}
T M(const V& n,const T& s,double xm,double ym,double zm) {
  const V sn=matvec(s,n);const double nsn=dot(n,sn);T out=zero_t();
  for(int a=0;a<3;++a)for(int b=0;b<3;++b){
    const double nn=n[a]*n[b],d=(a==b?1.:0.);
    out[3*a+b]=1.5*xm*(nn-d/3.)*nsn+ym*(n[a]*sn[b]+sn[a]*n[b]-2*nn*nsn)
      +.5*zm*(2*s[3*a+b]+(d+nn)*nsn-2*n[a]*sn[b]-2*sn[a]*n[b]);
  }return out;
}
V Gcontract(const V& n,const T& e,double xg,double yg) {
  const V en=matvec(e,n);return add(scale(en,2*yg),n,(xg-2*yg)*dot(n,en));
}
struct Pair { size_t i,j;V n;double r;std::array<double,22> c;bool lub; };

std::array<double,22> correction(double r) {
  // Distances themselves are tabulated; binary search avoids rounding a log
  // index onto the wrong interval. No hidden contact-floor regularization.
  if(r<HI_SD_TABLE::data[0][0])throw std::domain_error("SD lubrication gap below tabulated minimum 1e-4 a");
  size_t lo=0,hi=HI_SD_TABLE::count-1;
  while(hi-lo>1){const size_t mid=(lo+hi)/2;if(HI_SD_TABLE::data[mid][0]<=r)lo=mid;else hi=mid;}
  const double f=(r-HI_SD_TABLE::data[lo][0])/(HI_SD_TABLE::data[hi][0]-HI_SD_TABLE::data[lo][0]);
  std::array<double,22> c;
  for(size_t k=0;k<22;++k)c[k]=HI_SD_TABLE::data[lo][k+1]+f*(HI_SD_TABLE::data[hi][k+1]-HI_SD_TABLE::data[lo][k+1]);
  return c;
}

D grand_apply(size_t n,const std::vector<Pair>& pairs,const D& x) {
  D out(11*n,0.);std::vector<T> E(n,zero_t());
  // Self terms are explicit, including N=1 (upstream pair accumulation does
  // not supply an isolated-sphere self term in its open implementation).
  for(size_t i=0;i<n;++i){for(int a=0;a<3;++a){out[6*i+a]=x[6*i+a];out[6*i+3+a]=.75*x[6*i+3+a];}
    E[i]=add(E[i],stress_from5(x,6*n+5*i),.9);}
  for(const Pair& p:pairs){
    const double r=p.r,r2=r*r,r3=r*r2,r4=r2*r2,r5=r2*r3;
    const double xa=1.5/r-1/r3,ya=.75/r+.5/r3,yb=-.75/r2,xc=.75/r3,yc=-.375/r3;
    const double xm=-4.5/r3+10.8/r5,ym=2.25/r3-7.2/r5,zm=1.8/r5,xg=2.25/r2-3.6/r4,yg=1.2/r4,yh=-1.125/r3;
    const V fi=scale(getv(x,p.i,6),-1),fj=scale(getv(x,p.j,6),-1);
    const V ti=scale(getv(x,p.i,6,3),-1),tj=scale(getv(x,p.j,6,3),-1);
    const T si=stress_from5(x,6*n+5*p.i),sj=stress_from5(x,6*n+5*p.j);
    const V sin=matvec(si,p.n),sjn=matvec(sj,p.n);
    V ui=add(axial(p.n,fj,xa,ya),cross(tj,p.n),yb);
    V uj=add(axial(p.n,fi,xa,ya),cross(ti,p.n),-yb);
    ui=add(add(ui,p.n,(-xg+2*yg)*dot(p.n,sjn)),sjn,-2*yg);
    uj=add(add(uj,p.n,(xg-2*yg)*dot(p.n,sin)),sin,2*yg);
    V wi=add(add(axial(p.n,tj,xc,yc),cross(fj,p.n),yb),cross(sjn,p.n),2*yh);
    // Omega--S and E--T are even under reversing the line of centres.
    // The pinned upstream j-side signs differ; enforce tensor parity here.
    V wj=add(add(axial(p.n,ti,xc,yc),cross(fi,p.n),-yb),cross(sin,p.n),2*yh);
    addv(out,p.i,ui,6,0,-1);addv(out,p.j,uj,6,0,-1);
    addv(out,p.i,wi,6,3,-1);addv(out,p.j,wj,6,3,-1);
    E[p.i]=add(add(add(E[p.i],G(p.n,fj,xg,yg)),H(p.n,tj,yh),-1),M(p.n,sj,xm,ym,zm));
    E[p.j]=add(add(add(E[p.j],G(p.n,fi,xg,yg),-1),H(p.n,ti,yh),-1),M(p.n,si,xm,ym,zm));
  }
  for(size_t i=0;i<n;++i)add_dual(out,6*n+5*i,E[i]);
  return out;
}

D lub_apply(size_t n,const std::vector<Pair>& pairs,const D& u) {
  D f(6*n,0.);
  for(const Pair& p:pairs)if(p.lub){const auto& c=p.c;
    V ui=getv(u,p.i,6),uj=getv(u,p.j,6),wi=getv(u,p.i,6,3),wj=getv(u,p.j,6,3);
    V fi=add(axial(p.n,ui,c[0],c[2]),axial(p.n,uj,c[1],c[3]));
    fi=add(add(fi,cross(wi,p.n),-c[4]),cross(wj,p.n),c[5]);
    V fj=add(axial(p.n,uj,c[0],c[2]),axial(p.n,ui,c[1],c[3]));
    fj=add(add(fj,cross(wj,p.n),c[4]),cross(wi,p.n),-c[5]);
    V ti=add(add(axial(p.n,wi,c[6],c[8]),axial(p.n,wj,c[7],c[9])),cross(ui,p.n),c[4]);
    ti=add(ti,cross(uj,p.n),c[5]);
    V tj=add(add(axial(p.n,wj,c[6],c[8]),axial(p.n,wi,c[7],c[9])),cross(uj,p.n),-c[4]);
    tj=add(tj,cross(ui,p.n),-c[5]);
    addv(f,p.i,fi,6);addv(f,p.j,fj,6);addv(f,p.i,ti,6,3);addv(f,p.j,tj,6,3);
  }return f;
}
D lub_strain_rhs(size_t n,const std::vector<Pair>& pairs,const std::vector<T>& e) {
  D f(6*n,0.);
  for(const Pair& p:pairs)if(p.lub){const auto& c=p.c;
    V fi=add(scale(Gcontract(p.n,e[p.i],c[10],c[12]),-1),Gcontract(p.n,e[p.j],c[11],c[13]));
    V fj=add(Gcontract(p.n,e[p.j],c[10],c[12]),Gcontract(p.n,e[p.i],c[11],c[13]),-1);
    V ti=add(scale(cross(matvec(e[p.i],p.n),p.n),2*c[14]),cross(matvec(e[p.j],p.n),p.n),2*c[15]);
    V tj=add(scale(cross(matvec(e[p.j],p.n),p.n),2*c[14]),cross(matvec(e[p.i],p.n),p.n),2*c[15]);
    addv(f,p.i,fi,6);addv(f,p.j,fj,6);addv(f,p.i,ti,6,3);addv(f,p.j,tj,6,3);
  }return f;
}
std::vector<T> lub_stress(size_t n,const std::vector<Pair>& pairs,const D& u,const std::vector<T>& e) {
  std::vector<T> out(n,zero_t());
  for(const Pair& p:pairs)if(p.lub){const auto& c=p.c;
    const V ui=getv(u,p.i,6),uj=getv(u,p.j,6),wi=getv(u,p.i,6,3),wj=getv(u,p.j,6,3);
    T si=add(scale(G(p.n,ui,c[10],c[12]),-1),G(p.n,uj,c[11],c[13]),-1);
    T sj=add(G(p.n,uj,c[10],c[12]),G(p.n,ui,c[11],c[13]));
    si=add(add(si,H(p.n,wi,c[14]),-1),H(p.n,wj,c[15]),-1);
    sj=add(add(sj,H(p.n,wj,c[14]),-1),H(p.n,wi,c[15]),-1);
    si=add(add(si,M(p.n,e[p.i],c[16],c[18],c[20])),M(p.n,e[p.j],c[17],c[19],c[21]));
    sj=add(add(sj,M(p.n,e[p.j],c[16],c[18],c[20])),M(p.n,e[p.i],c[17],c[19],c[21]));
    out[p.i]=add(out[p.i],si);out[p.j]=add(out[p.j],sj);
  }return out;
}
void finite(const D& x,const std::string& name) { for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("SD nonfinite "+name); }
}

Result solve(const D& positions,const D& forces,const D& torques,const D& background_velocity,
             const D& background_angular_velocity,const D& strain,const Parameters& p) {
  const size_t n=positions.size()/3;
  if(n==0||positions.size()!=3*n||n>64)throw std::invalid_argument("SD needs1..64 3D spheres (dense reference)");
  if(forces.size()!=3*n||torques.size()!=3*n||background_velocity.size()!=3*n||background_angular_velocity.size()!=3*n||(strain.size()!=9&&strain.size()!=9*n))throw std::invalid_argument("SD vector/strain shape mismatch");
  finite(positions,"positions");finite(forces,"forces");finite(torques,"torques");finite(background_velocity,"ambient velocity");finite(background_angular_velocity,"ambient rotation");finite(strain,"strain");
  if(!(p.radius>0)||!std::isfinite(p.radius)||!(p.viscosity>0)||!std::isfinite(p.viscosity)||!std::isfinite(p.lubrication_cutoff)||p.lubrication_cutoff<=2||p.lubrication_cutoff>3.99)throw std::invalid_argument("SD invalid radius/viscosity/lubrication cutoff");
  std::vector<T> e(n,zero_t());
  for(size_t i=0;i<n;++i){const size_t o=(strain.size()==9?0:9*i);double norm=0;
    for(int k=0;k<9;++k){e[i][k]=strain[o+k];norm=std::max(norm,std::abs(e[i][k]));}
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)if(std::abs(e[i][3*a+b]-e[i][3*b+a])>1e-12*(1+norm))throw std::invalid_argument("SD strain must be symmetric");
    if(std::abs(e[i][0]+e[i][4]+e[i][8])>1e-12*(1+norm))throw std::invalid_argument("SD strain must be traceless");
  }
  const double pi=std::acos(-1.),drag=6*pi*p.viscosity*p.radius;
  if(!std::isfinite(drag*p.radius*p.radius)||drag*p.radius*p.radius==0)throw std::invalid_argument("SD unrepresentable unit scales");
  std::vector<Pair> pairs;double mingap=std::numeric_limits<double>::infinity();
  for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j){Pair q;q.c.fill(0.);q.i=i;q.j=j;q.n=scale(add(getv(positions,j),getv(positions,i),-1),1/p.radius);q.r=std::sqrt(dot(q.n,q.n));
    if(!std::isfinite(q.r)||q.r<=2)throw std::domain_error("SD overlapping/contact spheres are unsupported");
    mingap=std::min(mingap,q.r-2);q.n=scale(q.n,1/q.r);q.lub=p.lubrication&&q.r<p.lubrication_cutoff;
    if(q.lub)q.c=correction(q.r);
    pairs.push_back(q);
  }
  const size_t z=17*n;D a(z*z,0.),rhs(z,0.),x(z,0.);
  for(size_t col=0;col<11*n;++col){D v(11*n,0.);v[col]=1.;D w=grand_apply(n,pairs,v);for(size_t row=0;row<11*n;++row)a[row*z+col]=w[row];}
  for(size_t col=0;col<6*n;++col){D v(6*n,0.);v[col]=1.;D w=lub_apply(n,pairs,v);
    a[col*z+11*n+col]=1.;a[(11*n+col)*z+col]=1.;
    for(size_t row=0;row<6*n;++row)a[(11*n+row)*z+11*n+col]=-w[row];
  }
  const D fe=lub_strain_rhs(n,pairs,e);
  for(size_t i=0;i<n;++i){T minus_e=e[i];for(double& v:minus_e)v=-v;add_dual(rhs,6*n+5*i,minus_e);
    for(int d=0;d<3;++d){rhs[11*n+6*i+d]=-forces[3*i+d]/(drag*p.radius)+fe[6*i+d];rhs[11*n+6*i+3+d]=-torques[3*i+d]/(drag*p.radius*p.radius)+fe[6*i+3+d];}
  }
  D fact=a;gsl_matrix_view matrix=gsl_matrix_view_array(fact.data(),z,z);gsl_vector_view b=gsl_vector_view_array(rhs.data(),z),sol=gsl_vector_view_array(x.data(),z);
  gsl_permutation* permutation=gsl_permutation_alloc(z);int sign=0;
  // LU returns status; inspect diagonal before the solve rather than accepting
  // a singular configuration or adding mobility jitter.
  int status=gsl_linalg_LU_decomp(&matrix.matrix,permutation,&sign);
  if(status==GSL_SUCCESS)for(size_t k=0;k<z;++k)if(!std::isfinite(fact[k*z+k])||std::abs(fact[k*z+k])<1e-14){status=GSL_ESING;break;}
  if(status==GSL_SUCCESS)status=gsl_linalg_LU_solve(&matrix.matrix,permutation,&b.vector,&sol.vector);
  gsl_permutation_free(permutation);
  if(status!=GSL_SUCCESS)throw std::runtime_error("SD dense saddle solve failed: "+std::string(gsl_strerror(status)));
  finite(x,"solution");double rnorm=0,bnorm=0;
  for(size_t i=0;i<z;++i){double r=-rhs[i];for(size_t j=0;j<z;++j)r+=a[i*z+j]*x[j];rnorm+=r*r;bnorm+=rhs[i]*rhs[i];}
  Result result;result.relative_residual=std::sqrt(rnorm)/std::max(1.,std::sqrt(bnorm));
  if(result.relative_residual>1e-9)throw std::runtime_error("SD saddle residual exceeds1e-9");
  result.minimum_gap_over_radius=(n>1?mingap:std::numeric_limits<double>::max());
  result.translation.resize(3*n);result.rotation.resize(3*n);result.stresslet.resize(9*n);
  D u(x.begin()+11*n,x.end());const std::vector<T> stress_nf=lub_stress(n,pairs,u,e);
  for(size_t i=0;i<n;++i){T s=add(stress_nf[i],stress_from5(x,6*n+5*i),-1);
    for(int d=0;d<3;++d){result.translation[3*i+d]=background_velocity[3*i+d]+p.radius*u[6*i+d];result.rotation[3*i+d]=background_angular_velocity[3*i+d]+u[6*i+3+d];}
    for(int k=0;k<9;++k)result.stresslet[9*i+k]=drag*p.radius*p.radius*s[k];
  }
  finite(result.translation,"translation result");finite(result.rotation,"rotation result");finite(result.stresslet,"stresslet result");return result;
}
} // namespace HI_SD
