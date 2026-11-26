/* Local hard-sphere Metropolis MC, experiment-owned implementation.
 * No potential, forces, velocity or physical integration time.
 * PCG32 is separately retained, unmodified Apache-2.0 upstream source.
 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "hs_mc.h"
#include "pcg_basic.h"

struct HS {
    int n,m,nc,brute;
    double L,step;
    double *x,*disp;
    int *head,*next,*prev,*cell,*neighbors;
    pcg32_random_t rng;
};

static double dist2(const HS *s,const double *a,const double *b) {
    double r2=0;
    for(int k=0;k<3;k++) {
        double d=a[k]-b[k];
        if(d>0.5*s->L)d-=s->L;
        if(d< -0.5*s->L)d+=s->L;
        r2+=d*d;
    }
    return r2;
}
static int cellid(const HS *s,const double *x) {
    int c[3];
    for(int k=0;k<3;k++) {
        c[k]=(int)(x[k]/s->L*s->m);
        if(c[k]>=s->m)c[k]=s->m-1;
    }
    return (c[0]*s->m+c[1])*s->m+c[2];
}
static void insert(HS *s,int i,int cell) {
    s->cell[i]=cell;s->prev[i]=-1;s->next[i]=s->head[cell];
    if(s->head[cell]>=0)s->prev[s->head[cell]]=i;
    s->head[cell]=i;
}
static void remove_atom(HS *s,int i) {
    if(s->prev[i]>=0)s->next[s->prev[i]]=s->next[i];
    else s->head[s->cell[i]]=s->next[i];
    if(s->next[i]>=0)s->prev[s->next[i]]=s->prev[i];
}
int hs_allowed(const HS *s,int i,const double *trial) {
    if(!s || i<0 || i>=s->n)return 0;
    for(int k=0;k<3;k++)if(!isfinite(trial[k]) || trial[k]<0 || trial[k]>=s->L)return 0;
    if(s->brute) {
        for(int j=0;j<s->n;j++)if(j!=i && dist2(s,trial,s->x+3*j)<1.0)return 0;
    } else {
        int c=cellid(s,trial);
        for(int k=0;k<27;k++)for(int j=s->head[s->neighbors[27*c+k]];j>=0;j=s->next[j])
            if(j!=i && dist2(s,trial,s->x+3*j)<1.0)return 0;
    }
    return 1;
}
void hs_destroy(HS *s) {
    if(!s)return;
    free(s->x);free(s->disp);free(s->head);free(s->next);free(s->prev);
    free(s->cell);free(s->neighbors);free(s);
}
HS *hs_create(int n,double L,uint64_t seed,uint64_t stream,double step,int brute,const double *x) {
    if(!x || n<1 || n>4096 || !isfinite(L) || L<3.3 || L>100 || !isfinite(step) || step<=0 || step>=L)return NULL;
    HS *s=calloc(1,sizeof(HS));if(!s)return NULL;
    s->n=n;s->L=L;s->step=step;s->brute=brute;
    /* Cell width >=1.05: every collision and measured shell pair is in 27 cells. */
    s->m=(int)floor(L/1.05);s->nc=s->m*s->m*s->m;
    s->x=malloc(3*n*sizeof(double));s->disp=calloc(3*n,sizeof(double));
    s->head=malloc(s->nc*sizeof(int));s->next=malloc(n*sizeof(int));
    s->prev=malloc(n*sizeof(int));s->cell=malloc(n*sizeof(int));
    s->neighbors=malloc(27*s->nc*sizeof(int));
    if(!s->x || !s->disp || !s->head || !s->next || !s->prev || !s->cell || !s->neighbors){hs_destroy(s);return NULL;}
    memcpy(s->x,x,3*n*sizeof(double));
    for(int i=0;i<n;i++)for(int k=0;k<3;k++)if(!isfinite(x[3*i+k]) || x[3*i+k]<0 || x[3*i+k]>=L){hs_destroy(s);return NULL;}
    for(int c=0;c<s->nc;c++)s->head[c]=-1;
    for(int i=0;i<n;i++)insert(s,i,cellid(s,x+3*i));
    for(int a=0;a<s->m;a++)for(int b=0;b<s->m;b++)for(int c=0;c<s->m;c++) {
        int index=(a*s->m+b)*s->m+c,t=0;
        for(int da=-1;da<=1;da++)for(int db=-1;db<=1;db++)for(int dc=-1;dc<=1;dc++)
            s->neighbors[27*index+t++]=(((a+da+s->m)%s->m)*s->m+(b+db+s->m)%s->m)*s->m+(c+dc+s->m)%s->m;
    }
    /* Independent of cell acceleration: reject every initially overlapping pair. */
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(dist2(s,x+3*i,x+3*j)<1.0){hs_destroy(s);return NULL;}
    pcg32_srandom_r(&s->rng,seed,stream);return s;
}
uint64_t hs_sweeps(HS *s,int sweeps) {
    uint64_t accepted=0;
    for(int sweep=0;sweep<sweeps;sweep++)for(int attempt=0;attempt<s->n;attempt++) {
        int i=(int)pcg32_boundedrand_r(&s->rng,(uint32_t)s->n);
        double trial[3],increment[3];
        for(int k=0;k<3;k++) {
            double u=((double)pcg32_random_r(&s->rng)+0.5)*0x1p-32;
            increment[k]=s->step*(2*u-1);
            trial[k]=s->x[3*i+k]+increment[k];
            if(trial[k]<0)trial[k]+=s->L;
            if(trial[k]>=s->L)trial[k]-=s->L;
        }
        if(hs_allowed(s,i,trial)) {
            int c=cellid(s,trial);
            if(c!=s->cell[i]){remove_atom(s,i);insert(s,i,c);}
            for(int k=0;k<3;k++){s->x[3*i+k]=trial[k];s->disp[3*i+k]+=increment[k];}
            accepted++;
        }
    }
    return accepted;
}
void hs_copy(const HS *s,double *x,double *disp) {
    memcpy(x,s->x,3*s->n*sizeof(double));memcpy(disp,s->disp,3*s->n*sizeof(double));
}
void hs_measure(const HS *s,double *out) {
    const double widths[5]={.0025,.005,.01,.02,.04};
    double cuts[5];for(int k=0;k<5;k++)cuts[k]=(1+widths[k])*(1+widths[k]);
    for(int k=0;k<7;k++)out[k]=0;
    out[0]=1.05*1.05;
    for(int i=0;i<s->n;i++) {
        int c=s->cell[i];
        for(int k=0;k<27;k++)for(int j=s->head[s->neighbors[27*c+k]];j>=0;j=s->next[j])if(j>i) {
            double r2=dist2(s,s->x+3*i,s->x+3*j);
            if(r2<out[0])out[0]=r2;
            for(int t=0;t<5;t++)if(r2<cuts[t])out[1+t]+=1;
        }
        for(int k=0;k<3;k++)out[6]+=s->disp[3*i+k]*s->disp[3*i+k]/s->n;
    }
}
int hs_check_lists(const HS *s) {
    int *seen=calloc(s->n,sizeof(int));if(!seen)return 0;
    int good=1,count=0;
    for(int c=0;c<s->nc;c++) {
        int previous=-1,steps=0;
        for(int i=s->head[c];i>=0;i=s->next[i]) {
            if(i>=s->n || steps++>=s->n){good=0;break;}
            if(s->cell[i]!=c || s->prev[i]!=previous || cellid(s,s->x+3*i)!=c)good=0;
            seen[i]++;count++;previous=i;
        }
    }
    if(count!=s->n)good=0;
    for(int i=0;i<s->n;i++)if(seen[i]!=1)good=0;
    free(seen);return good;
}
