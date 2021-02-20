#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>
#include "rnorm.ox"
#include "dnorm.ox"
#include "rmvnorm.ox"
#include "runif.ox"
#include "BB2.ox"
#include "ivec.ox"

#include "sigmabar.ox"
#include "taubar.ox"
#include "DeltaK1.ox"
#include "DeltaK.ox"
#include "DeltaL.ox"
#include "DeltaKmin.ox"
#include "DeltaKmax.ox"
#include "igaus.ox"
#include "condminim.ox"
#include "condmaxim.ox"
#include "xminim.ox"
#include "xmaxim.ox"
#include "xminim2.ox"
#include "xmaxim2.ox"
#include "layerass.ox"
#include "xlayer.ox"
#include "xlayerass0.ox"

#include "DI.ox"

#include "phi.ox"
#include "dphi.ox"
#include "ddphi.ox"
#include "nrphi.ox"
#include "phi2.ox"
#include "dphi2.ox"
#include "ddphi2.ox"
#include "nrphi2.ox"
#include "Inti0.ox"
#include "xseed.ox"

#include "DeltaKK.ox"
#include "BBx.ox"
#include "XBarkerPAR_new.ox"

#include "Intii20.ox"
#include "GBarker2_new.ox"


main()
{
  decl time,Xt,th1,th2,g1,g2,T,dt,dtt,n,X0,M,gamaa,X1,X2,l,i,gama,s1,s2,r,cp,lims,MINI,MAXI,x0,xT,l_u,M1;

time=timer();
ranseed(1);
//ranseed(today());

decl conta,contb,cont1,cont2a,cont2b,cont2,temp,temp1;

conta=contb=cont1=0;

Xt=loadmat("X.mat");			// load data set


x0=Xt[0];
xT=Xt[rows(Xt)-1];

MINI=0.2;					// restricts the original WF diffusion to the interval (0.0099 , 0.9898) to avoid numerical problems - the observed data lies in (0.189 , 0.920) (or (0.900 , 2.570) under the Lamperti transform)
MAXI=2.94;

// Parameters

th1=4;
th2=4;
g1=th1+th2;
g2=th1/(th1+th2);

// Variables

T=50;			  // Final time
dt=1;			  // time length between consecutive observations
dtt=0.02;		  // interval refinement to apply the layered Brownian bridge algorithm

n=rows(Xt);

X0=Xt[:(n-2)]~Xt[1:];


M=10000;		  // number of MCMC iteration

l_u=10;			  // start updating the paramenters after l_u iterations of the chain

M1=5;			  // number of consecutive updates of gamma1 and of gamma2

gamaa=zeros(M+1,2);
gamaa[0][0]=g1;
gamaa[0][1]=g2;

cont2a=zeros(5*M,1);
cont2b=zeros(5*M,1);
cont2=zeros(M,1);

s1=0.3;		 // proposal for gamma1 -> U(gamma-s1 , gamma+s1)
s2=0.01;	 // proposal for gamma2 -> U(gamma-s2 , gamma+s2)


cp=nrphi(g1,g2,M_PI/2,10^(-10));


[X1,X2]=xseed(X0,dt,dtt,MINI,MAXI,cp,g1,g2);


decl beta1=0.98;  //0.95	Beta parameter for the Portkey Barker's of X
decl beta2=0.99; //0.98		Beta parameter for the Portkey Barker's of theta


for(l=1;l<=M;++l)
{


[X1,X2,temp,temp1]=XBarkerPAR(X0,X1,X2,T,dt,dtt,g1,g2,cp,MINI,MAXI,beta1);
cont1+=temp;
//cont2[l-1]=temp1;

if(l>l_u){
for(i=1;i<=M1;++i)
{
[X2,gama,cp,temp,temp1]=GBarker2(X1, X2, dtt, g1, g2, cp, s1, 0, MINI,MAXI,x0,xT,beta2);
if(g1!=gama[0]){conta+=1;}
//cont2a[M1*(l-1)+i-1]=temp;
g1=gama[0];
}
for(i=1;i<=M1;++i)
{
[X2,gama,cp,temp,temp1]=GBarker2(X1, X2, dtt, g1, g2, cp, 0, s2, MINI,MAXI,x0,xT,beta2);
if(g2!=gama[1]){contb+=1;}
//cont2b[M1*(l-1)+i-1]=temp;
g2=gama[1];
}
}

gamaa[l][]=g1~g2;
println(l~g1~g2);



}



println((conta/(M1*(M-l_u)))~(contb/(M1*(M-l_u))));	 // accetance rate of gamma1 and gamma2

println(cont1/(M));		 // average accetance rate of X

println("Time = ",timespan(time));
println("");



savemat("P:\\Statistics\\UFMG\\Pesquisa\\Dootika\\code3\\gama_new.mat",gamaa[10:][]);
savemat("P:\\Statistics\\UFMG\\Pesquisa\\Dootika\\code3\\cont2a_new.mat",cont2a);
savemat("P:\\Statistics\\UFMG\\Pesquisa\\Dootika\\code3\\cont2b_new.mat",cont2b);
savemat("P:\\Statistics\\UFMG\\Pesquisa\\Dootika\\code3\\cont2_new.mat",cont2);



}