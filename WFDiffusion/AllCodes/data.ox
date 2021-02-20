
#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>


main()
{
  decl time,th1,th2,dtt,dt,x0,T,r,n,nn,Xt,xt,i,nT,j,g1,g2,sigma;

time=timer();
//ranseed(1);
ranseed(today());


// Parameters

th1=4;
th2=4;
g1=th1+th2;
g2=th1/(th1+th2);
sigma=1;



x0=0.5;           // Starting value
T=0.1;
nT=2000;		  // ntxT = final observed time
dt=1;             // time interval between observations
dtt=0.000000001;  // time discretisation for the Euler scheme

r=dt/dtt;

n=(T*nT)/dt;
Xt=zeros(n+1,1);
Xt[0]=x0;

xt=x0;

nn=round(T/dtt);

for(j=1;j<=nT;++j)
{

for(i=1;i<=nn;++i)
{										   
xt=xt + ( -(g1/2)*(xt-g2))*dtt + sigma*sqrt(dtt*xt*(1-xt))*rann(1,1); if(fabs(xt-0.5)>=0.5){println(xt);} // if xt is printed, it means that xt went outside the interval (0,1)

}

if(fmod(j+10,10)==0){Xt[round(j/10)]=xt;println(round(j/10)~xt);}

}


println("Time = ",timespan(time));



savemat("X.mat",2*asin(sqrt(Xt))); // saves the Lamperti transform of the simulated WF diffusion


}