
deltaK1(t0,x0,T,xT,m,L)

{

decl u1,u2,u,c1,K,i,j,S,Sj1,Sj2,DI,Kh,zeta,xi,s,Kg,u1g,u2g,dn,sigma,tau;


s=T-t0;

K=fabs(L-m);

u1=fabs(x0-m);
u2=fabs(xT-m);

Kg=K/2;

u1g=u1-Kg;
u2g=u2-Kg;

dn=1-exp(-2*u1*u2/s);

u=ranu(1,1);

c1=0;

j=1;
sigma=tau=zeros(1,1);
S=0;
while(c1==0)
{
 sigma=sigmabar(s,u1g,u2g,Kg,j)+sigmabar(s,-u1g,-u2g,Kg,j);
 tau=taubar(s,u1g,u2g,Kg,j)+taubar(s,-u1g,-u2g,Kg,j);

  Sj2=1-S-sigma;
  S+=sigma-tau;
  Sj1=Sj2+tau;

 if(u<(Sj2/dn)){DI=1;c1=1;}
 else if(u>(Sj1/dn)){DI=0;c1=1;}
 else{j+=1;}
}

return(DI);
  
}	   
  