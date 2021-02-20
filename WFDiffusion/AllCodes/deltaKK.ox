
deltaKK(s,x,y,Li,Ls,indm)

{

decl a,b,c2,c3,u1s,u2s,u1,u2,K,sigma,tau,u,j,S,Sj1,Sj2,xbar,ybar;

u1s=indm*(x-Li)+(1-indm)*(Ls-x);
u2s=indm*(y-Li)+(1-indm)*(Ls-y);

a=min(x,y)-Li;
b=Ls-max(x,y);

u1=(x-y+a-b)/2;
u2=(y-x+a-b)/2;
K=(fabs(x-y)+a+b)/2;

u=ranu(1,1);

S=0;
c3=0;
j=1;
 sigma=tau=<0>;
 while(c3==0)
 {
  sigma=sigmabar(s,u1,u2,K,j)+sigmabar(s,-u1,-u2,K,j);
  tau=taubar(s,u1,u2,K,j)+taubar(s,-u1,-u2,K,j);
  Sj2=1-S-sigma;
  S+=sigma-tau;
  Sj1=Sj2+tau;
  if( u<(Sj2/(1-exp(-2*u1s*u2s/s))) ){c3=1;c2=1;}
  else if(u>(Sj1/(1-exp(-2*u1s*u2s/s)))){c3=1;c2=0;}
  else{j+=1;}
 }


return(c2);
  
}	   
  