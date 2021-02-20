
DI(s,x,y,Li,Ls)

{

decl a,c2,c3,u1,u2,K,sigma,tau,u,j,S,Sj1,Sj2,xbar,ybar;


xbar=min(x,y);
ybar=max(x,y);

u1=(x-y)/2+((xbar-Li)-(Ls-ybar))/2;
u2=(y-x)/2+((xbar-Li)-(Ls-ybar))/2;
K=fabs(y-x)/2+((xbar-Li)+(Ls-ybar))/2;

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
  if(u<Sj2){c3=1;c2=1;}
  else if(u>Sj1){c3=1;c2=0;}
  else{j+=1;}
 }


return(c2);
  
}	   
  