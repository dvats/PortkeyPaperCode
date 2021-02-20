
deltaKmax(t0,x0,T,XT,s0,Xt,m,aI)

{

decl a,t,x,y,c1,c2,c3,n,u1,u2,K,i,sigma,tau,zeta,xi,u,j,Sj1,Sj2,a1,DI,xbar,ybar,sg,u1g,u2g,Kg,s,Kh;

x=x0;
y=XT;
t=T-t0;
xbar=min(x,y);
ybar=max(x,y);

n=rows(s0);

s=sg=u1=u1g=u2=u2g=zeros(n-1,1);

K=m-xbar+aI;

for(i=1;i<=(n-1);++i)
{
 s[i-1]=s0[i]-s0[i-1];
 sg[i-1]=s[i-1];
 u1[i-1]=m-Xt[i-1];
 u1g[i-1]=u1[i-1]-K/2;
 u2[i-1]=m-Xt[i]; 
 u2g[i-1]=u2[i-1]-K/2;
}							  // println((K-u1)~(K-u2));

Kg=K/2;

u=ranu(1,1);
//println(u);

c1=0;

Kh=ceil(sqrt(t+K^2)/(2*K));


if(3*K^2>t || Kh==1)
{

j=1;
sigma=tau=zeta=xi=zeros(1,n-1);
while(c1==0)
{
 Sj1=Sj2=1;
 sigma=sigma|zeros(1,n-1);
 tau=tau|zeros(1,n-1);
 zeta=zeta|zeros(1,n-1);
 xi=xi|zeros(1,n-1);
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u2[i-1])*exp(-2*K*j*(K*j-u2[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u2[i-1])*exp(-2*K*j*(K*j+u2[i-1])/s[i-1]);
  }
  if(u2[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u1[i-1])*exp(-2*K*j*(K*j-u1[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u1[i-1])*exp(-2*K*j*(K*j+u1[i-1])/s[i-1]);
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   sigma[j][i-1]=(sigmabar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+sigmabar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
   tau[j][i-1]=(taubar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+taubar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
  }		  // println(xi);
 }
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   Sj2*=1-(1/u2[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);		//	println(u2[i-1]~K);
   Sj1*=1-(1/u2[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u2[i-1]==0)
  {
   Sj2*=1-(1/u1[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);			//			println(u1[i-1]~K);
   Sj1*=1-(1/u1[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   Sj2*=(1-sumc(sigma[0:(j-1)][i-1]-tau[0:(j-1)][i-1])-sigma[j][i-1])/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));
   Sj1*=(1-sumc(sigma[0:j][i-1]-tau[0:j][i-1]))/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));
  }				 
 }		//		 println(u~Sj1~Sj2);
 
 if(u<Sj2){DI=1;c1=1;}
 else if(u>Sj1){DI=0;c1=1;}
 else{j+=1;}
}

}


else
{


sigma=tau=zeta=xi=zeros(1,n-1);

for(j=1;j<=(Kh-1);++j)
{
 Sj1=Sj2=1;
 sigma=sigma|zeros(1,n-1);
 tau=tau|zeros(1,n-1);
 zeta=zeta|zeros(1,n-1);
 xi=xi|zeros(1,n-1);
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u2[i-1])*exp(-2*K*j*(K*j-u2[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u2[i-1])*exp(-2*K*j*(K*j+u2[i-1])/s[i-1]);
  }
  if(u2[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u1[i-1])*exp(-2*K*j*(K*j-u1[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u1[i-1])*exp(-2*K*j*(K*j+u1[i-1])/s[i-1]);
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   sigma[j][i-1]=(sigmabar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+sigmabar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
   tau[j][i-1]=(taubar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+taubar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
  }
 }
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   Sj2*=1-(1/u2[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);		//	println(u2[i-1]~K);
   Sj1*=1-(1/u2[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u2[i-1]==0)
  {
   Sj2*=1-(1/u1[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);			//	println(u1[i-1]~K);
   Sj1*=1-(1/u1[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   Sj2*=(1-sumc(sigma[0:(j-1)][i-1]-tau[0:(j-1)][i-1])-sigma[j][i-1])/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));			
   Sj1*=(1-sumc(sigma[0:j][i-1]-tau[0:j][i-1]))/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));
  }
 }		//							println(u~Sj1~Sj2);

} //for j

j=Kh;
while(c1==0)
{
 Sj1=Sj2=1;
 sigma=sigma|zeros(1,n-1);
 tau=tau|zeros(1,n-1);
 zeta=zeta|zeros(1,n-1);
 xi=xi|zeros(1,n-1);
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u2[i-1])*exp(-2*K*j*(K*j-u2[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u2[i-1])*exp(-2*K*j*(K*j+u2[i-1])/s[i-1]);
  }
  if(u2[i-1]==0)
  {
   zeta[j][i-1]=(2*K*j-u1[i-1])*exp(-2*K*j*(K*j-u1[i-1])/s[i-1]);
   xi[j][i-1]=(2*K*j+u1[i-1])*exp(-2*K*j*(K*j+u1[i-1])/s[i-1]);
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   sigma[j][i-1]=(sigmabar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+sigmabar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
   tau[j][i-1]=(taubar(sg[i-1],u1g[i-1],u2g[i-1],Kg,j)+taubar(sg[i-1],-u1g[i-1],-u2g[i-1],Kg,j));
  }
 }
 for(i=1;i<=(n-1);++i)
 {
  if(u1[i-1]==0)
  {
   Sj2*=1-(1/u2[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);		//	println(u2[i-1]~K);
   Sj1*=1-(1/u2[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u2[i-1]==0)
  {
   Sj2*=1-(1/u1[i-1])*(sumc(zeta[0:(j-1)][i-1]-xi[0:(j-1)][i-1])+zeta[j][i-1]);			//	println(u1[i-1]~K);
   Sj1*=1-(1/u1[i-1])*(sumc(zeta[0:(j)][i-1]-xi[0:(j)][i-1]));
  }
  if(u1[i-1]>0 && u2[i-1]>0)
  {
   Sj2*=(1-sumc(sigma[0:(j-1)][i-1]-tau[0:(j-1)][i-1])-sigma[j][i-1])/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));			
   Sj1*=(1-sumc(sigma[0:j][i-1]-tau[0:j][i-1]))/(1-exp(-2*u1[i-1]*u2[i-1]/s[i-1]));
  }
 }		//							println(u~Sj1~Sj2);
 
 if(u<Sj2){DI=1;c1=1;}
 else if(u>Sj1){DI=0;c1=1;}
 else{j+=1;}
}

} //else


return(DI);
  
}	   
  