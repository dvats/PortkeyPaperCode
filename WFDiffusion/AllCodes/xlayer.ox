#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>

xlayer(t0,x,T,y,a1,b1,tx,MINI,MAXI)

{

decl t,xbar,ybar,c1,ind1,Li,Ls,m,ai,aim,bi,bim,Fm,FM,FF,Xt,DI,LI,u,s0,XX,c2;

t=T-t0;
xbar=min(x,y);
ybar=max(x,y);


 aim=0;
 ai=a1;
 bim=0;
 bi=b1;
 Li=(xbar-a1)|(xbar);
 Ls=(ybar)|(ybar+b1);

Fm = 1 - exp( ( (y-x)^2 - (y+x-2*Li[0])^2 ) / ( 2*(T-t0) ) );
FM = 1 - exp( ( (y-x)^2 - (y+x-2*Ls[1])^2 ) / ( 2*(T-t0) ) );
FF=min(Fm,FM);

c1=0;
while(c1==0)
{
 ind1=ranbinomial(1,1,1,0.5);
// println(ind1);
 if(ind1==1)
 {
  m=condminim(t0,x,T,y,Li[0],Li[1]);
  c2=0;
  while(c2==0)
  {
   Xt=xminim(t0,x,T,y,m[1],m[0],tx);
   if(maxc(Xt)<MAXI){c2=1;}
  }
  s0=(t0|tx|m[0]|T);
  XX=(s0)~(x|Xt|m[1]|y);
  XX=sortbyc(XX,0);
  DI=deltaKmin(t0,x,T,y,XX[][0],XX[][1],m[1],bi);	  //println(DI);
  if(DI==1)
  {
   c1=ranbinomial(1,1,1,FF/(Fm+FM));
  }
 }

 if(ind1==0)
 {
  m=condmaxim(t0,x,T,y,Ls[0],Ls[1]);
  c2=0;
  while(c2==0)
  {
   Xt=xmaxim(t0,x,T,y,m[1],m[0],tx);
   if(minc(Xt)>MINI){c2=1;}
  }
  s0=(t0|tx|m[0]|T);
  XX=(s0)~(x|Xt|m[1]|y);
  XX=sortbyc(XX,0);
  DI=deltaKmax(t0,x,T,y,XX[][0],XX[][1],m[1],ai);				   //println(DI);
  if(DI==1)
  {
   c1=ranbinomial(1,1,1,FF/(Fm+FM));
  }
 }
}



return Xt;
  
}	   
  