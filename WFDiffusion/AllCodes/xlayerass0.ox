#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>

xlayerass0(t0,x,T,y,m1,m2,M1,M2)

{

decl xbar,ybar,c1,ind1,Li,Ls,m,ai,aim,Xt,DI,LI,u,s0,XX,bim,bi,BM,Bm,mB,c2;

xbar=min(x,y);
ybar=max(x,y);


Li=(m1)|(m2);
Ls=(M1)|(M2);

aim=xbar-m2;
ai=xbar-m1;
bim=M1-ybar;
bi=M2-ybar;


BM=exp(((y-x)^2-(y+x-2*(xbar-bim))^2)/(2*(T-t0)))-exp(((y-x)^2-(y+x-2*(xbar-bi))^2)/(2*(T-t0)));
Bm=exp(((y-x)^2-(y+x-2*(xbar-aim))^2)/(2*(T-t0)))-exp(((y-x)^2-(y+x-2*(xbar-ai))^2)/(2*(T-t0)));  // println(a1);	//println(x~y~xbar~ybar~aim~ai~bim~bi~T~t0); println(BM~Bm);
mB=min(BM,Bm);
																											 //  println(BM~Bm);
//println(b1);

c1=0;
while(c1==0)
{					 // println("try1");
 ind1=ranbinomial(1,1,1,0.5);
// println(ind1);
 if(ind1==1)
 {
  m=condminim(t0,x,T,y,Li[0],Li[1]);	  // println(m);

  DI=deltaK(m[1],m[0]-t0,x,Ls[1])*deltaK(m[1],T-m[0],y,Ls[1]);

  if(DI==1)
  {
   if(fabs(M1-ybar)<10^(-12)){LI=1;} else{LI=deltaL(m[1],m[0]-t0,x,Ls[0],Ls[1])*deltaL(m[1],T-m[0],y,Ls[0],Ls[1]);}
   if(LI==0){c1=ranbinomial(1,1,1,mB/BM);}
   if(LI==1){c1=ranbinomial(1,1,1,mB/(Bm+BM));}
  }
 }

 if(ind1==0)
 {				 // println("try2");
  m=condmaxim(t0,x,T,y,Ls[0],Ls[1]);		 //  println(m);
  
  DI=deltaK(m[1],m[0]-t0,x,Li[0])*deltaK(m[1],T-m[0],y,Li[0]);
  
  if(DI==1)
  {
   if(fabs(m2-xbar)<10^(-12)){LI=1;} else{LI=deltaL(m[1],m[0]-t0,x,Li[1],Li[0])*deltaL(m[1],T-m[0],y,Li[1],Li[0]);}
   if(LI==0){c1=ranbinomial(1,1,1,mB/Bm);}
   if(LI==1){c1=ranbinomial(1,1,1,mB/(Bm+BM));}
  }
 }
}



return ind1~m~(ind1*Ls[1]+(1-ind1)*Li[0]);
  
}	   
  