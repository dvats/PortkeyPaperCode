#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>

layerass(t0,x0,T,XT,a1,b1,dl,MINI,MAXI)

{

decl a,b,x,y,c1,c2,c3,s,u1,u2,K,i,sigma,tau,u,I,j,Sj1,Sj2,rtn,ind1,ind2;
decl xbar,ybar,bb1,nb,l,S3,S2,sigma2,tau2,Sj12,Sj22,sigma3,tau3,Sj13,Sj23;
decl aa1,na;

x=x0;
y=XT;

s=T-t0;
u1=(x-y)/2;
u2=(y-x)/2;
K=fabs(y-x)/2;

xbar=min(x,y);
ybar=max(x,y);

u=ranu(1,1);

ind1=ind2=0;

if((xbar-a1)<MINI){ind1=1;}
if((ybar+b1)>MAXI){ind2=1;}

/////////////////////////////////////////

if(ind1==1 && ind2==1)
{
rtn=((MINI)~(xbar))|((ybar)~(MAXI));
}

//////////////////////////////////////////

if(ind1==1 && ind2==0)
{

a=0|(xbar-MINI);
b=0|b1;
bb1=MAXI-ybar;

while(ybar+(maxc(b)+dl)<MAXI){b=b|(maxc(b)+dl);}b=b|bb1;
nb=rows(b);

c2=0;
i=1;

l=0; //counter for D's j
S3=0;

while(c2==0)
{
 j=1;
 c3=0;
 S2=0;
 while(c3==0)
 {
  sigma2=(sigmabar(s,u1+((a[1]-b[i])/2),u2+((a[1]-b[i])/2),K+((a[1]+b[i])/2),j)+sigmabar(s,-u1-((a[1]-b[i])/2),-u2-((a[1]-b[i])/2),K+((a[1]+b[i])/2),j));
  tau2=(taubar(s,u1+((a[1]-b[i])/2),u2+((a[1]-b[i])/2),K+(a[1]+b[i])/2,j)+taubar(s,-u1-((a[1]-b[i])/2),-u2-((a[1]-b[i])/2),K+((a[1]+b[i])/2),j));
  Sj22=1-S2-sigma2;
  S2+=sigma2-tau2;
  Sj12=Sj22+tau2;
				  
  if(l<j)
  {
   l+=1;
   sigma3=(sigmabar(s,u1+((a[1]-bb1)/2),u2+((a[1]-bb1)/2),K+((a[1]+bb1)/2),j)+sigmabar(s,-u1-((a[1]-bb1)/2),-u2-((a[1]-bb1)/2),K+((a[1]+bb1)/2),j));
   tau3=(taubar(s,u1+((a[1]-bb1)/2),u2+((a[1]-bb1)/2),K+(a[1]+bb1)/2,j)+taubar(s,-u1-((a[1]-bb1)/2),-u2-((a[1]-bb1)/2),K+((a[1]+bb1)/2),j));
   Sj23=1-S3-sigma3;
   S3+=sigma3-tau3;
   Sj13=Sj23+tau3;
  }

  if(u<(Sj22/Sj13)){c3=1;c2=1;rtn=((MINI)~(xbar))|((ybar+b[i-1])~(ybar+b[i]));}
  else if(u>(Sj12/Sj23)){c3=1; if(i<(nb-2)){i+=1;} if(i==(nb-2)){c2=1;rtn=((MINI)~(xbar))|((ybar+b[i])~(MAXI));}}
  else{j+=1;}

 }
}

}


//////////////////////////////////////////

if(ind1==0 && ind2==1)
{

a=<0>|a1;
b=<0>|(MAXI-ybar);
aa1=(xbar-MINI);

while(xbar-(maxc(a)+dl)>MINI){a=a|(maxc(a)+dl);}a=a|(xbar-MINI);
na=rows(a);

c2=0;
i=1;

l=0; //counter for D's j
S3=0;

while(c2==0)
{
 j=1;
 c3=0;
 S2=0;
 while(c3==0)
 {
  sigma2=(sigmabar(s,u1+((a[i]-b[1])/2),u2+((a[i]-b[1])/2),K+((a[i]+b[1])/2),j)+sigmabar(s,-u1-((a[i]-b[1])/2),-u2-((a[i]-b[1])/2),K+((a[i]+b[1])/2),j));
  tau2=(taubar(s,u1+((a[i]-b[1])/2),u2+((a[i]-b[1])/2),K+(a[i]+b[1])/2,j)+taubar(s,-u1-((a[i]-b[1])/2),-u2-((a[i]-b[1])/2),K+((a[i]+b[1])/2),j));
  Sj22=1-S2-sigma2;
  S2+=sigma2-tau2;
  Sj12=Sj22+tau2;
				  
  if(l<j)
  {
   l+=1;
   sigma3=(sigmabar(s,u1+((aa1-b[1])/2),u2+((aa1-b[1])/2),K+((aa1+b[1])/2),j)+sigmabar(s,-u1-((aa1-b[1])/2),-u2-((aa1-b[1])/2),K+((aa1+b[1])/2),j));
   tau3=(taubar(s,u1+((aa1-b[1])/2),u2+((aa1-b[1])/2),K+(aa1+b[1])/2,j)+taubar(s,-u1-((aa1-b[1])/2),-u2-((aa1-b[1])/2),K+((aa1+b[1])/2),j));
   Sj23=1-S3-sigma3;
   S3+=sigma3-tau3;
   Sj13=Sj23+tau3;
  }

  if(u<(Sj22/Sj13)){c3=1;c2=1;rtn=((xbar-a[i])~(xbar-a[i-1]))|((ybar)~(MAXI));}
  else if(u>(Sj12/Sj23)){c3=1; if(i<(na-2)){i+=1;} if(i==(na-2)){c2=1;rtn=((MINI)~(xbar-a[i]))|((ybar)~(MAXI));}}
  else{j+=1;}

 }
}

}


//////////////////////////////////////////

if(ind1==0 && ind2==0)
{

a=0|a1;
b=0|b1;
aa1=(xbar-MINI);
bb1=MAXI-ybar;

while((maxc(b)+dl)<MAXI){b=b|(maxc(b)+dl);}b=b|bb1;
nb=rows(b);

while(xbar-(maxc(a)+dl)>MINI){a=a|(maxc(a)+dl);}a=a|aa1;
na=rows(a);

c2=0;
I=0;
i=1;

l=0; //counter for D's j
S3=0;

while(c2==0)
{
 j=1;
 c3=0;
 S2=0;
 while(c3==0)
 {						 
  sigma2=(sigmabar(s,u1+((a[i]-b[i])/2),u2+((a[i]-b[i])/2),K+((a[i]+b[i])/2),j)+sigmabar(s,-u1-((a[i]-b[i])/2),-u2-((a[i]-b[i])/2),K+((a[i]+b[i])/2),j));
  tau2=(taubar(s,u1+((a[i]-b[i])/2),u2+((a[i]-b[i])/2),K+(a[i]+b[i])/2,j)+taubar(s,-u1-((a[i]-b[i])/2),-u2-((a[i]-b[i])/2),K+((a[i]+b[i])/2),j));
  Sj22=1-S2-sigma2;
  S2+=sigma2-tau2;
  Sj12=Sj22+tau2;
				  
  if(l<j)
  {
   l+=1;
   sigma3=(sigmabar(s,u1+((aa1-bb1)/2),u2+((aa1-bb1)/2),K+((aa1+bb1)/2),j)+sigmabar(s,-u1-((aa1-bb1)/2),-u2-((aa1-bb1)/2),K+((aa1+bb1)/2),j));
   tau3=(taubar(s,u1+((aa1-bb1)/2),u2+((aa1-bb1)/2),K+(aa1+bb1)/2,j)+taubar(s,-u1-((aa1-bb1)/2),-u2-((aa1-bb1)/2),K+((aa1+bb1)/2),j));
   Sj23=1-S3-sigma3;
   S3+=sigma3-tau3;
   Sj13=Sj23+tau3;
  }

  if(u<(Sj22/Sj13))
  {
   c3=1;c2=1;

   if(i<=(nb-2) && i<=(na-2)){rtn=((xbar-a[i])~(xbar-a[i-1]))|((ybar+b[i-1])~(ybar+b[i]));}
   if(i>(nb-2) && i<=(na-2)){rtn=((xbar-a[i])~(xbar-a[i-1]))|((ybar+b[nb-2])~(MAXI));}
   if(i>(na-2) && i<=(nb-2)){rtn=((MINI)~(xbar-a[na-2]))|((ybar+b[i-1])~(ybar+b[i]));}
  }
   
  else if(u>(Sj12/Sj23)){c3=1; if(i<(nb-2) && i<(na-2)){i+=1;} if(i>=(nb-2) && i<(na-2)){i+=1;b=b|maxc(b);} if(i>=(na-2) && i<(nb-2)){i+=1;a=a|maxc(a);}
		  if(i>=(na-2) && i>=(nb-2)){c2=1;rtn=((MINI)~(xbar-a[na-2]))|((ybar+b[nb-2])~(MAXI));}}
  else{j+=1;}

 }
}

}



return (rtn);
  
}	   
  