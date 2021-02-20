
BB2(const x0,const t0,const XT,const T,const t)

{
  decl Xt,i,n,tt;

n=rows(t);
if(n==0){n=1;}

Xt=zeros(n+1,1);
Xt[0]=x0;
tt=t0|t;

for(i=1;i<=n;++i)
{
 Xt[i]=rnorm(1,1,Xt[i-1]+(XT-Xt[i-1])*((tt[i]-tt[i-1])/(T-tt[i-1])),((tt[i]-tt[i-1])*(T-tt[i]))/(T-tt[i-1]));
}

Xt=Xt[1:n];

return Xt;

		  
}


  
