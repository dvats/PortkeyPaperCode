
xmaxim(t0,x0,T,XT,m,tm,t)

{

decl n,u33,Xt,i1,i2,w1,w2,w3,i,ts,delta,r;

n=rows(t);
if(n==0){n=1;}

Xt=x0|zeros(n,1)|XT;
i1=i2=0;
u33=t0|t|T;

for(i=1;i<=n;++i)
{
if(u33[i]<tm){i1+=1;}
if(u33[i]>tm){i2+=1;}
}
w1=w2=w3=0;
ts=0;
for(i=1;i<=i1;++i)
{
ts=ts|((tm-u33[i1-i+1])/(tm-t0));
w1=BB2(w1,ts[i-1],0,1,ts[i]);
w2=BB2(w2,ts[i-1],0,1,ts[i]);
w3=BB2(w3,ts[i-1],0,1,ts[i]);
delta=(m-x0)/(sqrt(tm-t0));
r=sqrt((delta*ts[i]+w1)^2+w2^2+w3^2);
Xt[i1-i+1]=-sqrt(tm-t0)*r+m;
}

w1=w2=w3=0;
ts=1;
for(i=1;i<=i2;++i)
{
ts=ts|((u33[n+1-i]-tm)/(T-tm));
w1=BB2(0,0,w1,ts[i-1],ts[i]);
w2=BB2(0,0,w2,ts[i-1],ts[i]);
w3=BB2(0,0,w3,ts[i-1],ts[i]);
delta=(m-XT)/(sqrt(T-tm));
r=sqrt((delta*ts[i]+w1)^2+w2^2+w3^2);
Xt[n+1-i]=-sqrt(T-tm)*r+m;
}

Xt=Xt[1:n];

return Xt;
  
}	   
  