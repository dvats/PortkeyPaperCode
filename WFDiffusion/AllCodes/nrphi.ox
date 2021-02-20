
nrphi(const g1, const g2, const x0, const eps)

{
  decl dif,f,df,x;


dif=1;
x=x0;

while(dif>eps)
{

f=dphi(g1,g2,x);
df=ddphi(g1,g2,x);

dif=f/df;

x=x-dif;

dif=fabs(dif);

//println(x);

}

  
return x;


}