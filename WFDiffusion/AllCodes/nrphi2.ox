
nrphi2(const g1n, const g2n, const g1o, const g2o, const x0, const eps)

{
  decl dif,f,df,x;


dif=1;
x=x0;

while(dif>eps)
{

f=dphi2(g1n,g2n,g1o,g2o,x);
df=ddphi2(g1n,g2n,g1o,g2o,x);

dif=f/df;

x=x-dif;

dif=fabs(dif);

//println(x);

}

  
return x;


}