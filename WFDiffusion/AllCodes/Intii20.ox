

Intii20(const X1, const dtt, const g1n, const g2n, const g1o, const g2o, const cp, const cphi)

{
  decl n,t,a1,dl,X11,x,y,xt,xtt,l,i,j,ind0,ind1,I,m,M,l1,l21,l22,mm,MM,MMM,indm,sindm,I1,I2,an,bn,ao,bo,xm,xM,rtn;

//t=range(0,dt,dtt)';
n=rows(X1);

//a1=0.25*sqrt(dtt*(1-dtt));
//dl=0.25*sqrt(dtt*(1-dtt));

an=-0.5*(g1n-1);
bn=g1n*g2n/2-1/4;

ao=-0.5*(g1o-1);
bo=g1o*g2o/2-1/4;

ind0=0;
if( ( (an+bn)/2+(an+bn)^2 - ((ao+bo)/2+(ao+bo)^2) )>0 ){ind0=1;}
ind1=0;
if( ( -bn/2+bn^2 - (-bo/2+bo^2) )>0 ){ind1=1;}


m=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
M=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];


if(ind0==1 && ind1==0)
{

xm=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
m=phi2(g1n,g2n,g1o,g2o,xm);

xM=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];
M=phi2(g1n,g2n,g1o,g2o,xM);

I=(-dtt*m);
MM=M-m;

}

else if(ind0==0 && ind1==1)
{

xm=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
xM=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];

m=phi2(g1n,g2n,g1o,g2o,xM);
M=phi2(g1n,g2n,g1o,g2o,xm);

I=(-dtt*m);
MM=M-m;

}

else if(ind0==1 && ind1==1)
{

xm=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
xM=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];

m=M=zeros(n,1);
for(i=1;i<=n;++i)
{

if(xm[i-1]>cp)
{
m[i-1]=phi2(g1n,g2n,g1o,g2o,xm[i-1]);
M[i-1]=phi2(g1n,g2n,g1o,g2o,xM[i-1]);
}

else if(xM[i-1]<cp)
{
m[i-1]=phi2(g1n,g2n,g1o,g2o,xM[i-1]);
M[i-1]=phi2(g1n,g2n,g1o,g2o,xm[i-1]);
}

else
{
m[i-1]=cphi;
M[i-1]=max(phi2(g1n,g2n,g1o,g2o,xM[i-1]),phi2(g1n,g2n,g1o,g2o,xm[i-1]));
}

}

I=(-dtt*m);
MM=M-m;

}

else
{

xm=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
xM=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];

m=M=zeros(n,1);
for(i=1;i<=n;++i)
{

if(xm[i-1]>cp)
{
m[i-1]=phi2(g1n,g2n,g1o,g2o,xM[i-1]);
M[i-1]=phi2(g1n,g2n,g1o,g2o,xm[i-1]);
}

else if(xM[i-1]<cp)
{
m[i-1]=phi2(g1n,g2n,g1o,g2o,xm[i-1]);
M[i-1]=phi2(g1n,g2n,g1o,g2o,xM[i-1]);
}

else
{
m[i-1]=min(phi2(g1n,g2n,g1o,g2o,xM[i-1]),phi2(g1n,g2n,g1o,g2o,xm[i-1]));
M[i-1]=cphi;
}

}

I=(-dtt*m);
MM=M-m;

}


rtn=m~M~MM~I;

return rtn;


}