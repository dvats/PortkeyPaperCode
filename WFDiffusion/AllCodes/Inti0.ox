
Inti0(const xtt, const dtt, const g1, const g2, const cp)

{
  decl ind0,ind1,cphi,m,M,xm,xM,MM,I,a,b,rtn,n,i;


a=-0.5*(g1-1);
b=g1*g2/2-1/4;

ind0=0;
if( ((a+b)/2+(a+b)^2)>0 ){ind0=1;}
ind1=0;
if( (-b/2+b^2)>0 ){ind1=1;}

m=M=MM=I=zeros(rows(xtt),1);

cphi=phi(g1,g2,cp);

n=rows(xtt);

if(ind0==1 && ind1==0)
{

xm=xtt[][0].*xtt[][2] + (1-xtt[][0]).*xtt[][3];
m=phi(g1,g2,xm);

xM=xtt[][0].*xtt[][3] + (1-xtt[][0]).*xtt[][2];
M=phi(g1,g2,xM);

I=(-dtt*m);
MM=M-m;

}

else if(ind0==0 && ind1==1)
{

xm=xtt[][0].*xtt[][2] + (1-xtt[][0]).*xtt[][3];
xM=xtt[][0].*xtt[][3] + (1-xtt[][0]).*xtt[][2];

m=phi(g1,g2,xM);
M=phi(g1,g2,xm);

I=(-dtt*m);
MM=M-m;

}

else if(ind0==1 && ind1==1)
{

xm=xtt[][0].*xtt[][2] + (1-xtt[][0]).*xtt[][3];
xM=xtt[][0].*xtt[][3] + (1-xtt[][0]).*xtt[][2];

m=M=zeros(n,1);
for(i=1;i<=n;++i)
{

if(xm[i-1]>cp)
{
m[i-1]=phi(g1,g2,xm[i-1]);
M[i-1]=phi(g1,g2,xM[i-1]);
}

else if(xM[i-1]<cp)
{
m[i-1]=phi(g1,g2,xM[i-1]);
M[i-1]=phi(g1,g2,xm[i-1]);
}

else
{
m[i-1]=cphi;
M[i-1]=max(phi(g1,g2,xM[i-1]),phi(g1,g2,xm[i-1]));
}

}

I=(-dtt*m);
MM=M-m;

}

else
{

xm=xtt[][0].*xtt[][2] + (1-xtt[][0]).*xtt[][3];
xM=xtt[][0].*xtt[][3] + (1-xtt[][0]).*xtt[][2];

m=M=zeros(n,1);
for(i=1;i<=n;++i)
{

if(xm[i-1]>cp)
{
m[i-1]=phi(g1,g2,xM[i-1]);
M[i-1]=phi(g1,g2,xm[i-1]);
}

else if(xM[i-1]<cp)
{
m[i-1]=phi(g1,g2,xm[i-1]);
M[i-1]=phi(g1,g2,xM[i-1]);
}

else
{
m[i-1]=min(phi(g1,g2,xM[i-1]),phi(g1,g2,xm[i-1]));
M[i-1]=cphi;
}

}

I=(-dtt*m);
MM=M-m;

}


rtn=m~M~MM~I;

return rtn;


}