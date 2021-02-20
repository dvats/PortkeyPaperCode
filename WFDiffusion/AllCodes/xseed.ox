
xseed(const X1, const dt, const dtt, const MINI, const MAXI, const cp, const g1, const g2)

{
  decl n,nn,t,dl,X11,X2,x,y,xt,xtt,l,i,j,a11,b11,c,xx,mm,tt,ind,mm1,mm2,X2na,X2nb;

n=rows(X1);
nn=dt/dtt;
t=range(dtt,dt-dtt,dtt)';


a11=b11=sqrt(dtt)/4;
dl=sqrt(dtt)/8;


X11=zeros(n*nn,12);
X2=zeros(2*n*nn,14);

for (i=1;i<=n;++i)
{
//			  println(i);
x=X1[i-1][0];
y=X1[i-1][1];

xt=zeros(nn,1);xt[0]=x;
tt=0|t;
for(j=1;j<=(nn-1);++j)
{
c=0;
while(c==0)
{
xt[j]=BB2(xt[j-1],tt[j-1],y,dt,tt[j]);
if(xt[j]>MINI && xt[j]<MAXI){ind=DI(tt[j]-tt[j-1],xt[j-1],xt[j],MINI,MAXI)*DI(dt-tt[j],xt[j],y,MINI,MAXI);if(ind==1){c=1;}}
}									   //println(xt);
}
xt=xt[1:];
xx=(x|xt)~(xt|y);
X11[((i-1)*nn):(i*nn-1)][:3]=constant(i-1,nn,1)~range(0,nn-1)'~xx;


xtt=zeros(nn,5);

for(j=1;j<=nn;++j)
{
 l=layerass(0,xx[j-1][0],dtt,xx[j-1][1],a11,b11,dl,MINI,MAXI);
 xtt[j-1][]=xlayerass0(0,xx[j-1][0],dtt,xx[j-1][1],l[0][0],l[0][1],l[1][0],l[1][1])~1;
}

X2na=constant(i-1,nn,1)~range(0,nn-1)'~zeros(nn,1)~xtt[][1]~xx[][0]~xtt[][2]~ones(nn,1)~zeros(nn,1)~zeros(nn,6);
X2nb=constant(i-1,nn,1)~range(0,nn-1)'~xtt[][1]~constant(dtt,nn,1)~xtt[][2]~xx[][1]~ones(nn,2)~zeros(nn,6);
X2[(2*(i-1)*nn):(2*i*nn-1)][]=sortbyc(X2na|X2nb,<1,2>);

X11[((i-1)*nn):(i*nn-1)][4:7]=xtt[][:3];
X11[((i-1)*nn):(i*nn-1)][8:]=Inti0(xtt,dtt,g1,g2,cp);

}

return {X11,X2};


}