

XBarkerPAR(const X0, X1, const X2, const T, const dt, const dtt, const g1, const g2, const cp, const MINI, const MAXI, const beta)

{
  decl nn,t,ind0,ind1,a11,b11,dl,cphi,c1,c2,c3,i,n,x,y,indX1,indX2,X1o,X2o,c,xt,xx,X1n,X2n,xtt,j,l,s0,s1,p,MMMo,MMMn,j1,j2,C1,np,xn,yn,mn,Mn,l1n,l21n,l22n,mmn,MMn,ln,un,indm,sindm,inde,nnp,cc,up,indd,phil,phiu,ind2,ind3,upp,inddd,x2,phix,k1,k2,xo,yo,mo,Mo,l1o,l21o,l22o,mmo,MMo,lo,uo,mm1,mm2,tt,ind,a,b,cppo,cppn,temp,S;

decl cont=0;	decl cont1=0;   decl ind4,X2na,X2nb;
  
nn=dt/dtt;
t=range(dtt,dt-dtt,dtt)';

a=-0.5*(g1-1);
b=g1*g2/2-1/4;

ind0=0;
if( ((a+b)/2+(a+b)^2)>0 ){ind0=1;}
ind1=0;
if( (-b/2+b^2)>0 ){ind1=1;}

a11=b11=sqrt(dtt)/4;
dl=sqrt(dtt)/8;


cphi=phi(g1,g2,cp);


n=rows(X0);

serial decl X22=zeros(1,14);

parallel for (i=1;i<=n;++i)
{
// println(i);

/////////////////Proposal////////////////

x=X0[i-1][0];
y=X0[i-1][1];

indX1=(((i-1)*nn)|(i*nn-1));
X1o=X1[indX1[0]:indX1[1]][];
indX2=vecindex(X2[][0],i-1);
indX2=(indX2[0]|indX2[rows(indX2)-1]);
X2o=X2[indX2[0]:indX2[1]][];

xt=zeros(nn,1);xt[0]=x;
tt=0|t;

																  //if(i>=0){println("a1");}
xt[1:(nn-1)]=xlayer(0,x,dt,y,min(x,y)-MINI,MAXI-max(x,y),t,MINI,MAXI);  //simulate xt conditional on asymetric layers
																  //if(i>=0){println("b1");}
xt=xt[1:];
xx=(x|xt)~(xt|y);
X1n=constant(i-1,nn,1)~range(0,nn-1)'~xx;

xtt=zeros(nn,5);
for(j=1;j<=nn;++j)
{													 //if(i>=0){println("a2");}
 l=layerass(0,xx[j-1][0],dtt,xx[j-1][1],a11,b11,dl,MINI,MAXI);		   //if(i>=0){println("b2");}
 xtt[j-1][]=xlayerass0(0,xx[j-1][0],dtt,xx[j-1][1],l[0][0],l[0][1],l[1][0],l[1][1]);
}

X2na=constant(i-1,nn,1)~range(0,nn-1)'~zeros(nn,1)~xtt[][1]~xx[][0]~xtt[][2]~ones(nn,1)~zeros(nn,1)~zeros(nn,6);
X2nb=constant(i-1,nn,1)~range(0,nn-1)'~xtt[][1]~constant(dtt,nn,1)~xtt[][2]~xx[][1]~ones(nn,2)~zeros(nn,6);
X2n=sortbyc(X2na|X2nb,<1,2>);

X1n=X1n~xtt[][:3];

X1n=X1n~Inti0(xtt,dtt,g1,g2,cp);

//if(i==35){println(X1n);println(X2n);}

////////////////////////////////////////////

s0=sumc(X1o[][11]);
s1=sumc(X1n[][11]);
p=1/(1+exp(s0-s1));
MMMo=X1o[][10];
MMMn=X1n[][10];

c=0;
j1=j2=0;				 //if(i==20){println(s0~s1~p);println(X1o);println(X2o);}
while(c==0)
{cont+=1;

S = ranbinomial(1,1,1,beta);
if(S == 0)
{
c = 1;
X22|=X2o;	cont1+=1; //println("coin 00"); println("X",i);
//code to reject proposed value. That is, new bridge for each i is the previous bridge. (*****)
}

else
{

 C1=ranbinomial(1,1,1,p);			//if(i==20){println(C1);}

 if(C1==1)
 {
  np=ranpoisson(nn,1,dtt*MMMn);	  //if(minc(dtt*MMMn)<0 || maxc(dtt*MMMn)>100){println((dtt*MMMn)~np);}
  if(sumc(np)==0){c=1;}
  else
  {
   c=1;

   inde=round(vecindex(np));
   np=np[inde];
   nnp=rows(np);

   j=1;
   cc=0;
  while(j<=nnp && cc==0)
  {
   indd=vecindex(X1n[][1],inde[j-1]);
   up=runif(np[j-1],2,zeros(np[j-1],2),(zeros(np[j-1],1)+dtt)~(zeros(np[j-1],1)+MMMn[indd]));
   up=sortbyc(up,0);

    inddd=vecindex(X2n[][1],inde[j-1]);		 //if(i>=0){println("a3");}//println(X1n[indd][]);println(X2n[inddd][]);println(up[][0]);println(dtt);}
	[x2,xt]=BBx(X1n[indd][],X2n[inddd][],up[][0],dtt);		  //if(i>=0){println("b3");}
//	xt= xt + ((dtt-up[][0])/dtt)*X1n[indd][2] + (up[][0]/dtt)*X1n[indd][3];

    X2n=sortbyc(dropr(X2n,inddd)|x2,<1,2>);		//	renew X2
	phix=phi(g1,g2,xt)-X1n[indd][8];
	if(up[][1]>phix){j+=1;}
	else{cc=1;c=0;}
   
  }

  } //else np=0

  if(c==1)
  {

   X22|=X2n;
   
   X1[indX1[0]:indX1[1]][]=X1n;
   
//   X1=insertr(X1,indX1[0],k2);
//   X1=dropr(X1,range((indX1[0]+k2),(indX1[0]+k2+k1-1),1));
//   X1[indX1[0]:(indX1[0]+k2-1)][]=X1n;
  }


 } //if C1=1




 else //(C1=0)
 {

  np=ranpoisson(nn,1,dtt*MMMo);
  if(sumc(round(np))==0){c=1;}
  else
  {
   c=1;

   inde=vecindex(np);
   np=np[inde];
   nnp=rows(np);

   j=1;
   cc=0;
  while(j<=nnp && cc==0)
  {
   indd=vecindex(X1o[][1],inde[j-1]);
   up=runif(np[j-1],2,zeros(np[j-1],2),(zeros(np[j-1],1)+dtt)~(zeros(np[j-1],1)+MMMo[indd]));
   up=sortbyc(up,0);

    inddd=vecindex(X2o[][1],inde[j-1]);		//if(i>=0){println("a4");println(inde[j-1],inddd);}
	[x2,xt]=BBx(X1o[indd][],X2o[inddd][],up[][0],dtt);		   //if(i>=0){println("b4");}
		
    X2o=sortbyc(dropr(X2o,inddd)|x2,<1,2>);
	phix=phi(g1,g2,xt)-X1o[indd][8];
	if(up[][1]>phix){j+=1;}
	else{cc=1;c=0;}

  }

  } //else np=0

    if(c==1)
  {

   X22|=X2o;
  }
 
 } //else (C1=0)

} //else (S=0)
 
} //while c=0


}//for i

X22=X22[1:][];

X22=sortbyc(X22,<0,1,2>);

return {X1,X22,cont/n,cont1};


}



