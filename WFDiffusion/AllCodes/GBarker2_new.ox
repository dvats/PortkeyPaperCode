
GBarker2(const X1, X2, const dtt, const g1o, const g2o, const cpo, const s1, const s2, const MINI, const MAXI, const X0, const XT, const beta)

{
  decl S,c,gamao,gaman,g1n,g2n,n,ind0,ind1,cp,cphi,c1n,c2n,c3n,c4n,c1o,c2o,c3o,c4o,c1,c2,c3,c4,XX1,s0,p,MMM,j,x,y,m,M,l1,l21,l22,mm,MM,l,u,C1,i,cB,np,X11,indm,sindm,up,phil,phiu,ind2,indX2,X2i,ind3,upp,xt,phix,k1,k2,gama,mm1,mm2,an,bn,ao,bo,cpp,a,b,dAn,dAo,cppp,temp;

decl cont=0;
decl cont1=0;
  
gamao=g1o|g2o;

//S=(s1~(r*sqrt(s1*s2)))|((r*sqrt(s1*s2))~s2);

//c=0;
//while(c==0)
//{
//gaman=rmvnorm(2,gamao,S);
gaman=((g1o|g2o)-(s1|s2)+2*(s1|s2).*ranu(2,1));			//runif(2,1,gamao-(s1|s2),gamao+(s1|s2));
//if(gaman[0]>2 && gaman[1]>0 && gaman[1]<1){c=1;}
//}

//gaman[0]=g1o;
g1n=gaman[0];
g2n=gaman[1];

n=rows(X1);

an=-0.5*(g1n-1);
bn=g1n*g2n/2-1/4;

ao=-0.5*(g1o-1);
bo=g1o*g2o/2-1/4;


ind0=0;
if( ( (ao+bo)/2+(ao+bo)^2 - (an+bn)/2-(an+bn)^2)>0 ){ind0=1;}
ind1=0;
if( ( -bo/2+bo^2 +bn/2-bn^2 )>0 ){ind1=1;}

//println(ind0~ind1);

if((ind0==0 && ind1==0) || (ind0==1 && ind1==1))
{
cp=nrphi2(g1o,g2o,g1n,g2n,M_PI/2,10^(-10));
cphi=phi2(g1o,g2o,g1n,g2n,cp);
}
else{cp=M_PI/2;cphi=0;}

XX1=Intii20(X1[][:7], dtt, g1o, g2o, g1n, g2n, cp, cphi);		// println("rate ",sumc(XX1[][1]*dtt));


dAn=2*an*log(1/cos(XT/2)) + 2*bn*log(tan(XT/2)) - 2*an*log(1/cos(X0/2)) - 2*bn*log(tan(X0/2));
dAo=2*ao*log(1/cos(XT/2)) + 2*bo*log(tan(XT/2)) - 2*ao*log(1/cos(X0/2)) - 2*bo*log(tan(X0/2));

decl prn=log(dnorm(g1n,8,1))+log(dnorm(g2n,0.5,0.1^2));
decl pro=log(dnorm(g1o,8,1))+log(dnorm(g2o,0.5,0.1^2));

s0=sumc(XX1[][3]);
p=(1+exp(dAo-dAn+s0+pro-prn))^(-1);	  //println(g1o~g1n~p);
MMM=XX1[][2];

c=0;
//j=zeros(n,1);
//x=y=m=M=l1=l21=l22=mm=MM=zeros(n,1);
//l=u=zeros(n,2);					//	println("p = ",p);

while(c==0)				 ////////// AQUI!!!
{cont+=1;					   //if(fmod(cont,100)==0){println("2-coin has reached = ",cont);}

S = ranbinomial(1,1,1,beta);
if(S == 0){c=1; gama=gamao; cpp=cpo; println("coin 0");}

else
{
 C1=ranbinomial(1,1,1,p);			  //println(C1);

 if(C1==1)
 {
  c=1;gama=gaman;cpp=nrphi(g1n,g2n,M_PI/2,10^(-10));
 } //if C1=1




 else //(C1=0)
 { 	   cont1+=1;

  i=1;
  cB=1;
  while (cB==1 && i<=n)		   ////////// AQUI!!!
  {
  									  //(dtt*MMM[i-1])<0 || (dtt*MMM[i-1])>500
  np=ranpoisson(1,1,dtt*MMM[i-1]); //if( fabs(dtt*MMM[i-1])>100 || fabs(np)>100){println((dtt*MMM[i-1])~np);}
  if(np==0){i+=1;}
  else
  {
   X11=X1[i-1][:7];

      
   up=runif(np,2,zeros(np,2),(zeros(np,1)+dtt)~(zeros(np,1)+MMM[i-1]));
   up=sortbyc(up,0);


   indX2=vecindex(X2[][0],X11[0][0]);
   indX2=indX2[0]+vecindex(X2[indX2][1],X11[0][1]);
   X2i=X2[indX2][];
   								// if(X2i[0][1]==78){println("antes",X11,X2i,up[][0],dtt);}						   
   [X2i,xt]=BBx(X11,X2i,up[][0],dtt);	//if(X2i[0][1]==78){println("depois",X11,X2i,up[][0],dtt);}

   phix=phi2(g1o,g2o,g1n,g2n,xt)-XX1[i-1][0];
   if(up[][1]>phix){i+=1;}   ////////// AQUI!!!
   else{cB=0;}

   k1=rows(indX2);	 ////////// AQUI!!!
   k2=rows(X2i);
   
   X2=insertr(X2,indX2[0],k2);
   X2=dropr(X2,range((indX2[0]+k2),(indX2[0]+k2+k1-1),1));
   X2[indX2[0]:(indX2[0]+k2-1)][]=X2i;
	


   } //else np=0
   
  }	//while(cB==1)

if(cB==1){c=1;gama=gamao;cpp=cpo;}

 
 } //else (C1=0)

 } //else(S=0)
 
} //while c=0



return {X2,gama,cpp,cont,cont1};


}



