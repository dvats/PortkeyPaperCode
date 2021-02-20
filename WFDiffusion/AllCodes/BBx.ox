
BBx(X1,X2,t,dtt)

{

decl x,y,indm,tm,m,L,n,i1,i2,i,c,xtm,c1,xtM,xt,x2,tt,ind,nn,cont,ind2,indmM,sindmM,j,x22,xttt,temp,c2,ttt,xm,xM,ww,xxx;

x=X1[][2];
y=X1[][3];

indm=X1[][4];
tm=X1[][5];
m=X1[][6];
L=X1[][7];


xm=indm*m+(1-indm)*L;
xM=indm*L+(1-indm)*m;

xttt=zeros(1,1);

n=rows(t);

if(n==0){n=1;}

if(rows(X2)==2)
{

x2=X2;

// if(X2[][6]==1)
// {
  i1=i2=0;
  for(i=1;i<=n;++i)
  {
   if(t[i-1]<tm){i1+=1;}
   else{i2+=1;}
  }

  if(i1>0)
  {
   c=0;
   while(c==0)			
   {
  	if(indm==1){c2=0;[xtm,ww]=xminim2(m,tm,X2[0][2:],t[:(i1-1)],dtt,x,y);if(xtm<L){c2=1;}}
	else{c2=0;[xtm,ww]=xmaxim2(m,tm,X2[0][2:],t[:(i1-1)],dtt,x,y);if(xtm>L){c2=1;}}
	if(c2==1)
	{
	 c=deltaK(m,tm-t[i1-1],xtm[i1-1],L);
	 i=1;
	 ttt=0|t[:(i1-1)]; xxx=x|xtm;  // size i1+1
	 while(c==1 && i<=i1){c=deltaK1(ttt[i1-i],xxx[i1-i],ttt[i1-i+1],xxx[i1-i+1],m,L);i+=1;}
	}
   }
   x2=((X2[0][0]+zeros(1+i1,1))~(X2[0][1]+zeros(1+i1,1))~ww)|(x2[1][]);
   xttt=xttt|xtm;
  }

  if(i2>0)
  {
   c=0;
   while(c==0)						  
   {
   	if(indm==1){c2=0;[xtM,ww]=xminim2(m,tm,X2[1][2:],t[i1:],dtt,x,y);if(xtM<L){c2=1;}}
	else{c2=0;[xtM,ww]=xmaxim2(m,tm,X2[1][2:],t[i1:],dtt,x,y);if(xtM>L){c2=1;}}
	if(c2==1)
	{
	 c=deltaK(m,t[i1]-tm,xtM[0],L);
	 i=1;
	 ttt=t[i1:]|dtt; xxx=xtM|y;  // size i2+1
	 while(c==1 && i<=i2){c=deltaK1(ttt[i-1],xxx[i-1],ttt[i],xxx[i],m,L);i+=1;}
	}
   }
    x2=x2[:(rows(x2)-2)][]|((X2[0][0]+zeros(1+i2,1))~(X2[0][1]+zeros(1+i2,1))~ww);
    xttt=xttt|xtM;
  }

}



else   //rows(X2)=1
{
 x2=zeros(1,14);

 tt=sortbyc((t|X2[][2]|dtt)~(ones(n,1)|zeros(rows(X2)+1,1)),0);
 ind=vecindex(1-tt[][1]);
 nn=rows(ind)-1;
 cont=zeros(nn,1);
 for(i=1;i<=nn;++i)
 {
  cont[i-1]=ind[i]-ind[i-1]-1;
 }
 ind=vecindex(cont); //row indexes of X2 where there are points to be simulated
 ind2=vecindex(cont,0); //row indexes of X2 where there are no points to be simulated
 cont=cont[ind];     //number of points to be simulated in each of these
// indmM=vecindex(indmM);

 cont=0|cont;
 
 for(i=1;i<=rows(ind);++i)
 {

  tt=t[(sumc(cont[:(i-1)])):(sumc(cont[:i])-1)];
  
  i1=i2=0;
  n=cont[i];
  for(j=1;j<=n;++j)
  {
   if(tt[j-1]<tm){i1+=1;}
   else{i2+=1;}
  }

  if(i1>0)
  {
   c=0;
   while(c==0)
   {
    if(indm==1){c2=0;[xtm,ww]=xminim2(m,tm,X2[ind[i-1]][2:],tt[:(i1-1)],dtt,x,y);if(xtm<L){c2=1;}}
	else{c2=0;[xtm,ww]=xmaxim2(m,tm,X2[ind[i-1]][2:],tt[:(i1-1)],dtt,x,y);if(xtm>L){c2=1;}}
	if(c2==1)
	{
	 if(X2[ind[i-1]][6]==1){c=deltaK(m,tm-tt[i1-1],xtm[i1-1],L);}
	 else{c=deltaK1(tt[i1-1],xtm[i1-1],X2[ind[i-1]][3],X2[ind[i-1]][5],m,L);}
	 j=1;
	 ttt=X2[ind[i-1]][2]|tt[:(i1-1)]; xxx=X2[ind[i-1]][4]|xtm;  // size i1+1
	 while(c==1 && j<=i1){c=deltaK1(ttt[i1-j],xxx[i1-j],ttt[i1-j+1],xxx[i1-j+1],m,L);j+=1;}
	}
   }
   x2=x2|((X2[0][0]+zeros(1+i1,1))~(X2[0][1]+zeros(1+i1,1))~ww);
   xttt=xttt|xtm;
  }

  if(i2>0)
  {
   c=0;
   while(c==0)
   {
    if(indm==1){c2=0;[xtM,ww]=xminim2(m,tm,X2[ind[i-1]][2:],tt[i1:],dtt,x,y);if(xtM<L){c2=1;}}
	else{c2=0;[xtM,ww]=xmaxim2(m,tm,X2[ind[i-1]][2:],tt[i1:],dtt,x,y);if(xtM>L){c2=1;}}
	if(c2==1)
	{
	 if(X2[ind[i-1]][6]==1){c=deltaK(m,tt[i1]-tm,xtM[0],L);}
	 else{c=deltaK1(X2[ind[i-1]][2],X2[ind[i-1]][4],tt[i1],xtM[0],m,L);}
	 j=1;
	 ttt=tt[i1:]|X2[ind[i-1]][3]; xxx=xtM|X2[ind[i-1]][5];  // size i1+1
	 while(c==1 && j<=i2){c=deltaK1(ttt[j-1],xxx[j-1],ttt[j],xxx[j],m,L);j+=1;}
	}
   }
    x2=x2|((X2[0][0]+zeros(1+i2,1))~(X2[0][1]+zeros(1+i2,1))~ww);
	xttt=xttt|xtM;
  }

 }	 //for

x2=x2|X2[ind2][];
x2=x2[1:][];
x2=sortbyc(x2,2);

}	//else rows(X2)=1


xttt=xttt[1:];

return {x2,xttt};

}	   
  