
deltaL(m,t,Xt,XM1,XM2)

{

decl u1,u,c1,K,L,i,j,SK,SL,Sj1K,Sj2K,Sj1L,Sj2L,DI,Kh,zetaK,xiK,zetaL,xiL;


K=fabs(XM1-m);

L=fabs(XM2-m);

u1=fabs(Xt-m);

u=ranu(1,1);

c1=0;

Kh=ceil(sqrt(t+K^2)/(2*K));


if(3*K^2>t || Kh==1)
{

j=1;
SK=SL=zeros(1,1);
while(c1==0)
{


   zetaK=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xiK=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2K=u1-(SK+zetaK);
   SK+=zetaK-xiK;
   Sj1K=u1-SK;

   zetaL=(2*L*j-u1)*exp(-2*L*j*(L*j-u1)/t);
   xiL=(2*L*j+u1)*exp(-2*L*j*(L*j+u1)/t);
   Sj2L=u1-(SL+zetaL);
   SL+=zetaL-xiL;
   Sj1L=u1-SL;

 
 if(u<(Sj2K/Sj1L)){DI=0;c1=1;}
 else if(u>(Sj1K/Sj2L)){DI=1;c1=1;}
 else{j+=1;}
}

}



else
{

SK=SL=zeros(1,1);
for(j=1;j<=(Kh-1);++j)
{
   zetaK=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xiK=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2K=u1-(SK+zetaK);
   SK+=zetaK-xiK;
   Sj1K=u1-SK;

   zetaL=(2*L*j-u1)*exp(-2*L*j*(L*j-u1)/t);
   xiL=(2*L*j+u1)*exp(-2*L*j*(L*j+u1)/t);
   Sj2L=u1-(SL+zetaL);
   SL+=zetaL-xiL;
   Sj1L=u1-SL;
}

c1=0;
j=Kh;
while(c1==0)
{

   zetaK=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xiK=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2K=u1-(SK+zetaK);
   SK+=zetaK-xiK;
   Sj1K=u1-SK;

   zetaL=(2*L*j-u1)*exp(-2*L*j*(L*j-u1)/t);
   xiL=(2*L*j+u1)*exp(-2*L*j*(L*j+u1)/t);
   Sj2L=u1-(SL+zetaL);
   SL+=zetaL-xiL;
   Sj1L=u1-SL;

 if(u<(Sj2K/Sj1L)){DI=0;c1=1;}
 else if(u>(Sj1K/Sj2L)){DI=1;c1=1;}
 else{j+=1;}
}


}


return(DI);
  
}	   
  