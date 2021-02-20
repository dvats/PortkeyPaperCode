
deltaK(m,t,Xt,XM)

{

decl u1,u,c1,K,i,j,S,Sj1,Sj2,DI,Kh,zeta,xi;


K=fabs(XM-m);

u1=fabs(Xt-m);

u=ranu(1,1);

c1=0;

Kh=ceil(sqrt(t+K^2)/(2*K));


if((3*K^2)>t || Kh==1)
{

j=1;
S=zeros(1,1);
while(c1==0)
{


   zeta=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xi=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2=1-(1/u1)*(S+zeta);
   S+=zeta-xi;
   Sj1=1-(1/u1)*S;

 
 if(u<Sj2){DI=1;c1=1;}
 else if(u>Sj1){DI=0;c1=1;}
 else{j+=1;}
}

}



else
{

S=zeros(1,1);
for(j=1;j<=(Kh-1);++j)
{
   zeta=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xi=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2=1-(1/u1)*(S+zeta);
   S+=zeta-xi;
   Sj1=1-(1/u1)*S;
}

j=Kh;
while(c1==0)
{

   zeta=(2*K*j-u1)*exp(-2*K*j*(K*j-u1)/t);
   xi=(2*K*j+u1)*exp(-2*K*j*(K*j+u1)/t);
   Sj2=1-(1/u1)*(S+zeta);
   S+=zeta-xi;
   Sj1=1-(1/u1)*S;

 if(u<Sj2){DI=1;c1=1;}
 else if(u>Sj1){DI=0;c1=1;}
 else{j+=1;}
}


}


return(DI);
  
}	   
  