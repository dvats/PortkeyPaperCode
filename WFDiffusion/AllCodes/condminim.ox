
condminim(t0,x0,T,XT,a,b)

{

decl E,Z1,Z2,c1,c2,u,I1,I2,V,m,tm,mm,uu,pa,pb;

pa=exp(((XT-x0)^2-(XT+x0-2*a)^2)/(2*(T-t0)));  
pb=exp(((XT-x0)^2-(XT+x0-2*b)^2)/(2*(T-t0)));

uu=runif(1,1,pa,pb);

E=-log(uu);
Z1=((XT-x0)-sqrt(2*(T-t0)*E+(XT-x0)^2))/2;
c1=(((XT-x0)-Z1)^2)/(2*(T-t0));
c2=(Z1^2)/(2*(T-t0)); 

u=ranu(1,1);
I1=igaus(sqrt(c1/c2),2*c1);	
I2=1/igaus(sqrt(c2/c1),2*c2);

if(u<((1+sqrt(c1/c2))^(-1))){V=I1;}
else{V=I2;}

Z2=(T-t0)/(1+V);

m=Z1+x0;
tm=t0+Z2;

mm=tm~m;

return mm;
  
}	   
  