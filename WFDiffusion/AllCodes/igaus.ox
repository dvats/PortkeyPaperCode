
igaus(a,b)

{

decl z,y,X1,X,u;

z=rann(1,1);
y=z^2;

X1=a+((a^2*y)/(2*b))-(a/(2*b))*sqrt(4*a*b*y+(a^2)*(y^2));

u=ranu(1,1);

if(u<=(a/(a+X1))){X=X1;}
else{X=(a^2)/(X1);}

return X;


}	   
  