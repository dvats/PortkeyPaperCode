
dphi(const g1, const g2, const x)

{
  decl a,b,c1,c2,c3,c4,f;

a=(1-g1)/2;
b=g1*g2/2-1/4;

c1=a^2;
c2=b^2;
c3=(a+b)/2 + 2*a*b;
c4=b/2;
  
f= c1*sin(x/2)/(cos(x/2))^3 + c2*((sin(x/2))^2-(cos(x/2))^2)/((cos(x/2))^3*(sin(x/2))^3) + c3*sin(x/2)/(cos(x/2))^3 + c4*cos(x/2)/(sin(x/2))^3;

return f;


}