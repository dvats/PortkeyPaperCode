
phi(const g1, const g2, const x)

{
  decl f,c1,c2,c3,c4,a,b;

a=-0.5*(g1-1);
b=g1*g2/2-1/4;

c1=a^2/2;
c2=b^2/2;
c3=(a+b)/4 + a*b;
c4=-b/4;

  
f= c1*(sin(x/2)).^2 ./(cos(x/2)).^2 + c2./( (sin(x/2)).^2 .* (cos(x/2)).^2 ) + c3./(cos(x/2)).^2 + c4./(sin(x/2)).^2;

return f;


}