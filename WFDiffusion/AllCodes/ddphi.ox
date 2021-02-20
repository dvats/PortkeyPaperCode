
ddphi(const g1, const g2, const x)

{
  decl a,b,c1,c2,c3,f;

a=(1-g1)/2;
b=g1*g2/2-1/4;

c1=a^2/2+(a+b)/4 + a*b;
c2=b^2/2;
c3=-b/2;

f= c1*(1/(cos(x/2))^2 +3*(sin(x/2))^2/(cos(x/2))^4  ) + c2*( -2/((sin(x/2))^2*(cos(x/2))^2) + 3/(cos(x/2))^4 + 3/(sin(x/2))^4 ) + c3*( 1/(sin(x/2))^2 -3*(cos(x/2))^2/(sin(x/2))^4);

return f;


}