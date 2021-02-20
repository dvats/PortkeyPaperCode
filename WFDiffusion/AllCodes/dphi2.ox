
dphi2(const g1n, const g2n, const g1o, const g2o, const x)

{
  decl f,c1n,c2n,c3n,c4n,an,bn,c1o,c2o,c3o,c4o,ao,bo,c1,c2,c3,c4;

an=-0.5*(g1n-1);
bn=g1n*g2n/2-1/4;

c1n=an^2;
c2n=bn^2;
c3n=(an+bn)/2 + 2*an*bn;
c4n=-bn/2;

ao=-0.5*(g1o-1);
bo=g1o*g2o/2-1/4;

c1o=ao^2;
c2o=bo^2;
c3o=(ao+bo)/2 + 2*ao*bo;
c4o=-bo/2;

c1=c1n-c1o;
c2=c2n-c2o;
c3=c3n-c3o;
c4=c4n-c4o;
  
f= c1*sin(x/2)/(cos(x/2))^3 + c2*((sin(x/2))^2-(cos(x/2))^2)/((cos(x/2))^3*(sin(x/2))^3) + c3*sin(x/2)/(cos(x/2))^3 + c4*cos(x/2)/(sin(x/2))^3;

return f;


}