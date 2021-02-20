
condmaxim(t0,x0,T,XT,a,b)

{

decl mi,tm1,m1,tm,m,mm;

mi=condminim(t0,x0,T,XT,x0+XT-b,x0+XT-a);  

tm1=mi[0];
m1=mi[1];

tm=T-(tm1-t0);
m=max(x0,XT)+(min(x0,XT)-m1);

mm=tm~m;

return mm;
  
}	   
  