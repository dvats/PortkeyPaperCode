
xmaxim2(m,tm,X2,t,dtt,x,y)

{

decl n,Xt,w1,w2,w3,ww1,ww2,ww3,ww,i,ts,delta,r,ind,ind2,tt,xx;
decl t0,T,x0,XT,w;

w=X2[][6:];
t0=X2[][0];
T=X2[][1];
x0=X2[][2];
XT=X2[][3];

n=rows(t);
if(n==0){n=1;}

ind=0;
if(T<=tm){ind=1;}

if(ind==1)
{
ts=((tm-t0)/(tm))|((tm-t)/(tm))|((tm-T)/(tm));
ts=ivec(ts);
w1=BB2(w[1],ts[0],w[0],ts[n+1],ts[1:n]);
w2=BB2(w[3],ts[0],w[2],ts[n+1],ts[1:n]);
w3=BB2(w[5],ts[0],w[4],ts[n+1],ts[1:n]);
// w1=(x01,t0,xT1,T)  w2=(x02,t0,xT2,T)  w3=(x03,t0,xT3,T)

w1=ivec(w1);
w2=ivec(w2);
w3=ivec(w3);
ts=ivec(ts);

ww1=w[0]|w1|w[1];
ww2=w[2]|w2|w[3];
ww3=w[4]|w3|w[5];

delta=(m-x)/(sqrt(tm));
r=sqrt((delta*ts[1:n]+w1).^2+w2.^2+w3.^2);
Xt=-sqrt(tm)*r+m;
}


else
{
ts=((t0-tm)/(dtt-tm))|((t-tm)/(dtt-tm))|((T-tm)/(dtt-tm));
w1=BB2(w[0],ts[0],w[1],ts[n+1],ts[1:n]);
w2=BB2(w[2],ts[0],w[3],ts[n+1],ts[1:n]);
w3=BB2(w[4],ts[0],w[5],ts[n+1],ts[1:n]);

ww1=w[0]|w1|w[1];
ww2=w[2]|w2|w[3];
ww3=w[4]|w3|w[5];

delta=(m-y)/(sqrt(dtt-tm));
r=sqrt((delta*ts[1:n]+w1).^2+w2.^2+w3.^2);
Xt=-sqrt(dtt-tm)*r+m;

}

tt=(t0|t)~(t|T);
xx=(x0|Xt)~(Xt|XT);

ww=ww1[:n]~ww1[1:]~ww2[:n]~ww2[1:]~ww3[:n]~ww3[1:];
ind2=0;
if(ind==1)
{
if(T==tm){ind2=1;}
ww=tt~xx~(zeros(n,1)|ind2)~zeros(n+1,1)~ww;
}
else
{
if(t0==tm){ind2=1;}
ww=tt~xx~(ind2|zeros(n,1))~ones(n+1,1)~ww;
}

return {Xt,ww};
  
}	   
