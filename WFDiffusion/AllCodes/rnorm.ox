

rnorm(const a,const b, const mu,const sigma2)

 {
  decl u;

  u = rann(a,b).*sqrt(sigma2)+mu;
  return u;
 }
