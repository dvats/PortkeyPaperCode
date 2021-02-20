
rmvnorm(const d, const mu,const sigma2)

 {
  decl u,A;

  A=choleski(sigma2);

  u = A*rann(d,1)+mu;
  
  return u;
 }
