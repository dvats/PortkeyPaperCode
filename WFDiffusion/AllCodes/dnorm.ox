#include <oxstd.h>
#include <oxfloat.h>

dnorm(const x, const mu,const sigma2)

 {
  decl u;

  u = ((sqrt(M_2PI*sigma2)).^(-1)).*exp(-((2*sigma2).^(-1)).*((x-mu).^2));

  return u;
 }