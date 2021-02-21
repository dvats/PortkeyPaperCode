set.seed(10)
library(mcmcse)
library(mvtnorm)
library(truncnorm)
library(invgamma)
library(portkey)

inv.portkey <- function(prop, curr, beta = 0.99, pf, Cprop, Ccurr, ...) 
{
    x <- NA
    loops <- 0
    # print(c(Ccurr, Cprop))
    C <- Ccurr/(Ccurr + Cprop)
    if(is.na(C))
    {
    	x <- curr
    	return(list(x = x, loops = loops))
    }
    S <- 1
    while (is.na(x)) 
    {
        loops <- loops + 1
        if (beta != 1) 
            S <- rbinom(1, 1, beta)
        if (S == 0) {
            x <- curr
            return(list(x = x, loops = loops))
        }
        else {
            C1 <- rbinom(1, 1, C)
            # print(c(C, C1))
            if (C1 == 1) {
                p1 <- pf(value = curr, ...)
                C2 <- rbinom(1, 1, p1)
                if (C2 == 1) {
                  x <- prop
                  return(list(x = x, loops = loops))
                }
            }
            else {
                p2 <- pf(value = prop, ...)
                C2 <- rbinom(1, 1, p2)
                if (C2 == 1) {
                  x <- curr
                  return(list(x = x, loops = loops))
                }
            }
        }
    }
}

lims <- function(i,j,R)
{
	R0 <- Rneg1 <- R1 <- R
	R1[j,i] <- R1[i,j] <- 1
	Rneg1[j,i] <- Rneg1[i,j] <- -1
	R0[j,i] <- R0[i,j] <- 0
	lims <- c(-1,1)

	f1 <- det(R1)
	fneg1 <- det(Rneg1)
	f0 <- det(R0)
	a <- (f1 + fneg1 - 2*f0)/2
	b <- (f1 - fneg1)/2
	c <- f0	
  if(delta(a,b,c) > 0){ # first case D>0
        x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
        x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
        result = c(x_1,x_2)
        return(result)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
        x = -b/(2*a)
        return(x)
  }
  else {return(c(-1,1))} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
      b^2-4*a*c
}

rij_like <- function(i,j,y, R, mu, sig2)
{
	n <- dim(y)[1]
	p <- dim(y)[2]

	B <- t(y) %*% y
	rij <- R[i,j]
	foo <- -(n/2)*determinant(R, logarithm = TRUE)$modulus - sum(diag((solve(R) %*%B)))/2
	foo <- foo - (rij - mu)^2/(2*sig2)
	return(foo)
}

drawR1 <- function(value, sig2, p)
{
	l <- p*(p-1)/2
	rs <- rtruncnorm(l, a = -1, b = 1, mean = value, sd = sqrt(sig2) ) #  rnorm(l, mean = value, sd = sqrt(sig2))
	R <- matrix(0, nrow = p, ncol = p)
	R[which(lower.tri(R))] <- rs
	R <- R + t(R)
	diag(R) <- 1
	rtn <- ifelse(min(eigen(R)$values) <= 0, 0, 1)
	return(rtn)
}
drawR2 <- function(value, mu, p)
{
	l <- p*(p-1)/2
	rs <- rtruncnorm(l, a = -1, b = 1, mean = mu, sd = sqrt(value) ) #rnorm(l, mean = mu, sd = sqrt(value))
	R <- matrix(0, nrow = p, ncol = p)
	R[which(lower.tri(R))] <- rs
	R <- R + t(R)
	diag(R) <- 1
	rtn <- ifelse(min(eigen(R)$values) <= 0, 0, 1)
	return(rtn)
}

Cbound.mu <- function(mu, mu.mean, mu.var, sig2, l)
{
	foo <- pnorm((1 - mu)/sqrt(sig2)) - pnorm((-1 - mu)/sqrt(sig2))^l
	 rtn <- log(foo) - dnorm(mu, mean = mu.mean, sd = sqrt(mu.var), log = TRUE)
	 return(exp(rtn))
}

Cbound.sig <- function(sig2, sig2.shape, sig2.param, mu, l)
{
	foo <- pnorm((1 - mu)/sqrt(sig2)) - pnorm((-1 - mu)/sqrt(sig2))^l
	 rtn <- log(foo) - dinvgamma(sig2, shape = sig2.shape, rate = sig2.param, log = TRUE)
	 return(exp(rtn))
}
CorMCMC <- function(y, tau2, alpha, beta, port = 1, nsim = 1e3, h = rep(.005, dim(y)[2]*(dim(y)[2] - 1)/2))
{
	n <- dim(y)[1]
	p <- dim(y)[2]
	R <- cor(y)#matrix(.9, nrow = p, ncol = p)
	diag(R) <- 1
	l <- p*(p-1)/2	

	mu <- sum(lower.tri(R)*R)/l
	sig2 <- 1
	loops <- matrix(0, nrow = nsim, ncol = 2)


	acc <- rep(0, l+2)
	out <- matrix(0, nrow = nsim, ncol = l+2)
	out[1, ] <- c(R[which(lower.tri(R) )], mu, sig2)	

	lower <- lower.tri(R)

	for(t in 2:nsim)
	{
		if(t%%(nsim/10) == 0) print(paste("nsim =", t))
		count <- 0
		for(i in 2:p)
		{
			for(j in 1:(i-1))
			{	
				count <- count+1
				R.prop <- R
				lim <- lims(i,j,R)
				# print(c(t,lim))
				prop <- R[i,j] + runif(1, min = -h[count], max = h[count]) #rtruncnorm(1, a=lim[1], b=lim[2], mean = R[i,j], sd = h)#runif(1, min = lim[1], max = lim[2])
				R.prop[i,j] <- R.prop[j,i] <- prop

				if( prop < lim[2] && prop > lim[1])
				{
					ratio <- rij_like(i,j, y, R.prop, mu, sig2) -  rij_like(i,j, y, R, mu, sig2)
					if(log(runif(1)) < ratio)
					{
						R[i,j] <- R[j,i] <- prop
						acc[count] <- acc[count] + 1
					}
				} 
			}
		}
		r.bar <- sum(lower*R)/l
		mu.mean <- tau2*l*r.bar/(sig2 + l*tau2) 
		mu.var <- sig2*tau2/(sig2 + l*tau2)
		

		# portkey
		mu.prop <- rnorm(1, mean = mu, sd = .55)
		Cprop <- Cbound.mu(mu = mu.prop, mu.mean, mu.var, sig2 = sig2, l = l)
		Ccurr <- Cbound.mu(mu = mu, mu.mean, mu.var, sig2 = sig2, l = l)
		# mu.prop <- rnorm(1, mean = mu.mean, sd = sqrt(mu.var))

		port1 <- inv.portkey(prop = mu.prop, curr = mu, beta = port, pf = drawR1, Cprop = Cprop, Ccurr = Ccurr, sig2 = sig2, p = p)
		if(port1$x == mu.prop) acc[count+1] = acc[count+1] + 1
		mu <- port1$x
		loops[t, 1] <- port1$loops

		sig2.param <- (2*beta + sum((R[which(lower )] - mu)^2)  )/2
		# sig2.prop <- 1/rgamma(1, shape = alpha + l/2, rate = sig2.param)

		sig2.prop <- rnorm(1, mean = sig2, sd = .40)
		
		if(sig2.prop > 0)
		{
			Cprop <- Cbound.sig(sig2.prop, sig2.shape = alpha + l/2, sig2.param = sig2.param, mu = mu, l = l)
			Ccurr <- Cbound.sig(sig2, sig2.shape = alpha + l/2, sig2.param = sig2.param, mu = mu, l = l)

			port2 <- inv.portkey(prop = sig2.prop, curr = sig2, beta = port, pf = drawR2, Cprop = Cprop, Ccurr = Ccurr, mu = mu, p = p)
			if(port2$x == sig2.prop) acc[count+2] = acc[count+2] + 1
			sig2 <- port2$x
			loops[t, 2] <- port2$loops		
		}
		out[t, ] <- c(R[which(lower )], mu, sig2)
	}
	print(colMeans(loops))
	print(apply(loops,2, max))
	print(acc/nsim)
	return(list("out" = out, "loops" = loops, "acc" = tail(acc,2)/nsim) )
}


y <- EuStockMarkets
y <- scale(y, center = FALSE, scale = apply(y, 2, sd))
tau2 <- 1
alpha <- 1
beta <- 1
nsim <- 1e4
h <- c(.0015, rep(.005,4), .01)
R <- cor(y)

# time1 <- system.time(chain1 <- CorMCMC(y = y, tau2 = tau2, alpha = alpha, beta = beta, port = 1, nsim = nsim, h = h))
# time9 <- system.time(chain9 <- CorMCMC(y = y, tau2 = tau2, alpha = alpha, beta = beta, port = .90, nsim = nsim, h = h))


B <- 1e1
time <- matrix(0, nrow = B, ncol = 2)
ess <- matrix(0, nrow = B, ncol = 2)
max1.l <- matrix(0, nrow = B, ncol = 2)
mean1.l <- matrix(0, nrow = B, ncol = 2)

max2.l <- matrix(0, nrow = B, ncol = 2)
mean2.l <- matrix(0, nrow = B, ncol = 2)

for(b in 1:B)
{

	time[b,1] <- system.time(chain1 <- CorMCMC(y = y, tau2 = tau2, alpha = alpha, beta = beta, port = 1, nsim = nsim, h = h))[3]
	time[b,2] <- system.time(chain99 <- CorMCMC(y = y, tau2 = tau2, alpha = alpha, beta = beta, port = .90, nsim = nsim, h = h))[3]


	ess[b,1] <- multiESS(chain1$out, size = "sqroot", r = 1)
	ess[b,2] <- multiESS(chain99$out, size = "sqroot", r = 1)

	max1.l[b,1] <- max(chain1$loops[,1])
	max1.l[b,2] <- max(chain99$loops[,1])

	max2.l[b,1] <- max(chain1$loops[,2])
	max2.l[b,2] <- max(chain99$loops[,2])	

	mean1.l[b,1] <- mean(chain1$loops[,1])
	mean1.l[b,2] <- mean(chain99$loops[,1])

	mean2.l[b,1] <- mean(chain1$loops[,2])
	mean2.l[b,2] <- mean(chain99$loops[,2])	

	print(b)

	save(time, ess, max1.l, max2.l, mean1.l, mean2.l, file = "correlation_sim")
}

rm(list = ls())
load("correlation_sim")
B <- 1e1
round(colMeans(ess))
round(apply(ess, 2, sd)/sqrt(B), 1)

round(colMeans(ess/time), 3)
round(apply(ess/time, 2, sd)/sqrt(B), 3)

round(colMeans(mean1.l), 2)
round(apply(mean1.l, 2, sd)/sqrt(B), 2)

round(colMeans(mean2.l), 2)
round(apply(mean2.l, 2, sd)/sqrt(B), 2)

round(colMeans(max1.l))
round(apply(max1.l, 2, sd)/sqrt(B), 2)

round(colMeans(max2.l))
round(apply(max2.l, 2, sd)/sqrt(B), 2)


