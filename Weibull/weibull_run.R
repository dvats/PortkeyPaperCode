
library(mcmcse)
library(ellipse)
library(coda)
set.seed(100)

cbound <- function(x, k)
{
	return(k/(exp(1)*x))
}

b2coin <- function(x.curr, x.prop, k, beta, a1, b1)
{
	foo <- 0

	c.prop <- cbound(x = x.prop, k = k)
	c.curr <- cbound(x = x.curr, k = k)
	loops <- 0

	while(foo == 0)
	{
		loops <- loops+1
		S <- rbinom(1, 1, prob = beta)
		if(S == 0)
		{
			return(c(x.curr, loops))
		} else {

			C1 <- rbinom(1, 1, c.prop/(c.curr + c.prop))
			if(C1 == 1)
			{
				lam1 <- rgamma(1, shape = a1, rate = b1)
				p1 <- dweibull(x.prop, shape = k, scale = lam1)/c.prop

				C2 <- rbinom(1, 1, p1)
				if(C2 == 1)
				{
					return(c(x.prop, loops))
				} 
			} else {
				lam1 <- rgamma(1, shape = a1, rate = b1)
				p2 <- dweibull(x.curr, shape = k, scale = lam1)/c.curr

				C2 <- rbinom(1, 1, p2)
				if(C2 == 1)
				{
					return(c(x.curr, loops))
				}
			}
		}
	}
}

mcmc <- function(N = 1e3, beta = .95, start, h = 4, k = 5, a1 = 2, b1 = 5)
{
	out <- numeric(length = N)
	loops <- numeric(length = N)
	out[1] <- start
	acc <- 0
	for(i in 2:N)
	{
		prop <- rnorm(1, out[i-1], sd = sqrt(h))
		# .3*rgamma(1, shape = a1, rate = b1) + .7*rgamma(1, shape = a2, rate = a2)

		if(prop < 0) 
		{
			out[i] <- out[i-1]
			next
		}

		interim <- b2coin(x.curr = out[i-1], x.prop = prop, k = k, beta = beta, a1 = a1, b1 = b1)
		out[i] <-  interim[1]
		if(out[i] != out[i-1]) acc <- acc+1
		loops[i] <- interim[2]
	}
	return(list("mcmc" = out, "loops" = loops, "accept" = acc/N))
}


a1 <- 10
b1 <- 100
k <- 10

true.mean <- (a1/b1 * gamma(1 + 1/k)) 
foo <- gamma(1 + 2/k) - (gamma(1 + 1/k))^2
true.var <-  ( (gamma(1 + 1/k))^2* a1/b1^2 + foo*(a1/b1^2 + a1^2/b1^2)) 

bark <- mcmc(N = 1e5, start = true.mean, beta = 1,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_99 <- mcmc(N = 1e5, start = true.mean, beta = .99,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_90 <- mcmc(N = 1e5, start = true.mean, beta = .90,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_75 <- mcmc(N = 1e5, start = true.mean, beta = .75,  k = k, a1 = a1, b1 = b1, h = true.var)

save(bark, stable_99, stable_90, stable_75, file = "weibull_one_run")



############### ESS and multiple reps ############
a1 <- 10
b1 <- 100
k <- 10

B <- 1000
time <- matrix(0, nrow = B, ncol = 4)
ess <- matrix(0, nrow = B, ncol = 4)
max.l <- matrix(0, nrow = B, ncol = 4)
mean.l <- matrix(0, nrow = B, ncol = 4)

for(b in 1:B)
{
	N <- 1e5

	time[b,1] <- system.time(bark <- mcmc(N = N, start = true.mean, beta = 1,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
	time[b,2] <- system.time(stable_99 <- mcmc(N = N, start = true.mean, beta = .99,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
	time[b,3] <- system.time(stable_90 <- mcmc(N = N, start = true.mean, beta = .90,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
	time[b,4] <- system.time(stable_75 <- mcmc(N = N, start = true.mean, beta = .75,  k = k, a1 = a1, b1 = b1, h = true.var))[3]

	ess[b,1] <- effectiveSize(bark$mcmc)
	ess[b,2] <- effectiveSize(stable_99$mcmc)
	ess[b,3] <- effectiveSize(stable_90$mcmc)
	ess[b,4] <- effectiveSize(stable_75$mcmc)


	max.l[b,1] <- max(bark$loops)
	max.l[b,2] <- max(stable_99$loops)
	max.l[b,3] <- max(stable_90$loops)
	max.l[b,4] <- max(stable_75$loops)

	mean.l[b,1] <- mean(bark$loops)
	mean.l[b,2] <- mean(stable_99$loops)
	mean.l[b,3] <- mean(stable_90$loops)
	mean.l[b,4] <- mean(stable_75$loops)

	print(b)

	save(time, ess, max.l, mean.l, file = "wei_sim")
}






