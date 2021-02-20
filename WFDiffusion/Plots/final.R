################################################
###### Making the plots for the paper
################################################
library(mcmcse)
# Gamma_1 for Barker's
g1_old <- t(t(read.table("output/cont2a_test4.mat"))) 
g2_old <- t(t(read.table("output/cont2b_test4.mat") ))

# Gamma_2 for portkey Barker's
g1_new <- t(t(read.table("output/cont2a_new4.mat") ))
g2_new <- t(t(read.table("output/cont2b_new4.mat") ))


pdf("WF_loops.pdf", height = 5, width = 5)
plot(g1_old, pch = 3, xlab = "Iteration", ylab = "Number of Bernoulli factory loops")
points(g1_new, col = "red", pch = 20)
legend("topleft", legend = c("two-coin", "Portkey two-coin"), pch = c(3,20), col = c("black", "red"))
dev.off()

pdf("WF_loops2.pdf", height = 5, width = 5)
plot(g2_old, pch = 3, xlab = "Iteration", ylab = "Number of Bernoulli factory loops")
points(g2_new, col = "red", pch = 20)
legend("topleft", legend = c("two-coin", "Portkey two-coin"), pch = c(3,20), col = c("black", "red"))
dev.off()


round(c(mean(g1_old), max(g1_old)))
round(c(mean(g1_new), max(g1_new)) )

round(c(mean(g2_old), max(g2_old)) )
round(c(mean(g2_new), max(g2_new) ) )


##########################################


ga_old <- read.table("output/gama_test4.mat")
ga_old <- t(t(ga_old))

ga_new <- read.table("output/gama_new4.mat")
ga_new <- t(t(ga_new))

pdf("WF_gamma1.pdf", height = 5, width = 5)
plot(density(ga_old[,1]), xlab = expression(gamma[1]), main = "", lty = 2)
lines(density(ga_new[,1]), col = "red")
legend("topleft", legend = c("Barker's", "Portkey Barker's"), lty = c(2,1), col = c("black", "red"))
dev.off()

pdf("WF_gamma2.pdf", height = 5, width = 5)
plot(density(ga_old[,2]), xlab = expression(gamma[2]), main = "", lty = 2, ylim = range(density(ga_new[,2])$y))
lines(density(ga_new[,2]), col = "red")
legend("topleft", legend = c("Barker's", "Portkey Barker's"), lty = c(2,1), col = c("black", "red"))
dev.off()

pdf("WF_acf.pdf", height = 5, width = 5)
par(mfrow = c(1,1))
plot(acf(ga_old[,1], main = expression(gamma[1]), plot = FALSE)$acf, type = 'l', ylab = expression(paste("ACF for ", gamma[1])))
lines(acf(ga_new[,1], main = expression(gamma[1]), plot = FALSE)$acf, lty = 2)
legend("topright", title = expression(beta), legend = c("1", ".99"), lty = 1:2, col = "black")
dev.off()

pdf("WF_acf2.pdf", height = 5, width = 5)
par(mfrow = c(1,1))
plot(acf(ga_old[,2], main = expression(gamma[2]), plot = FALSE)$acf, type = 'l', ylab = expression(paste("ACF for ", gamma[2])))
lines(acf(ga_new[,2], main = expression(gamma[2]), plot = FALSE)$acf, lty = 2)
legend("topright", title = expression(beta), legend = c("1", ".99"), lty = 1:2, col = "black")
dev.off()


## Effective sample size
multiESS(ga_old)
multiESS(ga_new)




# Running times from machine clock
#times

#34:35:18
t_old=34*60*60+35*60+18

#3:29:50
t_new=3*60*60+29*60+50

# t_old/ess1_old
# t_old/ess2_old

# t_new/ess1_new
# t_new/ess2_new

## ESSs
multiESS(ga_old, method = "lug")
multiESS(ga_new, method = "lug")

# ESS per hour
multiESS(ga_old, method = "lug")/ (t_old/3600)
multiESS(ga_new, method = "lug")/(t_new/3600)


### Plot data
dat <- read.table("data.txt")
pdf("data.pdf", height = 5, width = 7)
plot(0:50,dat$V1, type= "b", xlab = "t", ylab = "")
dev.off()

# Settings for the parameters for memory

#5 ->   g1 = 0.99995 , g2 = 0.9995 , X = 0.98

# old - a.r.    g1= 0.48378     g2= 0.46236
#average # X-loops = 1.69856

# old - a.r.    g1= 0.48914     g2= 0.38924
#average # X-loops = 1.64597


