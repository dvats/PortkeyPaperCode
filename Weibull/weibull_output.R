load("weibull_one_run")

pdf("wei_trace.pdf", height = 5, width = 7)
par(mfrow = c(2,2))
par(mar=c(3, 2, 1.5, 2) + 0.1)
plot.ts(tail(bark$mcmc, 1e4), main = expression(paste(beta, " = 1")) , ylab = "X", xlab = "", ylim = c(0, .26))
text(5000, .25, bark$accept, col = "red")
plot.ts(tail(stable_99$mcmc, 1e4), main = expression(paste(beta, " = .99")) ,ylab = "", xlab = "",  ylim = c(0, .26))
text(5000, .25, stable_99$accept, col = "red")
par(mar=c(3, 2, 1.5, 2) + 0.1)
plot.ts(tail(stable_90$mcmc, 1e4), main = expression(paste(beta, " = .90")) , ylab = "",  ylim = c(0, .26))
text(5000, .25, stable_90$accept, col = "red")
plot.ts(tail(stable_75$mcmc, 1e4), main = expression(paste(beta, " = .75")) ,ylab = "",  ylim = c(0, .26))
text(5000, .25, stable_75$accept, col = "red")
dev.off()

acf1 <- acf(bark$mcmc, lag.max = 50, plot = FALSE)$acf
acf99 <- acf(stable_99$mcmc, lag.max = 50, plot = FALSE)$acf
acf90 <- acf(stable_90$mcmc, lag.max = 50, plot = FALSE)$acf
acf75 <- acf(stable_75$mcmc, lag.max = 50, plot = FALSE)$acf

pdf("wei_acf.pdf", height = 4, width = 6)
plot(0:50, acf1, type= 'l', xlab = "Lags", ylab = "Autocorrelation")
lines(0:50, acf99, lty = 2)
lines(0:50, acf90, lty = 3)
lines(0:50, acf75, lty = 4)
legend("topright", title = expression(beta), legend = c("1", ".99", ".90", ".75"), lty = 1:4, col = "black")
dev.off()



load("wei_sim")

round(colMeans(ess))
round(apply(ess, 2, sd)/sqrt(B), 2)

round(colMeans(ess/time), 2)
round(apply(ess/time, 2, sd)/sqrt(B), 2)

round(colMeans(mean.l), 2)
round(apply(mean.l, 2, sd)/sqrt(B), 2)

round(colMeans(max.l))
round(apply(max.l, 2, sd)/sqrt(B), 2)