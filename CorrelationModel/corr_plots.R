
# load("one_run_acf")
load("compare_shadow")

a7.1 <- acf(chain1$out[,7], plot = FALSE)$acf
a8.1 <- acf(chain1$out[,8], plot = FALSE)$acf
a7.9 <- acf(chain9$out[,7], plot = FALSE)$acf
a8.9 <- acf(chain9$out[,8], plot = FALSE)$acf


pdf("corr_acf1.pdf", height = 5, width = 5)
plot(a7.1, col = "black", type = 'l', ylab = expression(paste("ACF for ", mu)) )
lines(a7.9, col = "black", lty = 2)
legend("topright", title = expression(beta), legend = c("1", ".90"), lty = 1:2, col = "black")
dev.off()

pdf("corr_acf2.pdf", height = 5, width = 5)
plot(a8.1, col = "black", type = 'l', ylab = expression(paste("ACF for ", sigma^2)))
lines(a8.9, col = "black", lty = 2)
legend("topright", title = expression(beta), legend = c("1", ".90"), lty = 1:2, col = "black")
dev.off()


chain1$out <- chain1$out[1:1e4, ]
chain9$out <- chain9$out[1:1e4, ]

# chain1$loops <- chain1$loops[1:1e4, ]
# chain9$loops <- chain9$loops[1:1e4, ]

pdf("corr_trace1.pdf", height = 5, width = 5)
plot.ts(chain1$out[,7], col = adjustcolor("black", alpha.f = .70), ylab = expression( mu), ylim = range(-.2, 2.5))
lines(chain9$out[,7], col = adjustcolor("red", alpha.f = .50), lty = 1)
legend("topright", title = expression(beta), legend = c("1", ".90"), lty = 1, col = c(adjustcolor("black", alpha.f = .70), adjustcolor("red", alpha.f = .50)) )
dev.off()

pdf("corr_trace2.pdf", height = 5, width = 5)
plot.ts(chain1$out[,8], col = adjustcolor("black", alpha.f = .70), ylab = expression( sigma^2) )
lines(chain9$out[,8], col = adjustcolor("red", alpha.f = .50), lty = 1)
legend("topright", title = expression(beta), legend = c("1", ".90"), lty = 1, col = c(adjustcolor("black", alpha.f = .70), adjustcolor("red", alpha.f = .50)) )
dev.off()

ind <- seq(1, 1e5, 1)
pdf("corr_loops1.pdf", height = 5, width = 5)
plot(ind, log(chain1$loops[ind,1]), col = adjustcolor("black", alpha.f = .30), ylab = expression(paste("Log of the number of loops for ", mu)), pch = 16)
points(ind, log(chain9$loops[ind,1]), col = adjustcolor("red", alpha.f = .30), lty = 1, pch = 18)
legend("topleft", horiz = TRUE, title = expression(beta), legend = c("1", ".90"), pch = c(20,18), col = c(adjustcolor("black", alpha.f = .70), adjustcolor("red", alpha.f = .70)) )
dev.off()

pdf("corr_loops2.pdf", height = 5, width = 5)
plot(ind,log(chain1$loops[ind,2]), col = adjustcolor("black", alpha.f = .30), ylab = expression(paste("Log of the number of loops for ", sigma^2)), pch = 16, ylim = c(0, 5))
points(ind, log(chain9$loops[ind,2]), col = adjustcolor("red", alpha.f = .30), lty = 1, pch = 18)
legend("topright", horiz = TRUE, title = expression(beta), legend = c("1", ".90"), pch = c(20,18), col = c(adjustcolor("black", alpha.f = .70), adjustcolor("red", alpha.f = .70)) )
dev.off()



## For Shadow only
load("compare_shadow")
pdf("mu_sig_compare.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(density(chain1$out[,7]), col = "black", xlim = range(c(chain$out[, 7], chain9$out[, 7], chain1$out[, 7])), lty = 1, xlab = expression(mu), main = "", ylim = c(0, 1.7))
lines(density(chain9$out[,7]), col = "black", lty = 2)
lines(density(chain$out[,7]), col = "black", lty = 3, lwd = 2.5)
legend("topright",legend = c("Barker's", "Portkey", "Shadow"), lty = 1:3, col = "black", lwd = c(1,1,2.5))

plot(density(chain1$out[,8]), col = "black", xlim = c(0,2), lty = 1, xlab = expression(sigma^2), main = "", ylim = c(0,2.6))
lines(density(chain9$out[,8]), col = "black", lty = 2)
lines(density(chain$out[,8]), col = "black", lty = 3, lwd = 2.5)
legend("topright",legend = c("Barker's", "Portkey", "Shadow"), lty = 1:3, col = "black", lwd = c(1,1,2.5))
dev.off()	



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



