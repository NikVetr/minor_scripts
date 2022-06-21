r <- 2.1
time = 1000
nrep <- 1E5
hist(replicate(sum(rpois(1000, lambda = r)), n = nrep), breaks = 100)
hist(rpois(nrep, lambda = r*time), breaks = 100, add = T)

nrep2 <- 1E3
sojourn <- unlist(replicate(diff(sort(runif(rpois(1, lambda = r*time), min = 0, max = time))), n = nrep2))
hist(sojourn, breaks = 5E2, freq = F)
h_range <- 0:(max(sojourn)*1000)/1000
lines(h_range, dexp(h_range, rate = r), col = 2, lwd = 3)
probs <- 0:100/100
plot(quantile(sojourn, probs = probs), qexp(probs, rate = r), type = "l", lwd = 2, col = 2, xlim = c(0,2), ylim = c(0,2), 
     main = "quantile-quantile plot", xlab = "empirically simulated quantiles", ylab = "analytically derived quantiles") 
abline(0,1,lty = 2)
