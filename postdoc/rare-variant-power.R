library(pwr)
library(parallel)
library(IHW)

maf2pwr <- function(maf, cohens_d = 2, desired_power = 0.8) {
  p <- 2
  x_var <- maf * (1-maf) * 2
  cov_xy <- cohens_d  * x_var
  y_var <- 1
  r2 <- cov_xy^2 / x_var / y_var
  pwr <- pwr.f2.test(u = 1, v = NULL, f2 = r2 / (1-r2), sig.level = 0.05, power = desired_power)
  n <- pwr$v + p
  n
}

mafs <- 1 / (3:500)^2
plot(-log10(mafs), sapply(mafs, maf2pwr), type = "l")


mafs <- 1 / 2^(1:15)
nrep <- 2E3
nsamp <- 2E4
b <- 0.1

sim <- function(maf, b, nsamp){
  x <- rbinom(nsamp, 2, maf)
  if(all(x == 0)){return(c(maf,0,1))} #invariant site
  y <- x * b + rnorm(nsamp)
  output <- summary(lm(y ~ x))$coefficients[2,c(1,4)]
  return(c(maf, output))
}

out <- mclapply(mafs, function(maf) t(replicate(nrep,sim(maf, b, nsamp))), mc.cores = length(mafs))
out <- as.data.frame(do.call(rbind, out))
colnames(out) <- c("maf", "bhat", "pval")
pvaladj <- ihw(out$pval, as.factor(out$maf), alpha = 0.05)
out$adjpval <- pvaladj@df$adj_pvalue
sigout <- out[out$adjpval < 0.05,]
plot(-log10(sigout$maf), sigout$bhat, col = adjustcolor(1, 0.25), pch = c(4, 19)[((sigout$bhat * sign(b)) > 0) + 1],
     main = paste0("true eff size = ", b, ", sample size = ", nsamp))
abline(h = b, col = "green", lty = 2)
mean_abs_bhats <- tapply(abs(sigout$bhat), sigout$maf, mean)
lines(-log10(mafs), rev(as.numeric(mean_abs_bhats)), col = "purple", lwd = 2)
legend(x = "bottomleft", legend = c("mean significant absolute estimate", "true parameter value"), 
       col = c("purple", "green"), lty = c(1,2), lwd = c(2,1))

