n <- 100
p <- 50
nrep <- 150

a <- 0
b <- matrix(0, p, 1)
b <- rnorm(p)
pvals <- t(replicate(nrep, {
  x <- matrix(rnorm(n*p), nrow = n, ncol = p)
  e <- rnorm(n, sd = sqrt(n))
  y <- a + x %*% b + e
  fit <- lm(y ~ 1 + x)
  summ <- summary(fit)
  summ$coefficients[,"Pr(>|t|)"]
}))

mean(pvals < 0.05)
mean(pvals < 0.05 / n)
mean(pvals < 0.05 / p / n)

mean(p.adjust(pvals, method = "BH") < 0.05)
mean(p.adjust(pvals, method = "BY") < 0.05)
mean(p.adjust(pvals, method = "fdr") < 0.05)

ihw_p <- IHW::ihw(c(pvals), as.factor(rep(1:nrep, p+1)), alpha = 0.05)
ihw_p_inv <- IHW::ihw(c(t(pvals)), as.factor(rep(1:(p+1), nrep)), alpha = 0.05)
mean(ihw_p@df$pvalue < 0.05)
mean(ihw_p_inv@df$pvalue < 0.05)
