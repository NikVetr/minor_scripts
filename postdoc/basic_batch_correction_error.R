nrep <- 1E5
fitlm <- function(n = 20){
  b1 <- 0 #true effect of x
  b2 <- 0 #true effect of batch
  a <- 0
  e <- rnorm(n)
  batch <- sample(c(rep(0, n/2), rep(1, n/2)))
  x <- numeric(n) 
  x[batch==0] <- sample(c(rep(0, n/4), rep(1, n/4))) #allocate treatments to batches
  x[batch==1] <- sample(c(rep(0, n/4), rep(1, n/4)))
  y <- a + b1 * x + b2 * batch + e
  fit1 <- lm(y ~ 1 + x + batch)
  y_c <- y - fit1$coefficients["batch"] * batch
  fit2 <- lm(y_c ~ 1 + x)
  return(c(p1 = summary(fit1)$coef["x", 4],
           p2 = summary(fit2)$coef["x", 4]))
}


pvals <- t(replicate(nrep, expr = fitlm()))
par(mfrow = c(2,1))
hist(pvals[,1], xlab = "p-value for coefficient on x", main = "batch included in model")
hist(pvals[,2], xlab = "p-value for coefficient on x", main = "batch regressed out from outcome")
