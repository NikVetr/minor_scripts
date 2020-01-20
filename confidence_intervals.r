
#sampling distribution of mean of normally distributed random variable -- symmetry argument
mu <- 7
sd <- 2

nrep <- 1E8
nsamp <- 1E2
x <- replicate(nrep, rnorm(nsamp, mu, sd))
t_97.5 <- qt(0.975, df = nsamp - 1)
t_2.5 <- qt(0.025, df = nsamp - 1)
upper <- sapply(1:nrep, function(rep) mean(x[,rep]) + t_97.5 * sd(x[,rep]) / sqrt(length(x[,rep])))
lower <- sapply(1:nrep, function(rep) mean(x[,rep]) + t_2.5 * sd(x[,rep]) / sqrt(length(x[,rep])))
sum(mu < upper & mu > lower) / nrep

#for a linear model y = a + bx
library(corpcor)
a <- 3
b <- 2

nrep <- 1E8
nsamp <- 10
x <- replicate(nrep, rnorm(nsamp, 5, 1))
y <- a + b * x + rnorm(nsamp)
A <- lapply(1:nrep, function(rep) cbind(rep(1, nsamp), x[,rep]))
coefs <- sapply(1:nrep, function(rep) solve(t(A[[rep]]) %*% A[[rep]]) %*% t(A[[rep]]) %*% y[,rep])
resids <- sapply(1:nrep, function(rep) y[,rep] - coefs[1,rep] + coefs[2,rep] * x[,rep])
Sx <- sapply(1:nrep, function(rep) sum(x[,rep]))
Sy <- sapply(1:nrep, function(rep) sum(y[,rep]))
Sxx <- sapply(1:nrep, function(rep) sum(x[,rep]^2))
Sxy <- sapply(1:nrep, function(rep) sum(x[,rep]*y[,rep]))
Syy <- sapply(1:nrep, function(rep) sum(y[,rep]^2))
Sr <- sapply(1:nrep, function(rep) 
  sqrt(1 / (nsamp * (nsamp - 2)) * (nsamp * Syy[rep] - Sy[rep]^2 - coefs[2,rep]^2 * (nsamp * Sxx[rep] - Sx[rep]^2))))
Sb <- sapply(1:nrep, function(rep) sqrt(nsamp * Sr[rep]^2 / (nsamp * Sxx[rep] - Sx[rep]^2)))
Sa <- sapply(1:nrep, function(rep) sqrt(Sb[rep]^2 / nsamp * Sxx[rep]))

t_97.5 <- qt(0.975, df = nsamp - 2)
t_2.5 <- qt(0.025, df = nsamp - 2)

upper <- sapply(1:nrep, function(rep) coefs[2,rep] + t_97.5 * Sb[rep])
lower <- sapply(1:nrep, function(rep) coefs[2,rep] + t_2.5 * Sb[rep])
sum(b < upper & b > lower) / nrep

upper <- sapply(1:nrep, function(rep) coefs[1,rep] + t_97.5 * Sa[rep])
lower <- sapply(1:nrep, function(rep) coefs[1,rep] + t_2.5 * Sa[rep])
sum(a < upper & a > lower) / nrep
