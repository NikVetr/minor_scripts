#lkj density
rlkj <- function (d, eta = 1) {
  alpha <- eta + (d - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, d, d)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (d > 2) 
    for (m in 2:(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

lkj_normalizing_constant <- function(d, eta){
  2^(sum(sapply(1:(d-1), function(k) (2 * (eta - 1) + d - k) * (d - k)))) *
    prod(sapply(1:(d-1), function(k) beta(a = eta + (d - k - 1) / 2, b = eta + (d - k - 1) / 2)^(d-k)))
}

dlkj <- function(C, eta){
  d <- dim(C)[1]
  normalizing_constant <- lkj_normalizing_constant(d, eta)
  propto_dens <- det(C)^(eta-1)
  return(normalizing_constant * propto_dens)
}


#MVN monte carlo integral
p <- 2
r <- 0.5
R <- diag(p)*(1-r) + r
nrunif <- 1E6
bound <- 10
x <- matrix(runif(nrunif*p, min = -bound, max = bound), nrunif, p)
mean(mvtnorm::dmvnorm(x, sigma = R, mean = rep(0, p))) * (2 * bound)^p

#lkj monte carlo integral
p <- 3
eta <- 3
nrunif <- 1E4
runif_corrmats <- replicate(nrunif, rlkj(p, eta = 1))
densities <- sapply(1:nrunif, function(i) dlkj(runif_corrmats[,,i], eta = eta))
mean(densities)# * lkj_normalizing_constant(d = p, eta = 1)
