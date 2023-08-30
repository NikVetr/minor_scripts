
lmpois <- function(n = 1E3, p = 3, b = 0){
  e <- rnorm(n)
  x <- replicate(p, rbinom(n, 1, 0.5))
  bs <- matrix(rep(b,p), p, 1)
  l <- exp(1 + x %*% bs  + e)
  y <- rpois(n, l)
  summary(lm(y ~ 1 + x))$coefficients[-1,4]
}

glmpois <- function(n = 1E3, p = 3, b = 0){
  e <- rnorm(n)
  x <- replicate(p, rbinom(n, 1, 0.5))
  bs <- matrix(rep(b,p), p, 1)
  l <- exp(1 + x %*% bs  + e)
  y <- rpois(n, l)
  summary(glm(y ~ 1 + x, family = poisson))$coefficients[-1,4]
}


hist(unlist(replicate(1E3, lmpois(n = 50, p = 10))))
hist(unlist(replicate(1E3, glmpois(n = 50, p = 10))))
