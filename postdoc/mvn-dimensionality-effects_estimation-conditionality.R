conditional_mvn <- function(mu, sigma, indices, observed_values) {
  # Check input dimensions
  p <- length(mu)
  if (!all(dim(sigma) == c(p, p))) {
    stop("sigma matrix must be p x p.")
  }
  if (length(indices) != length(observed_values)) {
    stop("Length of indices must match length of observed values.")
  }
  
  # Identify indices for unobserved and observed variables
  all_indices <- 1:p
  unobserved_indices <- setdiff(all_indices, indices)
  
  # Partition mu vector and sigma matrix
  mu1 <- mu[unobserved_indices]   # mu of unobserved variables
  mu2 <- mu[indices]              # mu of observed variables
  
  Sigma11 <- sigma[unobserved_indices, unobserved_indices] # sigma of unobserved
  Sigma12 <- sigma[unobserved_indices, indices]            # Cross-sigma
  Sigma21 <- sigma[indices, unobserved_indices]
  Sigma22 <- sigma[indices, indices]                       # sigma of observed
  Sigma22_inv <- solve(Sigma22) # Inverse of Sigma22
  
  # Compute conditional mu and sigma
  conditional_mu <- mu1 + Sigma12 %*% Sigma22_inv %*% (observed_values - mu2)
  conditional_sigma <- Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
  
  # Return result as a list
  list(mu = as.vector(conditional_mu), sigma = conditional_sigma)
}


n <- 50
p <- 40
mu <- rnorm(p)
r <- 0.75
R <- diag(p) + r - diag(p) * r
sds <- diag(rexp(p, 0.1)^1)
sigma <-  sds %*% R %*% sds
L <- t(chol(sigma))
x <- t(L %*% matrix(rnorm(p*n), p, n) + mu)
est_mu <- apply(x, 2, mean)
est_sigma <- cov(x)
est_sigma <- corpcor::cov.shrink(cov(x))

#induce outliers 
prop_outliers <- 0.1
n_assess <- 20
for(i in 1:n_assess){
  x[i,] <- x[i,] - (x[i,] - mu) * sample(c(rep(1, round(p * prop_outliers)), rep(0, p - round(p * prop_outliers)))) * 2
}

zscores <- data.frame(do.call(rbind, parallel::mclapply(1:n_assess, function(j){ 
  do.call(rbind, lapply(1:p, function(i){
    true_cond <- conditional_mvn(mu = mu, sigma = sigma, indices = (1:p)[-i], observed_values = x[j,-i])
    est_cond <- conditional_mvn(mu = est_mu, sigma = est_sigma, indices = (1:p)[-i], observed_values = x[j,-i])
    true_cz <- (x[j,i] - true_cond$mu) / sqrt(true_cond$sigma)
    est_cz <- (x[j,i] - est_cond$mu) / sqrt(est_cond$sigma)
    true_mz <- (x[j,i] - mu[i]) / sqrt(sigma[i,i])
    est_mz <- (x[j,i] - est_mu[i]) / sqrt(est_sigma[i,i])
    return(c(true_cz = true_cz, est_cz = est_cz,
             true_mz = true_mz, est_mz = est_mz))
  }))
}, mc.cores = 12)))

par(mfrow = c(2,4))
plot(zscores$true_cz, zscores$est_cz)
abline(0,1)
plot(zscores$true_mz, zscores$est_mz)
abline(0,1)

plot(zscores$true_mz, zscores$true_cz)
abline(0,1)
plot(zscores$est_mz, zscores$est_cz)
abline(0,1)

for(i in 1:4){
  hist(zscores[,i], main = colnames(zscores)[i], breaks = 100)
}
