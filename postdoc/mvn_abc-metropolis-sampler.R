# Set up the parameters
n_iter <- 1E5
p_burnin <- 0.2
n_burnin <- n_iter * p_burnin
n <- 200  
a_obs <- rnorm(n)  
sigma <- matrix(c(1, 1, 1, 2), ncol=2)
thin <- 10 

# Precompute the inverse upper cholesky factor of the covariance matrix
sigma_chol_inv <- solve(chol(sigma))

# Initialize b as a matrix where each row is a sampled vector at each iteration
b <- matrix(NA, nrow=n_iter, ncol=n)
b[1, ] <- rnorm(n)  # Initial values for b
c_obs <- a_obs + b[1, ]  # Initial values of c

# Function to calculate the log of the target density for a vector of b's
log_target_density <- function(b_vec, a_obs, sigma_chol_inv) {
  c_vec <- a_obs + b_vec
  total_log_density <- sum(dnorm(cbind(a_obs, c_vec) %*% sigma_chol_inv, log = T))
  return(total_log_density)
}

# Metropolis sampling
log_target_densities <- numeric(n_iter)
log_target_densities[1] <- log_target_density(b[1,], a_obs, sigma_chol_inv)
accept_reject <- logical(n_iter)
for (i in 2:n_iter) {
  #progress
  if(i %% round(n_iter/20) == 0) print(i)
  
  # Symmetric proposal distribution for each element of b
  n_proposal <- sample(1:min(n, 10), 1)
  inds_proposal <- sample(1:n, n_proposal)
  b_proposal <- b[i-1, ]
  b_proposal[inds_proposal] <- b_proposal[inds_proposal] + rnorm(n_proposal, mean=0, sd = 5/n_proposal)  # Tune the sd as needed
  b_proposal_log_target_density <- log_target_density(b_proposal, a_obs, sigma_chol_inv)
  
  # Calculate acceptance ratio
  log_accept_ratio <- b_proposal_log_target_density - log_target_densities[i-1]
  
  # Accept or reject the proposal
  if (log(runif(1)) < log_accept_ratio) {
    b[i, ] <- b_proposal
    log_target_densities[i] <- b_proposal_log_target_density
    accept_reject[i] <- T
  } else {
    b[i, ] <- b[i-1, ]
    log_target_densities[i] <- log_target_densities[i-1]
  }
}
mean(accept_reject)

#remove burnin
mean(accept_reject)
b_samps <- b[seq((n_burnin+1), n_iter, by = 10),]
plot(b_samps[,2], type = "l")

par(mfrow = c(2,1))
hist(rowMeans(b_samps), xlab = "posterior means of individual b_i", freq = F, main = ""); abline(v=1, col=2)
hist(apply(b_samps, 1, var), xlab = "posterior variances of all b_i", freq = F, main = "")
b_means <- colMeans(b_samps)
c_means <- colMeans(t(t(b_samps) + a))

pairs(cbind(a=a, b=b_means, c=c_means), 
      diag.panel = panel.hist, 
      lower.panel = panel.pts, upper.panel = NULL, cex.labels = 3)
