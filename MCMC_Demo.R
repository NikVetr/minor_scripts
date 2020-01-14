#quick demo of different MCMC algorithms

#gibbs
condMVN <- function(means, cov, obs, inds){
  A <- cov[-inds, -inds]
  C <- cov[inds, inds]
  B <- cov[-inds, inds]
  condCov <- A - B%*%solve(C)%*%t(B)
  condMean <- as.vector(means[-inds] + B%*%solve(C)%*%(obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}

cov <- matrix(c(2,0.7,0.7,1), 2, 2)
mean <- c(3, -1.5)
n_samps <- 5000
samps <- matrix(rep(NA, 2*n_samps), n_samps, 2)
samps[1,] <- c(0,0)
for(i in 2:n_samps){
  if(i %% (n_samps / 100) == 0){print(paste0(i / n_samps * 100, "%"))}
  curr_param <- ifelse(i %% 2 == 0, 1, 2)
  prev_param <- setdiff(c(1,2), curr_param)
  prev_samp <- samps[i-1,]
  cond <- condMVN(mean, cov, prev_samp[prev_param], prev_param)
  curr_samp <- c(NA, NA)
  curr_samp[curr_param] <- rnorm(1, cond$means, sqrt(cond$cov))
  curr_samp[prev_param] <- samps[i-1, prev_param]
  samps[i,] <- curr_samp
}

apply(samps, 2, mean)
cov(samps)
par(mfrow = c(2,1))
plot(samps[,1], type = "l"); abline(h = mean[1], col = 2); plot(samps[,2], type = "l"); abline(h = mean[2], col = 2)

#met-hastings
n_samps <- 1E5
samps <- matrix(rep(NA, 2*n_samps), n_samps, 2)
window_size <- 1
samps[1,] <- c(0,0)
library(mvtnorm)

curr_log_dens <- dmvnorm(x = samps[1,], mean = mean, sigma = cov, log = T)
for(i in 2:n_samps){
  if(i %% (n_samps / 100) == 0){print(paste0(i / n_samps * 100, "%"))}
  param <- sample(x = 1:dim(cov)[1], size = 1)
  curr_samp <- samps[i-1,]
  prop_samp <- curr_samp
  prop_samp[param] <- prop_samp[param] + runif(min = -window_size/2, max = window_size/2, n = 1)
  log_prop_ratio <- 0
  prop_log_dens <- dmvnorm(x = prop_samp, mean = mean, sigma = cov, log = T)
  log_dens_ratio <- prop_log_dens - curr_log_dens 
  accept_prob <- exp(log_dens_ratio + log_prop_ratio)
  accept_roll <- runif(1)
  if(accept_prob > accept_roll){
    curr_samp <- prop_samp
    curr_log_dens <- prop_log_dens
  }
  samps[i,] <- curr_samp
}

apply(samps, 2, mean)
cov(samps)
par(mfrow = c(2,1))
plot(samps[,1], type = "l"); abline(h = mean[1], col = 2); plot(samps[,2], type = "l"); abline(h = mean[2], col = 2)
