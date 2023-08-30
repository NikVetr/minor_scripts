# There is a single position to fill.
# 
# There are n applicants for the position, and the value of n is known.
# 
# The applicants, if all seen together, can be ranked from best to worst 
# unambiguously.
# 
# The applicants are interviewed sequentially in random order, with each order 
# being equally likely.
# 
# Immediately after an interview, the interviewed applicant is either accepted or 
# rejected, and the decision is irrevocable.
# 
# The decision to accept or reject an applicant can be based only on the relative 
# ranks of the applicants interviewed so far.
# 
# The objective of the general solution is to have the highest probability of selecting 
# the best applicant of the whole group. This is the same as maximizing the expected payoff, 
# with payoff defined to be one for the best applicant and zero otherwise.

sim_apps_canon <- function(){
  mu <- rnorm(1, sd = 100)
  sigma <- invgamma::rinvgamma(1, 1, 1)
  n <- 1E2
  apps <- rnorm(n, mu, sigma)
  cutoff <- floor(n / exp(1))
  thresh <- max(apps[1:cutoff])
  winner <- which(apps[(cutoff+1):n] > thresh)
  
  if(length(winner) > 0){
    winner_id = min(winner) + cutoff
    max_id <- which.max(apps)
    return(as.numeric(winner_id == max_id))
  } else {
    return(0)
  }
}

mean(replicate(1E4, sim_apps_canon()))


#alternatively, why not learn the distribution at finer grain
#and accept a candidate once we are >50% sure they are the top candidate in the sample
#by learning eg the mean and sd of the sample
#finding the posterior distribution of each new candidate's value
#and evaluating whether a more impressive candidate is likely to exist in the
#rest of the sample using rank statistics?

#can also extend it to having interviewed past secretaries before, 
#and retaining diffuse prior from those experiences

#first let's try the omniscient case, where we know the true parameter values
sim_apps_omnisc <- function(){
  n <- 1E2
  apps <- runif(n, 0, 1)
  thresh <- 0.5
  prob_best_remainder <- sapply(1:n, function(i) 1 - apps[i]^(n-i))
  winner <- which(prob_best_remainder < thresh)
  
  if(length(winner) > 0){
    winner_id = min(winner)
    max_id <- which.max(apps)
    return(as.numeric(winner_id == max_id))
  } else {
    return(0)
  }
}

mean(unlist(parallel::mclapply(1:1E4, function(i) sim_apps_omnisc(), mc.cores = 14)))

#now let's try the case where we do not, using sample statistics & the known dist
sim_apps_sample <- function(){
  mu <- rnorm(1, sd = 100)
  sigma <- invgamma::rinvgamma(1, 1, 1)
  n <- 1E3
  apps <- rnorm(n, mu, sigma)
  mu_samples <- sapply(1:n, function(i) mean(apps[1:i]))
  sd_samples <- sapply(1:n, function(i) sd(apps[1:i]))
  app_qs <- sapply(1:n, function(i) pnorm(apps[i], mu_samples[i], sd_samples[i]))
  app_qs[1] <- 0
  thresh <- 0.5
  prob_best_remainder <- sapply(1:n, function(i) 1 - app_qs[i]^(n-i))
  winner <- which(prob_best_remainder < thresh)
  
  if(length(winner) > 0){
    winner_id = min(winner)
    max_id <- which.max(apps)
    return(as.numeric(winner_id == max_id))
  } else {
    return(0)
  }
}

mean(unlist(parallel::mclapply(1:1E3, function(i) sim_apps_sample(), mc.cores = 14)))


#now let's try the case where we do not, using a KDE
sim_apps_kde <- function(){
  mu <- rnorm(1, sd = 100)
  sigma <- invgamma::rinvgamma(1, 1, 1)
  n <- 1E4
  apps <- rnorm(n, mu, sigma)
  r_apps <- range(apps) + c(-1,1) * diff(range(apps)) / 10
  kdes <- c(NA, lapply(2:n, function(i) 
    density(apps[1:i], from = r_apps[1], to = r_apps[2], n = 2048)
  ))
  dx <- diff(kdes[[i]]$x[1:2])
  app_qs <- c(0, sapply(2:n, function(i){
    hits <- kdes[[i]]$x < apps[i]
    sum(dx * (kdes[[i]]$y[hits][-sum(hits)] + 
                diff(kdes[[i]]$y[hits]) / 2))
  }))
  thresh <- 0.5
  prob_best_remainder <- sapply(1:n, function(i) 1 - app_qs[i]^(n-i))
  winner <- which(prob_best_remainder < thresh)
  
  if(length(winner) > 0){
    winner_id = min(winner)
    max_id <- which.max(apps)
    return(as.numeric(winner_id == max_id))
  } else {
    return(0)
  }
}

mean(unlist(parallel::mclapply(1:1E3, function(i) sim_apps_kde(), mc.cores = 14)))

#now let's using a Bayesian update from the true prior
sim_apps_bayes <- function(){
  mu <- rnorm(1, sd = 100)
  sigma <- invgamma::rinvgamma(1, 1, 1)
  n <- 1E2
  apps <- rnorm(n, mu, sigma)
  cutoff <- floor(n / exp(1))
  thresh <- max(apps[1:cutoff])
  winner <- which(apps[(cutoff+1):n] > thresh)
  
  if(length(winner) > 0){
    winner_id = min(winner) + cutoff
    max_id <- which.max(apps)
    return(as.numeric(winner_id == max_id))
  } else {
    return(0)
  }
}

#now let's using a Bayesian update from a more diffuse prior


#now let's using a Bayesian update from a more informative, historical, "false" prior

#now let's try a "distribution free" approach, leveraging chebyshev's inequality
