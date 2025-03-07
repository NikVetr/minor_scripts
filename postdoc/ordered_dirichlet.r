truncated_beta_inverse_sampling <- function(size, shape1, shape2, lb, ub) {
  Fa <- pbeta(lb, shape1, shape2)  # CDF at truncation point
  Fb <- pbeta(ub, shape1, shape2)  # CDF at truncation point
  U <- runif(size, Fa, Fb)            # Uniform samples in the truncated range
  X <- qbeta(U, shape1, shape2)   # Inverse CDF
  return(X)
}

ordered_dirichlet <- function(K, alpha = NA) {
  x <- numeric(K)
  remaining_mass <- 1.0
  if(all(is.na(alpha))){
    alpha <- rep(1, K)
  }
  if(length(alpha) != K){
    alpha <- rep(alpha, K)
  }
  mass_excluded <- numeric(K-1)
  for (k in 1:(K - 1)) {
    lb <- remaining_mass / (K - k + 1)
    ub <- min(ifelse(k > 1, x[k - 1], 1.0), remaining_mass)
    beta_bounds <- c(lb, ub) / remaining_mass
    
    # transform uniform sample to truncated Beta by sampling a quantile
    shapes <- c(alpha[k], sum(alpha[(k+1):K]))
    beta_quantiles <- pbeta(beta_bounds, shapes[1], shapes[2])
    
    #sample quantile uniformly
    u <- runif(1)
    
    #maybe if we weigh the draw in [0,1] by the flipped shape?
    u <- rbeta(1, shapes[2], shapes[1])
    
    #what if we stretch this out from just the region of truncation? 
    #equiv to weighing by the density in the truncated region
    # u <- truncated_beta_inverse_sampling(size = 1, shape1 = shapes[2], shape2 = shapes[1], lb = lb, ub = ub)
    # u <- (u - lb) / (ub - lb) #now u is sampled in [0,1]
    
    #no this also does not work. What if we increment target (in Stan) by how much the remaining space shrinks?
    #ie how much of the mass of the remaining Betas gets removed by each successive element? 
    
    #what if we use order statistics for the r'th ranked shape in the remainder?
    # u <- rbeta(1, K-k+1, k)
    # this is the marginal for the element, not the broken stick
    # this one is more broken-stick-like (rank of the remaining dirichlet of the rest of the broken stick)
    u <- rbeta(1, K-k+1, 1)
    
    #what if we take the quantile from *just* within the valid region?
    
    beta_q <- beta_quantiles[1] + u * diff(beta_quantiles)
    truncated_beta_sample <- qbeta(beta_q, shapes[1], shapes[2])
    
    # update
    mass_excluded[k] <- diff(beta_quantiles)
    mass_excluded[k] <- remaining_mass
    x[k] <- truncated_beta_sample * remaining_mass
    remaining_mass <- remaining_mass - x[k]
  }
  x[K] <- remaining_mass
  return(list(x = x, mass_excluded = mass_excluded))
}

ordered_dirichlet_2 <- function(K, alpha = NA) {
  x <- numeric(K)
  remaining_mass <- 1.0
  if(all(is.na(alpha))){
    alpha <- rep(1, K)
  }
  if(length(alpha) != K){
    alpha <- rep(alpha, K)
  }
  mass_excluded <- numeric(K-1)
  for (k in 1:(K - 1)) {
    lb <- remaining_mass / (K - k + 1)
    ub <- min(ifelse(k > 1, x[k - 1], 1.0), remaining_mass)
    beta_bounds <- c(lb, ub) / remaining_mass
    
    # marginal distribution of that element
    shapes <- c(alpha[k], sum(alpha[-k]))
    beta_quantiles <- pbeta(beta_bounds, shapes[1], shapes[2])
    
    #quantile of that element
    u <- runif(1, 0, 1)
    beta_q_shapes <- c(K-k+1, k)
    beta_q_quantiles <- pbeta(beta_quantiles, beta_q_shapes[1], beta_q_shapes[2])
    beta_q_q <- beta_q_quantiles[1] + u * diff(beta_q_quantiles)
    
    #value of that element
    beta_q <- beta_quantiles[1] + beta_q_q * diff(beta_quantiles)
    truncated_beta_sample <- qbeta(beta_q, shapes[1], shapes[2])
    
    # update
    mass_excluded[k] <- diff(beta_quantiles)
    x[k] <- truncated_beta_sample * remaining_mass
    remaining_mass <- remaining_mass - x[k]
  }
  x[K] <- remaining_mass
  return(list(x = x, mass_excluded = mass_excluded))
}


ordered_dirichlet_emp <- function(K, alpha = NA, bpar = NA) {
  x <- numeric(K)
  remaining_mass <- 1.0
  if(all(is.na(alpha))){
    alpha <- rep(1, K)
  }
  if(length(alpha) != K){
    alpha <- rep(alpha, K)
  }
  mass_excluded <- numeric(K-1)
  for (k in 1:19) {
    lb <- remaining_mass / (K - k + 1)
    ub <- min(ifelse(k > 1, x[k - 1], 1.0), remaining_mass)
    beta_bounds <- c(lb, ub) / remaining_mass
    
    # transform uniform sample to truncated Beta by sampling a quantile
    shapes <- c(alpha[k], sum(alpha[(k+1):K]))
    beta_quantiles <- pbeta(beta_bounds, shapes[1], shapes[2])
    
    #sample quantile uniformly
    beta_q_rm <- rbeta(1, shape1 = bpar$shape1[k], shape2 = bpar$shape2[k])
    lb_rm <- lb / remaining_mass
    ub_rm <- ub / remaining_mass
    betaq_lb <- pbeta(lb_rm, shapes[1], shapes[2])
    betaq_ub <- pbeta(ub_rm, shapes[1], shapes[2])
    beta_q <- beta_q_rm * (betaq_ub - betaq_lb) + betaq_lb
    
    # update
    x[k] <- qbeta(beta_q, shapes[1], shapes[2]) * remaining_mass
    remaining_mass <- remaining_mass - x[k]
  }
  x[K] <- remaining_mass
  return(list(x = x, mass_excluded = mass_excluded))
}

unordered_dirichlet <- function(K, alpha = NA) {
  x <- numeric(K)
  remaining_mass <- 1.0
  if(all(is.na(alpha))){
    alpha <- rep(1, K)
  }
  if(length(alpha) != K){
    alpha <- rep(alpha, K)
  }
  for (k in 1:(K - 1)) {
    beta_bounds <- c(0,1)
    
    # transform uniform sample to truncated Beta
    shapes <- c(alpha[k], sum(alpha[(k+1):K]))
    beta_quantiles <- pbeta(beta_bounds, shapes[1], shapes[2])
    u <- runif(1)
    beta_q <- beta_quantiles[1] + u * diff(beta_quantiles)
    truncated_beta_sample <- qbeta(beta_q, shapes[1], shapes[2])

    # update
    x[k] <- truncated_beta_sample * remaining_mass
    remaining_mass <- remaining_mass - x[k]
  }
  
  x[K] <- remaining_mass
  # x2[K] <- remaining_mass
  return(x)
}


rejection_sample_dirichlet <- function(K, alpha, n_samples = 1000) {

  if(all(is.na(alpha))){
    alpha <- rep(1, K)
  }
  if(length(alpha) != K){
    alpha <- rep(alpha, K)
  }
  
  accepted_samples <- list()
  while (length(accepted_samples) < n_samples) {
    sample <- as.numeric(extraDistr::rdirichlet(1, alpha))
    if (all(diff(sample) <= 0)) {
      accepted_samples <- append(accepted_samples, list(sample))
    }
    # if (sample[1] >= 1/K) {
    #   accepted_samples <- append(accepted_samples, list(sample))
    # }
  }
  
  return(do.call(rbind, accepted_samples))
}

# Parameters
K <- 20
alpha <- rep(1, K)
n_samples <- 2E3

# Generate samples
ordered_samples <- lapply(1:n_samples, function(i) ordered_dirichlet(K, alpha))
ordered_samples <- lapply(1:n_samples, function(i) ordered_dirichlet_emp(K, alpha, bfit_par)) #bfit_par from below
ordered_ME <- do.call(rbind, lapply(ordered_samples, function(x) x$mass_excluded))
ordered_samples <- do.call(rbind, lapply(ordered_samples, function(x) x$x))
# unordered_samples <- t(replicate(n_samples, unordered_dirichlet(K, alpha)))
# filtered_samples <- unordered_samples[apply(unordered_samples, 1, function(y) all(diff(y) <= 0)),]
# rejected_samples <- rejection_sample_dirichlet(K, alpha, n_samples)
sorted_samples <- t(apply(extraDistr::rdirichlet(n_samples, alpha), 1, sort, decreasing = T))

#reweigh samples by some property?
ordered_lME <- rowSums(log(1-ordered_ME[,-1]))
ordered_lME <- rowSums(log(ordered_ME))
lME_weight <- exp(ordered_lME - mean(ordered_lME))
# ordered_samples <- ordered_samples[sample(1:n_samples, n_samples, replace = T, lME_weight),]

#plot?
i <- 2
# breaks_range <- range(c(ordered_samples[,i], filtered_samples[,i], rejected_samples[,i], sorted_samples[,i]))
ifelse2 <- function(test, yes, no) if(test) yes else no
breaks_range <- range(c(ordered_samples[,i], sorted_samples[,i]))
breaks <- seq(from = breaks_range[1], to = breaks_range[2], length.out = 20)
hist(ordered_samples[,i], freq = F, breaks = breaks, 
     xlim = ifelse2(diff(breaks_range) > 0.25, c(0,1), breaks_range), 
     main = paste0("marginal distributions for element ", i))
hist(sorted_samples[,i], add = T, col = adjustcolor(3,0.5), breaks = breaks, freq = F)

#ok, so rejection sampling of later entries influences acceptance of earlier ones
#if we rejection sample only on the first element, we do get the correct distribution
#but if the first element is too small (close to the bound) the later elements are likely
#to be too large, thus leading to rejection. Maybe we can sample the uniform RV from a
#non-uniform Beta, or add some sort of jacobian adjustment to fix this? eg increment target
#by the log (or negative log) Beta density or something?


#plot the marginal means and variances
ordered_means <- colMeans(ordered_samples)
ordered_vars <- apply(ordered_samples, 2, var)
sorted_means <- colMeans(sorted_samples)
sorted_vars <- apply(sorted_samples, 2, var)

# Rescale function to map values to [0, 1]
rescale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Rescale both datasets to [0, 1] range
ordered_means_rescaled <- rescale(ordered_means)
sorted_means_rescaled <- rescale(sorted_means)
ordered_vars_rescaled <- rescale(ordered_vars)
sorted_vars_rescaled <- rescale(sorted_vars)

# plot on same graph
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
lines(ordered_means_rescaled, sorted_means_rescaled, type = "l", lwd = 2, col = 3)
lines(ordered_vars_rescaled, sorted_vars_rescaled, type = "l", lwd = 2, col = 4)
abline(0,1,lty=2,lwd=2,col=adjustcolor(1,0.2))

# Add axes with original variable scales
axis(1, at = rescale(pretty(ordered_means)), labels = pretty(ordered_means), col.axis = 3)
axis(2, at = rescale(pretty(sorted_means)), labels = pretty(sorted_means), col.axis = 3)
axis(3, at = rescale(pretty(ordered_vars)), labels = pretty(ordered_vars), col.axis = 4)
axis(4, at = rescale(pretty(sorted_vars)), labels = pretty(sorted_vars), col.axis = 4)

title(main = paste0("Moments of Ordered Dirichlet Approximation (K = ", K, ")"), 
      xlab = "Approximation (truncated stick breaking)", 
      ylab = "True Constraint (truncated stick breaking)", line = 3)
legend("topleft", legend = c("Means", "Variances", "1-to-1 line"), 
       col = c(3, 4, adjustcolor(1, 0.2)), lwd = 2, lty = c(1,1,2))

#### inverse approach ####
n_samples <- 1E6
K <- 20
ks <- 1:(K-1)
alpha <- rep(1, K)
sorted_samples <- t(apply(extraDistr::rdirichlet(n_samples, alpha), 1, sort, decreasing = T))
sorted_prop_samples_all <- lapply(ks, function(k){
  
  #total mass available
  if(k == 1){
    remaining_mass <- 1
  } else {
    remaining_mass <- 1-rowSums(sorted_samples[,1:(k-1), drop = F])  
  }
  
  
  #get lower bound
  if(k == K){
    lb <- 0
  } else {
    lb <- remaining_mass/(K-k+1)
  }
  
  #get upper bound
  if(k == 1){
    ub <- 1
  } else {
    ubs <- cbind(sorted_samples[,k-1], remaining_mass)
    ub <- ubs[cbind(1:n_samples, (ubs[,1] > ubs[,2]) + 1)]
  }
  
  #rescale to valid interval to get a proportion
  prop <- (sorted_samples[,k] - lb) / (ub - lb)
  prop_rm <- sorted_samples[,k] / remaining_mass
  
  #find quantile in Beta for these
  shapes <- c(alpha[k], sum(alpha[(k+1):K]))
  qs <- pbeta(prop, shapes[1], shapes[2])
  
  #find quantile in truncated Beta for these
  qs_rm <- pbeta(prop_rm, shapes[1], shapes[2])
  lb_rm <- lb / remaining_mass
  ub_rm <- ub / remaining_mass
  betaq_lb <- pbeta(lb_rm, shapes[1], shapes[2])
  betaq_ub <- pbeta(ub_rm, shapes[1], shapes[2])
  qs_rm <- (qs_rm - betaq_lb) / (betaq_ub - betaq_lb)
  
  #return output
  out <- list(prop = prop,
              prop_rm = prop_rm,
              qs = qs,
              qs_rm = qs_rm)
  return(out)
  
})

extract_lol <- function(l, vn){
  lapply(l, function(x) x[[vn]])
}

vari <- 4
variable <- c("prop", "prop_rm", "qs", "qs_rm")[vari]
sorted_prop_samples <- do.call(cbind, extract_lol(sorted_prop_samples_all, variable))

#fit Beta distributions to these
mom_beta <- function(x) {
  mean_x <- mean(x)
  var_x <- var(x)
  alpha <- mean_x * ((mean_x * (1 - mean_x) / var_x) - 1)
  beta <- (1 - mean_x) * ((mean_x * (1 - mean_x) / var_x) - 1)
  return(list(estimate = c(shape1 = alpha, shape2 = beta)))
}

bfits <- lapply(ks, function(k){
  cat(paste0(k, " "))
  # init <- as.list(c(shape1 = mean(sorted_prop_samples[,i]),
  #                   shape2 = 1 - mean(sorted_prop_samples[,i])) * sum(alpha[i:K]))
  # init <- list(shape1=1, shape2=1)
  # bfit <- MASS::fitdistr(sorted_prop_samples[,k], "beta", start = init)
  bfit <- mom_beta(sorted_prop_samples[,k])
  return(bfit)
})
bfit_par <- data.frame(do.call(rbind, lapply(bfits, function(foo) foo$estimate)))
plot(ks, bfit_par$shape1, type = "l")
plot(ks, bfit_par$shape2, type = "l")
plot(ks, bfit_par$shape1 + bfit_par$shape2, type = "l")
plot(ks, bfit_par$shape1 / (bfit_par$shape1 + bfit_par$shape2), type = "l")
plot(bfit_par, type = "l")
text(x = bfit_par$shape1, y = bfit_par$shape2, labels = ks, col = 2)

#### Stan approach ####
try_Stan <- F
if(try_Stan){
  
library(cmdstanr)
library(posterior)
mod_path <- "~/scripts/minor_scripts/postdoc/ordered_simplex.stan"
mod <- cmdstan_model(mod_path)
dat <- list(K = 10)
fit <- mod$sample(dat = dat, chains = 4, parallel_chains = 4)

summ <- fit$summary()
print(summ[order(summ$ess_bulk), c("variable", 
                                   "rhat", "ess_bulk", "ess_tail", 
                                   "mean", "sd")])
print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                               "rhat", "ess_bulk", "ess_tail", 
                                               "mean", "sd")])
samps <- data.frame(as_draws_df(fit$draws()))
hist(samps$x.1.)
hist(samps$y.1., add = T, col = adjustcolor(2, 0.5))
hist(samps$x.10.)
hist(samps$y.10., add = T, col = adjustcolor(2, 0.5))


}