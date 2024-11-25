library(cmdstanr)
library(posterior)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#does iterative EM-flavored EB converge on Full Bayes MCMC?
#or maybe I can call this "maximum marginal likelihood"

#### functions ####
logit <- function(p) log(p/(1-p))
inv_logit <- function(x) exp(x)/(1+exp(x))

#### simulate data ####

#specify data structure
total <- 500
n_obs <- 50
n_groups_1 <- 5
n_groups_2 <- 10
n <- n_obs * n_groups_1 * n_groups_2
n_groups_3 <- 10 #extra grouping

#specify hyperparameters
sd_a_1 <- 0.8
sd_a_2 <- 0.5
sd_b_1 <- 0.7
sd_b_2 <- 1.2
sd_c_1 <- 0.9
sd_c_2 <- 0.6
sd_g <- 1

#sample parameters
a_1 <- rnorm(n_groups_1, sd = sd_a_1)
a_2 <- do.call(rbind, lapply(a_1, function(i) rnorm(n_groups_2, mean = i, sd = sd_a_2)))
b_1 <- rnorm(n_groups_1, sd = sd_b_1)
b_2 <- do.call(rbind, lapply(b_1, function(i) rnorm(n_groups_2, mean = i, sd = sd_b_2)))
c_1 <- rnorm(n_groups_1, sd = sd_c_1)
c_2 <- do.call(rbind, lapply(c_1, function(i) rnorm(n_groups_2, mean = i, sd = sd_c_2)))
g <- rnorm(n_groups_3, sd = sd_g)

#specify predictors (can do discrete predictors but need to vary group assignment)
sd_x <- 3
x_b <- rnorm(n) * sd_x
x_c <- rnorm(n) * sd_x

#specify group assignment
base_groups <- expand.grid(group_1 = 1:n_groups_1, 
                           group_2 = 1:n_groups_2, 
                           group_3 = 1:n_groups_3)
d <- data.frame(base_groups, 
                x_b = x_b, 
                x_c = x_c)
group_vec_ref <- matrix(1:(n_groups_1 * n_groups_2), nrow = n_groups_1, ncol = n_groups_2)
d$group_1x2 <- group_vec_ref[cbind(d$group_1, d$group_2)]

#simulate data
d$logit_p <- a_2[as.matrix(d[,c("group_1", "group_2")])] +
             b_2[as.matrix(d[,c("group_1", "group_2")])] * d$x_b +
             c_2[as.matrix(d[,c("group_1", "group_2")])] * d$x_c +
             g[d$group_3]
d$p <- inv_logit(d$logit_p)
d$count <- rbinom(n = n, size = total, prob = d$p)
d$total <- total

hist(d$count / d$total)

#### specify Stan models ####

# full Bayesian model
stan_program_full <- "
data {
  int<lower=1> n;
  int<lower=1> n_obs;
  int<lower=1> n_groups_1;
  int<lower=1> n_groups_2;
  int<lower=1> n_groups_3;
  array[n] int<lower=0> count;
  array[n] int<lower=1> total;
  vector[n] x_b;
  vector[n] x_c;
  array[n] int<lower=1, upper=n_groups_1> group_1;
  array[n] int<lower=1, upper=n_groups_2> group_2;
  array[n] int<lower=1, upper=n_groups_1 * n_groups_2> group_1x2;
  array[n] int<lower=1, upper=n_groups_3> group_3;
}

parameters {
  // location parameters and hyperparameters
  real a;
  vector[n_groups_1] a_1;
  vector[n_groups_1 * n_groups_2] a_2;
  vector[n_groups_1] b_1;
  vector[n_groups_1 * n_groups_2] b_2;
  vector[n_groups_1] c_1;
  vector[n_groups_1 * n_groups_2] c_2;
  vector[n_groups_3] g;
  
  // scale hyperparameters
  real<lower=0> sd_a_1;
  real<lower=0> sd_a_2;
  real<lower=0> sd_b_1;
  real<lower=0> sd_b_2;
  real<lower=0> sd_c_1;
  real<lower=0> sd_c_2;
  real<lower=0> sd_g;
}

model {
  // priors
  a ~ std_normal();
  a_1 ~ std_normal();
  a_2 ~ std_normal();
  b_1 ~ std_normal();
  b_2 ~ std_normal();
  c_1 ~ std_normal();
  c_2 ~ std_normal();
  g ~ std_normal();
  
  sd_a_1 ~ std_normal();
  sd_a_2 ~ std_normal();
  sd_b_1 ~ std_normal();
  sd_b_2 ~ std_normal();
  sd_c_1 ~ std_normal();
  sd_c_2 ~ std_normal();
  sd_g ~ std_normal();
  
  // model
  vector[n] logit_p = a + a_1[group_1] * sd_a_1 + a_2[group_1x2] * sd_a_2 +
                      (b_1[group_1] * sd_b_1 + b_2[group_1x2] * sd_b_2) .* x_b +
                      (c_1[group_1] * sd_c_1 + c_2[group_1x2] * sd_c_2) .* x_c +
                      g[group_3] * sd_g;
  count ~ binomial_logit(total, logit_p);
}
"

stan_program_fixed <- "
data {
  int<lower=1> n;
  int<lower=1> n_groups_3;
  array[n] int<lower=0> count;
  array[n] int<lower=1> total;
  real<lower=0> sd_a;
  real<lower=0> sd_b;
  real<lower=0> sd_c;
  vector[n] x_b;
  vector[n] x_c;
  array[n] int<lower=1, upper=n_groups_3> group_3;
}

parameters {
  // location parameters and hyperparameters
  real a;
  real b;
  real c;
  vector[n_groups_3] g;
  
  // scale hyperparameters
  real<lower=0> sd_g;
}

model {
  // priors
  a ~ std_normal();
  b ~ std_normal();
  c ~ std_normal();
  g ~ std_normal();
  sd_g ~ std_normal();
  
  // model
  vector[n] logit_p = a * sd_a +
                      b * sd_b * x_b +
                      c * sd_c * x_c +
                      g[group_3] * sd_g;
  count ~ binomial_logit(total, logit_p);
}
"

stan_program_errprop <- "
data {
  int<lower=1> n;
  vector[n] x_est;
  vector<lower=0>[n] x_err_sd;
}
parameters {
  real x_mu;
  real<lower=0> x_sigma;
  vector[n] x_raw;
}
transformed parameters {
  vector[n] x = x_raw * x_sigma + x_mu;
}
model {
  x_mu ~ std_normal();
  x_sigma ~ std_normal();
  x_raw ~ std_normal();
  x_est ~ normal(x, x_err_sd);
}
"

stan_program_errprop_multivariate <- "
data {
  int<lower=1> n;                                    // Number of observations
  int<lower=1> k;                                    // Dimension of each observation
  array[n] vector[k] x_est;                          // Observed estimates
  array[n] cholesky_factor_cov[k] x_err_L;           // Cholesky factors of error covariance matrices
}

parameters {
  vector[k] x_mu;                                    // Mean vector
  vector<lower=0>[k] x_sigma;                        // Scale parameters (standard deviations)
  cholesky_factor_corr[k] L_Omega;                   // Cholesky factor of the correlation matrix
  array[n] vector[k] x_raw;                          // Latent variables
}

transformed parameters {
  array[n] vector[k] x;                              // Transformed latent variables
  cholesky_factor_cov[k] x_sigma_L;                  // Cholesky factor of the covariance matrix
  x_sigma_L = diag_pre_multiply(x_sigma, L_Omega);
  for (i in 1:n)
    x[i] = x_mu + x_sigma_L * x_raw[i];
}

model {
  // Priors
  
  x_mu ~ normal(0, 1);                               // Prior for the mean vector
  x_sigma ~ normal(0, 1);                            // Prior for the scale parameters
  L_Omega ~ lkj_corr_cholesky(1);                    // Prior for the correlation matrix
  
  for (i in 1:n){
    x_raw[i] ~ std_normal();                         // Standard normal prior for latent variables
  }
  
  // Likelihood
  for (i in 1:n)
    x_est[i] ~ multi_normal_cholesky(x[i], x_err_L[i]);
}
"


stan_program_errprop_multivariate_empirical <- "
data {
  int<lower=1> n;                                    // Number of observations
  int<lower=1> k;                                    // Dimension of each observation
  array[n] vector[k] x_est;                          // Observed estimates
  array[n] cholesky_factor_cov[k] x_err_L;           // Cholesky factors of error covariance matrices
  vector[k] mean_x_mu;                               // Empirical hyperprior for x_mu (mean)
  cholesky_factor_cov[k] L_x_mu;                                  // Empirical hyperprior for x_mu (covariance)
}

parameters {
  vector[k] x_mu;                                    // Mean vector
  vector<lower=0>[k] x_sigma;                        // Scale parameters (standard deviations)
  cholesky_factor_corr[k] L_Omega;                   // Cholesky factor of the correlation matrix
  array[n] vector[k] x_raw;                          // Latent variables
}

transformed parameters {
  array[n] vector[k] x;                              // Transformed latent variables
  cholesky_factor_cov[k] x_sigma_L;                  // Cholesky factor of the covariance matrix
  x_sigma_L = diag_pre_multiply(x_sigma, L_Omega);
  for (i in 1:n)
    x[i] = x_mu + x_sigma_L * x_raw[i];
}

model {
  
  // Priors
  x_mu ~ multi_normal_cholesky(mean_x_mu, L_x_mu);   // Prior for the mean vector
  x_sigma ~ normal(0, 1);                            // Prior for the scale parameters
  L_Omega ~ lkj_corr_cholesky(1);                    // Prior for the correlation matrix
  
  for (i in 1:n){
    x_raw[i] ~ std_normal();                         // Standard normal prior for latent variables
  }
  
  // Likelihood
  for (i in 1:n)
    x_est[i] ~ multi_normal_cholesky(x[i], x_err_L[i]);
}
"

stan_program_errprop_multilevel <- "
data {
  int<lower=1> n;
  vector[n] x_est;
  vector<lower=0>[n] x_err_sd;
  int<lower=1> n_g;
  array[n] int<lower=0, upper=n_g> g;
}
parameters {
  // first level
  real x1_mu;
  real<lower=0> x1_sigma;
  vector[n_g] x1_raw;
  
  // second level
  real<lower=0> x2_sigma;
  vector[n] x2_raw;
}
transformed parameters {
  vector[n] x =  x1_mu + x1_raw[n_g] * x1_sigma + x2_raw * x2_sigma;
}
model {
  x1_mu ~ std_normal();
  x1_sigma ~ std_normal();
  x1_raw ~ std_normal();
  x2_sigma ~ std_normal();
  x2_raw ~ std_normal();
  x_est ~ normal(x, x_err_sd);
}
"

mod_full <- cmdstan_model(write_stan_file(stan_program_full))
mod_fixed <- cmdstan_model(write_stan_file(stan_program_fixed))
mod_errprop <- cmdstan_model(write_stan_file(stan_program_errprop))
mod_errprop_multilevel <- cmdstan_model(write_stan_file(stan_program_errprop_multilevel))
mod_errprop_multivariate <- cmdstan_model(write_stan_file(stan_program_errprop_multivariate))
mod_errprop_multivariate_empirical <- cmdstan_model(write_stan_file(stan_program_errprop_multivariate_empirical))

#### fit main models ####
fixed_scale_priors <- 1

# put data in list for Stan
dat_full <- list(
  n = n,
  n_obs = n_obs,
  n_groups_1 = n_groups_1,
  n_groups_2 = n_groups_2,
  n_groups_3 = n_groups_3,
  count = d$count,
  total = d$total,
  x_b = d$x_b,
  x_c = d$x_c,
  group_1 = d$group_1,
  group_2 = d$group_2,
  group_1x2 = d$group_1x2,
  group_3 = d$group_3
)
  
# specify mcmc parameters
nchains <- 4
niter <- 2E3
adapt_delta <- 0.85
max_treedepth <- 10
thin <- 1
init_mcmc <- 1

#full model
fit_full <- mod_full$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                              data = dat_full, parallel_chains = nchains, adapt_delta = adapt_delta,
                              refresh = 10, max_treedepth = max_treedepth,
                              thin = thin, init = init_mcmc)

#examine mcmc diagnostics
summ_full <- fit_full$summary()
summ_full[order(summ_full$ess_bulk),]
summ_full[order(summ_full$rhat, decreasing = T),]
samps_full <- data.frame(as_draws_df(fit_full$draws()))

#some super quick checks
apply(samps_full[,grepl("sd_", colnames(samps_full))], 2, mean)
utri <- function(x) x[upper.tri(x)]
utri(cor(samps_full[,grepl("sd_", colnames(samps_full))]))

#one group at a time
dats_fixed <- lapply(1:length(group_vec_ref), function(i){
  dsub <- d[d$group_1x2 == i,]
  dat_sub <- list(
    n = nrow(dsub),
    n_groups_3 = n_groups_3,
    count = dsub$count,
    total = dsub$total,
    x_b = dsub$x_b,
    x_c = dsub$x_c,
    group_3 = dsub$group_3,
    sd_a = fixed_scale_priors,
    sd_b = fixed_scale_priors,
    sd_c = fixed_scale_priors
  )
  return(dat_sub)
})

fits_fixed <- lapply(1:length(group_vec_ref), function(i){
  cat(paste0("(", i, "/", length(group_vec_ref),") "))
  sink(tempfile())
  fit_sub <- mod_fixed$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                              data = dats_fixed[[i]], parallel_chains = nchains, adapt_delta = adapt_delta,
                              refresh = 1E3, max_treedepth = max_treedepth,
                              thin = thin, init = init_mcmc)
  sink()
  return(fit_sub)
})

summs_fixed <- do.call(rbind, lapply(1:length(group_vec_ref), function(i){
  cat(paste0("(", i, "/", length(group_vec_ref),") "))
  summ_sub <- fits_fixed[[i]]$summary()
  c(ess_bulk = min(summ_sub$ess_bulk), rhat = max(summ_sub$rhat))
}))
summs_fixed

samps_fixed <- lapply(1:length(group_vec_ref), function(i){
  cat(paste0("(", i, "/", length(group_vec_ref),") "))
  samps <- data.frame(as_draws_df(fits_fixed[[i]]$draws()))
  samps
})

moments_fixed <- do.call(rbind, lapply(1:length(group_vec_ref), function(i){
  cat(paste0("(", i, "/", length(group_vec_ref),") "))
  samps_sub <- samps_fixed[[i]][,c("a", "b", "c", "sd_g")]
  samps_sub$log_sd_g <- log(samps_sub$sd_g)
  moms <- data.frame(apply(samps_sub, 2, function(x) c(mean = mean(x),
                                                       sd = sd(x), 
                                                       skew = e1071::skewness(x),
                                                       kurtosis = e1071::kurtosis(x),
                                                       gr0 = mean(x>0)), 
                           simplify = F))
  moms$g <- i
  moms$group_1 <- which(group_vec_ref == i, arr.ind = T)[1]
  moms$group_2 <- which(group_vec_ref == i, arr.ind = T)[2]
  moms$type <- rownames(moms)
  return(moms)
}))

param_names <- setNames(c("a", "b", "c", "log_sd_g"), 
                        c("a", "b", "c", "log_sd_g"))
moments_fixed <- lapply(param_names, function(i){
  psub <- data.frame(val = moments_fixed[,i], 
                     type = moments_fixed$type, 
                     group_1 = moments_fixed$group_1, 
                     group_2 = moments_fixed$group_2, 
                     g = moments_fixed$g)
  psub <- split(psub, psub$type)
  return(psub)
})

#### fit error models ####
dats_errprop_multilevel <- lapply(param_names, function(i){
  list(
    n = nrow(moments_fixed[[i]]$mean),
    x_est = moments_fixed[[i]]$mean$val,
    x_err_sd = moments_fixed[[i]]$sd$val,
    n_g = n_groups_1,
    g = moments_fixed[[i]]$mean$group_1
  )
})

nchains <- 4
niter <- 5E3
adapt_delta <- 0.95
max_treedepth <- 12
thin <- 1
init_mcmc <- 1
fits_errprop_multilevel <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  
  sink(tempfile())
  fit_errprop_multilevel <- mod_errprop_multilevel$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                                          data = dats_errprop_multilevel[[i]], parallel_chains = nchains, 
                                                          adapt_delta = adapt_delta,
                                                          refresh = 100, max_treedepth = max_treedepth,
                                                          thin = thin, init = init_mcmc)
  sink()
  
  return(fit_errprop_multilevel)
})

summs_errprop_multilevel <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  summ_sub <- fits_errprop_multilevel[[i]]$summary()
  c(ess_bulk = min(summ_sub$ess_bulk), rhat = max(summ_sub$rhat))
})

samps_errprop_multilevel <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  samps <- data.frame(as_draws_df(fits_errprop_multilevel[[i]]$draws()))
  samps[,grepl("_mu|_sigma", colnames(samps))]
})

#compare fits

#recovery of top level sigmas?
hist(samps_full$sd_a_1, freq = F)
hist(samps_errprop_multilevel$a$x1_sigma, freq = F); curve(dnorm(x) * 2, add = T)

hist(samps_full$sd_b_2, freq = F)
hist(samps_errprop_multilevel$b$x2_sigma, freq = F); curve(dnorm(x) * 2, add = T)

hist(samps_full$sd_c_1, freq = F)
hist(samps_errprop_multilevel$c$x1_sigma, freq = F); curve(dnorm(x) * 2, add = T)

#partitioned variance fails

#but combined does well!
par(mfrow = c(3,1))
breaks = 0:100/10
hist(sqrt(samps_full$sd_c_2^2 + samps_full$sd_c_1^2), freq = F, breaks = breaks)
hist(sqrt(samps_errprop_multilevel$c$x2_sigma^2 + samps_errprop_multilevel$c$x1_sigma),
     freq = F, add = T, breaks = breaks, col = adjustcolor(1,0.5))
mean(sqrt(samps_full$sd_c_2^2 + samps_full$sd_c_1^2))
mean(sqrt(samps_errprop_multilevel$c$x2_sigma^2 + samps_errprop_multilevel$c$x1_sigma))

hist(sqrt(samps_full$sd_a_2^2 + samps_full$sd_a_1^2), freq = F, breaks = breaks)
hist(sqrt(samps_errprop_multilevel$a$x2_sigma^2 + samps_errprop_multilevel$a$x1_sigma),
     freq = F, add = T, breaks = breaks, col = adjustcolor(1,0.5))
mean(sqrt(samps_full$sd_a_2^2 + samps_full$sd_a_1^2))
mean(sqrt(samps_errprop_multilevel$a$x2_sigma^2 + samps_errprop_multilevel$a$x1_sigma))

hist(sqrt(samps_full$sd_b_2^2 + samps_full$sd_b_1^2), freq = F, breaks = breaks)
hist(sqrt(samps_errprop_multilevel$b$x2_sigma^2 + samps_errprop_multilevel$b$x1_sigma),
     freq = F, add = T, breaks = breaks, col = adjustcolor(1,0.5))
mean(sqrt(samps_full$sd_b_2^2 + samps_full$sd_b_1^2))
mean(sqrt(samps_errprop_multilevel$b$x2_sigma^2 + samps_errprop_multilevel$b$x1_sigma))


#fit single level partitioned variance error prop models
dats_errprop <- lapply(param_names, function(i){
  submoms_mean <- split(moments_fixed[[i]]$mean, moments_fixed[[i]]$mean$group_1)
  submoms_sd <- split(moments_fixed[[i]]$sd, moments_fixed[[i]]$sd$group_1)
  lapply(1:n_groups_1, function(j){
    list(
      n = nrow(submoms_mean[[j]]),
      x_est = submoms_mean[[j]]$val,
      x_err_sd = submoms_sd[[j]]$val
    )  
  })
})

nchains <- 4
niter <- 2E3
adapt_delta <- 0.9
max_treedepth <- 10
thin <- 1
init_mcmc <- 1
fits_errprop <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  
  fits_errprop_param <- lapply(1:n_groups_1, function(j){
    cat(paste0("(", j, ") "))
    
    sink(tempfile())
    fit_errprop_param <- mod_errprop$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                                            data = dats_errprop[[i]][[j]], parallel_chains = nchains, 
                                                            adapt_delta = adapt_delta,
                                                            refresh = 100, max_treedepth = max_treedepth,
                                                            thin = thin, init = init_mcmc)
    sink()
    
    return(fit_errprop_param)
  })
  
  return(fits_errprop_param)
})

summs_errprop <- do.call(rbind, lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  summs_sub <- do.call(rbind, lapply(1:n_groups_1, function(j){
    summ_sub <- fits_errprop[[i]][[j]]$summary()
    c(ess_bulk = min(summ_sub$ess_bulk), rhat = max(summ_sub$rhat))
  }))
  summs_sub
}))

samps_errprop <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  samps_errprop_param <- lapply(setNames(1:n_groups_1, paste0("gr_", 1:n_groups_1)), function(j){
    samps <- data.frame(as_draws_df(fits_errprop[[i]][[j]]$draws()))
    samps[,grepl("_mu|_sigma", colnames(samps))]
  })
})

mean(samps_full$sd_a_2)
mean(samps_errprop$a$gr_1$x_sigma)
mean(samps_errprop$a$gr_2$x_sigma)
mean(samps_errprop$a$gr_3$x_sigma)
mean(samps_errprop$a$gr_4$x_sigma)
mean(samps_errprop$a$gr_5$x_sigma)

mean(samps_full$sd_b_2)
mean(samps_errprop$b$gr_1$x_sigma)
mean(samps_errprop$b$gr_2$x_sigma)
mean(samps_errprop$b$gr_3$x_sigma)
mean(samps_errprop$b$gr_4$x_sigma)

mean(samps_full$sd_c_2)
mean(samps_errprop$c$gr_1$x_sigma)
mean(samps_errprop$c$gr_2$x_sigma)
mean(samps_errprop$c$gr_3$x_sigma)
mean(samps_errprop$c$gr_4$x_sigma)

#let's do it again! up the chain


#fit single level partitioned variance error prop models
dats_errprop_recurs <- lapply(param_names, function(i){
  
  dtsub <- do.call(rbind, lapply(1:n_groups_1, function(j){
    data.frame(
      x_est = mean(samps_errprop[[i]][[j]]$x_mu),
      x_err_sd = sd(samps_errprop[[i]][[j]]$x_mu)
    )
  }))
  
  datsub <- list(
    n = nrow(dtsub),
    x_est = dtsub$x_est,
    x_err_sd = dtsub$x_err_sd
  )
  
  return(datsub)
  
})

nchains <- 4
niter <- 2E3
adapt_delta <- 0.9
max_treedepth <- 10
thin <- 1
init_mcmc <- 1
fits_errprop_recurs <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  
  sink(tempfile())
  fit_errprop_param <- mod_errprop$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                          data = dats_errprop_recurs[[i]], parallel_chains = nchains, 
                                          adapt_delta = adapt_delta,
                                          refresh = 100, max_treedepth = max_treedepth,
                                          thin = thin, init = init_mcmc)
  sink()
  
  return(fit_errprop_param)
})


samps_errprop_recurs <- lapply(param_names, function(i){
  cat(paste0("(", i, ") "))
  samps <- data.frame(as_draws_df(fits_errprop_recurs[[i]]$draws()))
  samps[,grepl("_mu|_sigma", colnames(samps))]
})

#compare to recovery of top-level params
par(mfrow = c(3,1), mar = c(5,5,3,3))
breaks = 0:100/10
hist(samps_full$sd_a_1, breaks = breaks, freq = F, ylim = c(0,2))
hist(samps_errprop_recurs$a$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))

hist(samps_full$sd_b_1, breaks = breaks, freq = F, ylim = c(0,2))
hist(samps_errprop_recurs$b$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))

hist(samps_full$sd_c_1, breaks = breaks, freq = F, ylim = c(0,2))
hist(samps_errprop_recurs$c$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))

#works often enough, but not perfectly!

# try multivariate error model

multiv_input <- lapply(1:n_groups_1, function(i){
  
  cat(paste0("(", i, "/", n_groups_1,") "))
  
  out <- data.frame(do.call(rbind, lapply(1:n_groups_2, function(j){
      
    #get sumstats
    samps_sub <- samps_fixed[[group_vec_ref[i,j]]][,c("a", "b", "c")]
    ss_cov <- cov(samps_sub)
    ss_cor <- cor(samps_sub)
    ss_L <- t(chol(ss_cov))
    ss_mu <- apply(samps_sub, 2, mean)
    
    return(list(
      ss_L = ss_L,
      ss_mu = ss_mu,
      g1i = i,
      g2i = j,
      cors = utri(ss_cor)
    ))
    
  })))
  
  return(out)
  
})

#quick check of corrs
hist(sapply(1:n_groups_1, function(i) unlist(multiv_input[[i]]$cors)), breaks = 100)

#prepare data for Stan
dats_errprop_multivariate <- lapply(1:n_groups_1, function(i){
  
  cat(paste0("(", i, "/", n_groups_1,") "))
  
  list(
    n = nrow(multiv_input[[i]]),
    k = length(multiv_input[[i]]$ss_mu[[1]]),
    x_est = multiv_input[[i]]$ss_mu,
    x_err_L = multiv_input[[i]]$ss_L
  )
  
})

#fit error propagation model
nchains <- 4
niter <- 5E3
adapt_delta <- 0.95
max_treedepth <- 12
thin <- 1
init_mcmc <- 1
fits_errprop_multivariate <- lapply(1:n_groups_1, function(i){
    cat(paste0("(", i, "/", n_groups_1,") "))
  
    sink(tempfile())
    fit <- mod_errprop_multivariate$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                                        data = dats_errprop_multivariate[[i]], parallel_chains = nchains, 
                                                        adapt_delta = adapt_delta,
                                                        refresh = 100, max_treedepth = max_treedepth,
                                                        thin = thin, init = init_mcmc)
    sink()
    
    return(fit)
    
})

#examine mcmc diagnostics
summs_errprop_multivariate <- lapply(1:n_groups_1, function(i){
  cat(paste0("(", i, "/", n_groups_1,") "))
  summ_sub <- fits_errprop_multivariate[[i]]$summary()
  c(ess_bulk = min(summ_sub$ess_bulk, na.rm = T), rhat = max(summ_sub$rhat, na.rm = T))
})
summs_errprop_multivariate

#extract samples and process
samps_errprop_multivariate <- lapply(1:n_groups_1, function(i){
    cat(paste0("(", i, "/", n_groups_1,") "))
    samps <- data.frame(as_draws_df(fits_errprop_multivariate[[i]]$draws()))
    samps[,grepl("_mu", colnames(samps))]
})

multiv_input_recurs <- lapply(1:n_groups_1, function(i){

  cat(paste0("(", i, "/", n_groups_1,") "))

  list(
    x_est = apply(samps_errprop_multivariate[[i]], 2, mean),
    x_err_L = t(chol(cov(samps_errprop_multivariate[[i]])))
  )

})

dat_errprop_multivariate_recurs <- list(
  x_est = lapply(multiv_input_recurs, function(foo) foo$x_est),
  x_err_L = lapply(multiv_input_recurs, function(foo) foo$x_err_L),
  n = length(multiv_input_recurs),
  k = length(multiv_input_recurs[[1]]$x_est)
)

fit_errprop_multivariate_recurs <- mod_errprop_multivariate$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                       data = dat_errprop_multivariate_recurs, parallel_chains = nchains, 
                                       adapt_delta = adapt_delta,
                                       refresh = 100, max_treedepth = max_treedepth,
                                       thin = thin, init = init_mcmc)

summ_errprop_multivariate_recurs <- fit_errprop_multivariate_recurs$summary()
summ_errprop_multivariate_recurs[order(summ_errprop_multivariate_recurs$ess_bulk),]
summ_errprop_multivariate_recurs[order(summ_errprop_multivariate_recurs$rhat, decreasing = T),]
samps_errprop_multivariate_recurs <- data.frame(as_draws_df(fit_errprop_multivariate_recurs$draws()))

#compare recovery of top-level params
par(mfrow = c(3,2), mar = c(5,5,3,3))
breaks = 0:40/10
hist(samps_full$sd_a_1, breaks = breaks, freq = F, ylim = c(0,3), main = "multivariate retrieval for <<a>> param")
hist(samps_errprop_multivariate_recurs$x_sigma.1., freq = F, add = T, breaks = breaks, 
     col = adjustcolor(4,0.5))

hist(samps_full$sd_a_1, breaks = breaks, freq = F, ylim = c(0,3), main = "univariate retrieval for <<a>> param")
hist(samps_errprop_recurs$a$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))

hist(samps_full$sd_b_1, breaks = breaks, freq = F, ylim = c(0,3), main = "multivariate retrieval for <<b>> param")
hist(samps_errprop_multivariate_recurs$x_sigma.2., freq = F, add = T, breaks = breaks, 
     col = adjustcolor(4,0.5))

hist(samps_full$sd_b_1, breaks = breaks, freq = F, ylim = c(0,3), main = "univariate retrieval for <<b>> param")
hist(samps_errprop_recurs$b$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))

hist(samps_full$sd_c_1, breaks = breaks, freq = F, ylim = c(0,3), main = "multivariate retrieval for <<c>> param")
hist(samps_errprop_multivariate_recurs$x_sigma.3., freq = F, add = T, breaks = breaks, 
     col = adjustcolor(4,0.5))

hist(samps_full$sd_c_1, breaks = breaks, freq = F, ylim = c(0,3), main = "univariate retrieval for <<c>> param")
hist(samps_errprop_recurs$c$x_sigma, freq = F, add = T, breaks = breaks, 
     col = adjustcolor(2,0.5))


#check full samples for correlation between params
# cor(samps_full[,c("sd_a_1", "sd_b_1", "sd_c_1")])
# cor_unscaled <- cor(cbind(
#   samps_full[,colnames(samps_full)[grepl("a_1.", colnames(samps_full))]],# * samps_full$sd_a_1,
#   samps_full[,colnames(samps_full)[grepl("b_1.", colnames(samps_full))]],# * samps_full$sd_b_1,
#   samps_full[,colnames(samps_full)[grepl("c_1.", colnames(samps_full))]]# * samps_full$sd_c_1
# ))
# cor_scaled <- cor(cbind(
#   samps_full[,colnames(samps_full)[grepl("a_1.", colnames(samps_full))]] * samps_full$sd_a_1,
#   samps_full[,colnames(samps_full)[grepl("b_1.", colnames(samps_full))]] * samps_full$sd_b_1,
#   samps_full[,colnames(samps_full)[grepl("c_1.", colnames(samps_full))]] * samps_full$sd_c_1
# ))
# 
# unscaled_cors <- which(abs(cor_unscaled) > 0.1, arr.ind = T)
# scaled_cors <- which(abs(cor_scaled) > 0.1, arr.ind = T)
# 
# unscaled_cors <- matrix(rownames(cor_unscaled)[unscaled_cors], nrow = nrow(unscaled_cors), ncol = 2)
# scaled_cors <- matrix(rownames(cor_scaled)[scaled_cors], nrow = nrow(scaled_cors), ncol = 2)
# 
# unscaled_cors <- unscaled_cors[apply(unscaled_cors, 1, function(x) strsplit(x[1], "_")[[1]][1] != strsplit(x[2], "_")[[1]][1]),]
# scaled_cors <- scaled_cors[apply(scaled_cors, 1, function(x) strsplit(x[1], "_")[[1]][1] != strsplit(x[2], "_")[[1]][1]),]

