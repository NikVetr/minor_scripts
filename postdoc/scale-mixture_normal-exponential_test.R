library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#first try with just generating a univariate random variable
n <- 1E5
x <- matrix(rnorm(n * 2), ncol = 2)
e <- rexp(n, 1)
l1 <- x * sqrt(2*e)
l2 <- e * (rbinom(n, 1, 0.5)*2-1)
breaks <- seq(min(c(l1,l2)), max(c(l1,l2)), length.out = 100)
hist(l1, freq = F, breaks = breaks, col = adjustcolor(4,0.5))
hist(l2, add = T, freq = F, breaks = breaks, col = adjustcolor(2,0.5))
ks.test(l1, l2)

#try with a bivariate random variable?
r <- 0.6
R <- diag(2) + r - diag(2) * r
L <- t(chol(R))
mx <- t(L %*% t(x))
ml <- t(L %*% t(x * sqrt(e)))

cov(mx)
cov(ml)

hist(ml[,1], breaks = breaks, col = adjustcolor(2,0.5), freq = F)
hist(mx[,1], breaks = breaks, col = adjustcolor(4,0.5), add = T, freq = F)

#test it in Stan
p <- 50
n <- 5
x_mu <- rnorm(p) * 1.4 + 3
sigma <- 1
x <- do.call(cbind, lapply(1:p, function(i) rnorm(n, x_mu, sigma)))

dat <- list(n=n, p=p, x=x)

laplace_model <- "
data{
    int<lower = 1> p; // number of parameters 
    int<lower = 1> n; // number of obs 
    matrix[n, p] x;
}

parameters{
  
  // MVN parameters
  vector[p] x_mu_raw;
  real x_mu_mu;
  real<lower=0> x_sigma;
  real<lower=0> sigma;
  
}

transformed parameters{
  vector[p] x_mu = x_mu_raw * x_sigma + x_mu_mu;
}

model{
  x_mu_raw ~ double_exponential(0, 1);
  x_mu_mu ~ std_normal();
  x_sigma ~ std_normal();
  sigma ~ std_normal();
  for(i in 1:n){
    x[i,] ~ normal(x_mu, sigma);
  }
}
"

normal_mixture_model <- "
data{
    int<lower = 1> p; // number of parameters 
    int<lower = 1> n; // number of obs 
    matrix[n, p] x;
}

parameters{
  
  // MVN parameters
  //real<lower=0> w;
  vector<lower=0>[p] w;
  vector[p] x_mu_raw;
  real x_mu_mu;
  real<lower=0> x_sigma;
  real<lower=0> sigma;
  
}

transformed parameters{
  vector[p] x_mu = x_mu_raw * x_sigma .* sqrt(w) + x_mu_mu;
}

model{
  w ~ exponential(1);
  x_mu_raw ~ std_normal();
  x_mu_mu ~ std_normal();
  x_sigma ~ std_normal();
  sigma ~ std_normal();
  for(i in 1:n){
    x[i,] ~ normal(x_mu, sigma);
  }
}
"

#fit model
mod <- cmdstan_model(write_stan_file(laplace_model))
mod_alt <- cmdstan_model(write_stan_file(normal_mixture_model))

fit <- mod$sample(chains = 4, iter_sampling = 1E4, iter_warmup = 1E4,
                  data = dat, adapt_delta = 0.85, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 2)
fit_alt <- mod_alt$sample(chains = 4, iter_sampling = 1E4, iter_warmup = 1E4,
                  data = dat, adapt_delta = 0.85, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 2)


#check convergence
# summ <- fit$summary()
# print(summ[order(summ$ess_bulk),])
# print(summ[order(summ$rhat, decreasing = T),])
# 
# summ_alt <- fit_alt$summary()
# print(summ_alt[order(summ_alt$ess_bulk),])
# print(summ_alt[order(summ_alt$rhat, decreasing = T),])


#### inspect posterior  ####

#extract samples and inspect
samps <- data.table(as_draws_df(fit$draws()))
samps_alt <- data.table(as_draws_df(fit_alt$draws()))
xmu_samps <- munge_samps("x_mu", subset_samps("x_mu", samps))
xmu_samps_alt <- munge_samps("x_mu", subset_samps("x_mu", samps_alt))
xmu_samps <- do.call(rbind, xmu_samps)
xmu_samps_alt <- do.call(rbind, xmu_samps_alt)

plot(apply(xmu_samps, 2, mean), apply(xmu_samps_alt, 2, mean))
abline(0,1)

#run ks tests
hist(sapply(1:p, function(i){
  ks.test(xmu_samps[,i], xmu_samps_alt[,i])$p.value
}), breaks = 0:20/20, xlab = "ks.test p-value", main = "")

#hmm weirdly distinct at the tails...
#ah it is bc I need p independent auxiliary variables (w), not just one scalar w
#fixing now...


