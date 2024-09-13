library(tidyverse)
library(cmdstanr)
library(posterior)

# Generative model ####
# Model parameters (we will try to estimate these later)
mean_x <- 1
sd_x <- 0.2
mean_y <- 2
sd_y <- 0.4

# Generate data
n <- 100
# Non-observable data
x <- rnorm(n, mean = mean_x, sd = sd_x)
y <- rnorm(n, mean = mean_y, sd = sd_y)

# Observable data
w <- runif(n)
V <- w*exp(x) - (1-w)*exp(y)

# log((V - w * exp(x)) / -(1-w)) #- y
# x = log((V + (1-w) * exp(y)) / w)

# for this to be valid, (V + (1-w) * exp(y)) / w needs to be > 0
# so (V + (1-w) * exp(y)) / w > 0
# V + (1-w) * exp(y) > 0
# (1-w) * exp(y) > -V
# exp(y) > -V / (1-w)
lb_y <- log(-V / (1-w))
unconstrained_y <- is.na(lb_y)
lb_y[unconstrained_y] <- -1E3
dat <- list(n = n, V = V, w = w)

#check derivatives for jacobian adjustment
D(expression(log((V + (1 - w) * exp(y)) / w)), "y") #forward transform
D(expression(log( (w * exp(x) - V) / (1 - w) )), "x") #reverse transformation


#### normal model ####

stan_model <- "

functions {
  real vector_product(vector x) {
    real product = 1;
    for (n in 1:size(x)) {
      product *= x[n];
    }
    return product;
  }
}

data {
  int<lower=0> n;
  vector<lower=0,upper=1>[n] w;
  vector[n] V;
}

transformed data {
  vector[n] lb_y;
  for(i in 1:n){
    lb_y[i] = V[i] > 0 ? negative_infinity() : log(-V[i] / (1-w[i])); 
  }
}

parameters {
  real mean_x;
  real mean_y;
  real<lower=0> sd_x;
  real<lower=0> sd_y;
  vector<lower=lb_y>[n] y;
}
transformed parameters {
  vector[n] x = log( ( V + ( 1 - w ) .* exp( y ) ) ./ w );
}

model {
  // priors on hyperparameters
  mean_x ~ std_normal();
  sd_x ~ std_normal();
  mean_y ~ std_normal();
  sd_y ~ std_normal();
  
  // centered distributions on latent variables
  x ~ normal(mean_x, sd_x);
  y ~ normal(mean_y, sd_y);

  // Jacobian adjustment (these are equivalent up to a constant)
  //target += log( abs( vector_product ( exp(x) ./ (exp(x) - V ./ w) ) ) );
  target += log( abs( vector_product ( (1 - w) .* exp(y) ./ w ./ ( (V + (1 - w) .* exp(y)) ./ w) ) ) );
  

}

"

#### lognormal model ####

stan_model_lognormal <- "
data {
  int<lower=0> n;
  vector<lower=0,upper=1>[n] w;
  vector[n] V;
}

transformed data {
  vector[n] lb_y;
  for(i in 1:n){
    lb_y[i] = V[i] > 0 ? 0 : -V[i] / (1-w[i]); 
  }
}

parameters {
  real mean_x;
  real mean_y;
  real<lower=0> sd_x;
  real<lower=0> sd_y;
  vector<lower=lb_y>[n] y;
}

transformed parameters {
  vector[n] x = ( V + ( 1 - w ) .* y ) ./ w;
}

model {
  // priors on hyperparameters
  mean_x ~ std_normal();
  sd_x ~ std_normal();
  mean_y ~ std_normal();
  sd_y ~ std_normal();
  
  // centered distributions on latent variables
  x ~ lognormal(mean_x, sd_x);
  y ~ lognormal(mean_y, sd_y);

  // Jacobian adjustment not necessary
}

"

#### just x lognormal model ####

stan_model_lognormal_x <- "
data {
  int<lower=0> n;
  vector<lower=0,upper=1>[n] w;
  vector[n] V;
}

transformed data {
  vector[n] lb_y;
  for(i in 1:n){
    lb_y[i] = V[i] > 0 ? negative_infinity() : log(-V[i] / (1-w[i])); 
  }
}

parameters {
  real mean_x;
  real mean_y;
  real<lower=0> sd_x;
  real<lower=0> sd_y;
  vector<lower=lb_y>[n] y;
}

transformed parameters {
  vector[n] x = ( V + ( 1 - w ) .* exp(y) ) ./ w;
}

model {
  // priors on hyperparameters
  mean_x ~ std_normal();
  sd_x ~ std_normal();
  mean_y ~ std_normal();
  sd_y ~ std_normal();
  
  // centered distributions on latent variables
  x ~ lognormal(mean_x, sd_x);
  y ~ normal(mean_y, sd_y);

  // Jacobian adjustment
  // target += sum(y);

}

"

#### fits ####

mod <- cmdstan_model(write_stan_file(stan_model))
fit <- mod$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3,
                  data = dat, adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 10, init = 0.1)
summ <- fit$summary()
print(summ[order(summ$rhat, decreasing = T),])

#fit lognormal form of the model
mod_logn <- cmdstan_model(write_stan_file(stan_model_lognormal))
fit_logn <- mod_logn$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3,
                            data = dat, adapt_delta = 0.9, parallel_chains = 4,
                            refresh = 100, max_treedepth = 10, 
                            thin = 10, init = 0.1)
summ_logn <- fit_logn$summary()
print(summ_logn[order(summ_logn$rhat, decreasing = T),])

#fit lognormal_x form of the model
mod_logn_x <- cmdstan_model(write_stan_file(stan_model_lognormal_x))
fit_logn_x <- mod_logn_x$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3,
                            data = dat, adapt_delta = 0.9, parallel_chains = 4,
                            refresh = 100, max_treedepth = 10, 
                            thin = 10, init = 0.1)
summ_logn_x <- fit_logn_x$summary()
print(summ_logn[order(summ_logn$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws()))
samps_logn <- data.frame(as_draws_df(fit_logn$draws()))
samps_logn_x <- data.frame(as_draws_df(fit_logn_x$draws()))

#do some plotting

breaks = 0:300/100
hist(samps$mean_x, col = adjustcolor("red", 0.5), freq = F,
     xlim = range(samps$mean_x, samps_logn$mean_x), breaks = breaks)
hist(samps_logn$mean_x, add = T, col = adjustcolor("blue", 0.5), freq = F, breaks = breaks)
hist(samps_logn_x$mean_x, add = T, col = adjustcolor("orange", 0.5), freq = F, breaks = breaks)
abline(v = mean_x, col = 3, lwd = 5)

hist(samps$mean_y, col = adjustcolor("red", 0.5), freq = F, 
     xlim = range(samps$mean_y, samps_logn$mean_y), breaks = breaks)
hist(samps_logn$mean_y, add = T, col = adjustcolor("blue", 0.5), freq = F, breaks = breaks)
hist(samps_logn_x$mean_y, add = T, col = adjustcolor("orange", 0.5), freq = F, breaks = breaks)
abline(v = mean_y, col = 3, lwd = 5)

hist(samps$sd_x, col = adjustcolor("red", 0.5), freq = F,
     xlim = range(samps$sd_x, samps_logn$sd_x), breaks = breaks)
hist(samps_logn$sd_x, add = T, col = adjustcolor("blue", 0.5), freq = F, breaks = breaks)
hist(samps_logn_x$sd_x, add = T, col = adjustcolor("orange", 0.5), freq = F, breaks = breaks)
abline(v = sd_x, col = 3, lwd = 5)

hist(samps$sd_y, col = adjustcolor("red", 0.5), freq = F,
     xlim = range(samps$sd_y, samps_logn$sd_y), breaks = breaks)
hist(samps_logn$sd_y, add = T, col = adjustcolor("blue", 0.5), freq = F, breaks = breaks)
hist(samps_logn_x$sd_y, add = T, col = adjustcolor("orange", 0.5), freq = F, breaks = breaks)
abline(v = sd_y, col = 3, lwd = 5)

