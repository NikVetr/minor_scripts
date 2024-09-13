library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

m <- 1000
x <- rnorm(m)
y <- rbeta(1,2,2)
e <- rexp(1)
z <- sqrt(y) * x + e 

stan_model <- "
data {
  int<lower=1> m;
  vector[m] z;
}

parameters {
  real<lower=0, upper=1> y;
  real<lower=0> e;
}

transformed parameters {
  vector[m] x = (z - e) / sqrt(y);
}

model {
  x ~ std_normal();
  y ~ beta(2, 2);
  e ~ exponential(1);

  // jacobian adjustment
  target += -(m / 2.0) * log(y);  // Jacobian adjustment for transformation.

}
"


stan_model_err <- "
data {
  int<lower=1> m;
  vector[m] z;
  real<lower=0> z_err;
}

parameters {
  vector[m] x;
  real<lower=0, upper=1> y;
  real<lower=0> e;
}

transformed parameters {
  vector[m] true_z = sqrt(y) * x + e;
}

model {
  x ~ std_normal();
  y ~ beta(2, 2);
  e ~ exponential(1);

  z ~ normal(true_z, z_err);
}
"

mod <- cmdstan_model(write_stan_file(stan_model))
dat <- list(m=m, z=z)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, 
                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                  refresh = 100, init = 0.1, thin = 2)
# summ <- fit$summary()
# summ

mod_err <- cmdstan_model(write_stan_file(stan_model_err))
dat_err <- list(m=m, z=z, z_err=0.01)
fit_err <- mod_err$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat_err, 
                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                  refresh = 100, init = 0.1, thin = 2)
# summ_err <- fit_err$summary()
# summ_err

#correct model
samps <- as.data.frame(as_draws_df(fit$draws()))
samps <- samps[,!(colnames(samps) %in% c(".chain", ".iteration", ".draw"))]

x.samps <- samps[,paste0("x[", 1:m, "]")]
y.samps <- samps[,paste0("y")]
e.samps <- samps[,paste0("e")]

xm <- apply(x.samps, 2, mean)
ym <- mean(y.samp)
em <- mean(e.samp)

#95% CI
CI_prob <- 0.95
xCI <- apply(x.samps, 2, quantile, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))
yCI <- quantile(y.samps, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))
eCI <- quantile(e.samps, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))


plot(x, xm, ylim = range(xCI))
abline(0,1)
for(i in 1:m){
  segments(x0 = x[i], x1 = x[i], y0 = xCI[1,i], y1 = xCI[2,i])
}

#incorrect model
samps_err <- as.data.frame(as_draws_df(fit_err$draws()))
samps_err <- samps_err[,!(colnames(samps_err) %in% c(".chain", ".iteration", ".draw"))]

x.samps_err <- samps_err[,paste0("x[", 1:m, "]")]
y.samps_err <- samps_err[,paste0("y")]
e.samps_err <- samps_err[,paste0("e")]

xm_err <- apply(x.samps_err, 2, mean)
ym_err <- mean(y.samps_err)
em_err <- mean(e.samps_err)

#95% CI
CI_prob <- 0.95
xCI_err <- apply(x.samps_err, 2, quantile, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))
yCI_err <- quantile(y.samps_err, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))
eCI_err <- quantile(e.samps_err, prob = c((1-CI_prob)/2, CI_prob + (1-CI_prob)/2))


plot(x, xm_err, ylim = range(xCI_err))
abline(0,1)
for(i in 1:m){
  segments(x0 = x[i], x1 = x[i], y0 = xCI_err[1,i], y1 = xCI_err[2,i])
}

#compare point estimates
(ym - y) / y * 100
yCI
y
(em - e) / e * 100
eCI
e

(ym_err - y) / y * 100
yCI_err
y
(em_err - e) / e * 100
eCI_err
e
