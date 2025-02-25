# Load necessary libraries
library(cmdstanr)
library(posterior)  # For working with draws and computing ESS

# data simulation
n_g <- 20
n_x_per_g <- sample(c(1,10), n_g, replace = T)
g <- rep(1:n_g, n_x_per_g)
n_x <- sum(n_x_per_g)
mu <- rnorm(n_g)
sigma <- exp(rnorm(n_g))
x <- rnorm(n_x, mu[g], sigma[g])
dat <- list(n_x = n_x,
            n_g = n_g,
            x = x,
            g = g)

#munge for mixed model
centered <- n_x_per_g == 50
n_x_c <- sum(centered[g])
n_x_nc <- sum(!centered[g])

n_g_c <- sum(centered)
n_g_nc <- sum(!centered)

x_c <- x[centered[g]]
x_nc <- x[!centered[g]]

g_c <- g[centered[g]]
g_nc <- g[!centered[g]]

g_c <- match(g_c, sort(unique(g_c)))
g_nc <- match(g_nc, sort(unique(g_nc)))

dat_mixed <- list(
  n_x_c = n_x_c,
  n_x_nc = n_x_nc,
  n_g_c = n_g_c,
  n_g_nc = n_g_nc,
  x_c = x_c,
  x_nc = x_nc,
  g_c = g_c,
  g_nc = g_nc
)

# Define the Stan models as strings
centered_model <- "
data {
  int<lower=1> n_x; //number of observations
  int<lower=1> n_g; //number of groups
  vector[n_x] x;  //observations
  array[n_x] int<lower=1, upper=n_g> g;  //group index of observations
}
parameters {
  vector[n_g] mu;
  vector[n_g] log_sigma;
  
  real mean_mu;
  real<lower=0> sd_mu;
  
  real mean_log_sigma;
  real<lower=0> sd_log_sigma;
}
transformed parameters {
  vector[n_g] sigma = exp(log_sigma);
}
model {
  //hyperpriors
  mean_mu ~ std_normal();
  sd_mu ~ std_normal();
  mean_log_sigma ~ std_normal();
  sd_log_sigma ~ std_normal();
  
  //priors
  mu ~ normal(mean_mu, sd_mu);
  log_sigma ~ normal(mean_log_sigma, sd_log_sigma);
  
  //likelihood
  x ~ normal(mu[g], sigma[g]);     // Prior on x given y (centered)
}
"

noncentered_model <- "
data {
  int<lower=1> n_x; //number of observations
  int<lower=1> n_g; //number of groups
  vector[n_x] x;  //observations
  array[n_x] int<lower=1, upper=n_g> g;  //group index of observations
}
parameters {
  vector[n_g] raw_mu;
  vector[n_g] raw_log_sigma;
  
  real mean_mu;
  real<lower=0> sd_mu;
  
  real mean_log_sigma;
  real<lower=0> sd_log_sigma;
}
transformed parameters {
  vector[n_g] mu = mean_mu + raw_mu * sd_mu;
  vector[n_g] log_sigma = mean_log_sigma + raw_log_sigma * sd_log_sigma;
  vector[n_g] sigma = exp(log_sigma);
}
model {
  //hyperpriors
  mean_mu ~ std_normal();
  sd_mu ~ std_normal();
  mean_log_sigma ~ std_normal();
  sd_log_sigma ~ std_normal();
  
  //priors
  raw_mu ~ std_normal();
  raw_log_sigma ~ std_normal();
  
  //likelihood
  x ~ normal(mu[g], sigma[g]);     // Prior on x given y (centered)
}
"

mixed_model <- "
data {
  int<lower=1> n_x_c; //number of observations
  int<lower=1> n_g_c; //number of groups
  vector[n_x_c] x_c;  //observations
  array[n_x_c] int<lower=1, upper=n_g_c> g_c;  //group index of observations
  
  int<lower=1> n_x_nc; //number of observations
  int<lower=1> n_g_nc; //number of groups
  vector[n_x_nc] x_nc;  //observations
  array[n_x_nc] int<lower=1, upper=n_g_nc> g_nc;  //group index of observations
  
}
parameters {
  vector[n_g_c] mu_c;
  vector[n_g_c] log_sigma_c;
  
  vector[n_g_nc] raw_mu_nc;
  vector[n_g_nc] raw_log_sigma_nc;
  
  real mean_mu;
  real<lower=0> sd_mu;
  
  real mean_log_sigma;
  real<lower=0> sd_log_sigma;
}
transformed parameters {
  vector[n_g_nc] mu_nc = mean_mu + raw_mu_nc * sd_mu;
  vector[n_g_nc] log_sigma_nc = mean_log_sigma + raw_log_sigma_nc * sd_log_sigma;
  vector[n_g_nc] sigma_nc = exp(log_sigma_nc);
  vector[n_g_c] sigma_c = exp(log_sigma_c);
}
model {
  //hyperpriors
  mean_mu ~ std_normal();
  sd_mu ~ std_normal();
  mean_log_sigma ~ std_normal();
  sd_log_sigma ~ std_normal();
  
  //priors
  raw_mu_nc ~ std_normal();
  raw_log_sigma_nc ~ std_normal();
  mu_c ~ normal(mean_mu, sd_mu);
  log_sigma_c ~ normal(mean_log_sigma, sd_log_sigma);
  
  //likelihood
  x_c ~ normal(mu_c[g_c], sigma_c[g_c]);     // Prior on x given y (centered)
  x_nc ~ normal(mu_nc[g_nc], sigma_nc[g_nc]);     // Prior on x given y (centered)
}
"


# Compile the Stan models
centered_mod <- cmdstan_model(write_stan_file(centered_model))
noncentered_mod <- cmdstan_model(write_stan_file(noncentered_model))
mixed_mod <- cmdstan_model(write_stan_file(mixed_model))

# Fit models
niter <- 1E3
nsamp <- 1E3
nref <- 20
fit_centered <- centered_mod$sample(
  data = dat, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

fit_noncentered <- noncentered_mod$sample(
  data = dat, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

fit_mixed <- mixed_mod$sample(
  data = dat_mixed, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

#compute mcmc diagnostics
summ_centered <- fit_centered$summary()
summ_centered[order(summ_centered$ess_bulk),]
summ_centered[order(summ_centered$rhat, decreasing = T),]

summ_noncentered <- fit_noncentered$summary()
summ_noncentered[order(summ_noncentered$ess_bulk),]
summ_noncentered[order(summ_noncentered$rhat, decreasing = T),]

summ_mixed <- fit_mixed$summary()
summ_mixed[order(summ_mixed$ess_bulk),]
summ_mixed[order(summ_mixed$rhat, decreasing = T),]

plot(summ_noncentered$ess_bulk[grepl(pattern = "mu\\[", summ_noncentered$variable)][1:n_g],
     n_x_per_g)
plot(summ_centered$ess_bulk[grepl(pattern = "mu\\[", summ_centered$variable)][1:n_g],
     n_x_per_g)

# Extract draws
samps_centered <- data.frame(as_draws_df(fit_centered$draws()))
samps_noncentered <- data.frame(as_draws_df(fit_noncentered$draws()))
samps_mixed <- data.frame(as_draws_df(fit_mixed$draws()))

mu_centered <- data.frame(as_draws_df(fit_centered$draws("mu")))[,-c((n_g+1):(n_g+4))]
mu_noncentered <- data.frame(as_draws_df(fit_noncentered$draws("mu")))[,-c((n_g+1):(n_g+4))]

sigma_centered <- data.frame(as_draws_df(fit_centered$draws("sigma")))[,-c((n_g+1):(n_g+4))]
sigma_noncentered <- data.frame(as_draws_df(fit_noncentered$draws("sigma")))[,-c((n_g+1):(n_g+4))]



