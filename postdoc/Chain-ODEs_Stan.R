
# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
library(deSolve)      # for ODE simulation
library(cmdstanr)     # for fitting Stan model
library(posterior)    # for examining results

set.seed(123)

# ------------------------------------------------------------
# Define a 20-dimensional ODE system
# 
# We'll have a chain-like system:
# dy[1]/dt = -k[1]*y[1]
# dy[i]/dt = k[i-1]*y[i-1] - k[i]*y[i] for i = 2,...,20
#
# Parameters: k[1], ..., k[20]
# Initial conditions: y0[1] = 1, y0[2..20] = 0
# ------------------------------------------------------------

num_states <- 10
num_params <- num_states
times <- seq(0, 9, by=1)   # 10 timepoints: 0 through 9
n_time <- length(times)
n_reps <- 5                # 5 replicate measurements per timepoint

# True parameters for simulation
k_true <- runif(num_params, 0, 0.5)  # random rate constants
y0 <- rnorm(num_states)       # initial condition

ode_system <- function(t, y, parms) {
  k <- parms
  dy <- numeric(num_states)
  # first equation
  dy[1] <- -k[1]*y[1]
  # chain equations
  for (i in 2:num_states) {
    dy[i] <- k[i-1]*y[i-1] - k[i]*y[i]
  }
  return(list(dy))
}

# Simulate "true" trajectories (no noise) and add noise
out <- ode(y=y0, times=times, func=ode_system, parms=k_true)
true_states <- out[,-1]  # matrix: rows=times, cols=states
ts_sds <- apply(true_states, 2, sd)
true_sigma <- exp(mean(log(ts_sds)))
obs_data <- array(NA, dim=c(n_time, n_reps, num_states))
for (i in 1:n_time) {
  for (r in 1:n_reps) {
    obs_data[i,r,] <- rnorm(num_states, mean=true_states[i,], sd=true_sigma)
  }
}

#quick visualization
state_i <- 4
plot(true_states[,state_i], type = "l", ylim = range(obs_data[,,state_i]))
points(true_states[,state_i])
for(i in 1:n_reps){
  lines(obs_data[,i,state_i], col = 2)
}


# ------------------------------------------------------------
# Write Stan model to a file
# ------------------------------------------------------------
stan_model_code <- '
functions {
  // ODE system: 
  // dy[1]/dt = -theta[1]*y[1]
  // dy[i]/dt = theta[i-1]*y[i-1] - theta[i]*y[i], for i = 2,...,M
  vector ode_system(
    real t,
    vector y,
    vector theta,
    array[] real x_r,
    array[] int x_i
  ) {
    int M = x_i[1];
    vector[M] dy;
    dy[1] = -theta[1]*y[1];
    for (i in 2:M) {
      dy[i] = theta[i-1]*y[i-1] - theta[i]*y[i];
    }
    return dy;
  }
}
data {
  int<lower=1, upper=2> solver; //1 for stiff solver, 2 for non-stiff solver
  int<lower=1> N;           // number of timepoints
  int<lower=1> R;           // number of replicates
  int<lower=1> M;           // number of states
  real t0;                  // initial time
  array[N] real ts;         // observation times
  array[N, R, M] real y_obs;// observed data
  vector[M] y0;             // initial conditions as a vector
  array[0] real x_r;              // empty real data array
  array[1] int x_i;               // holds M
}
parameters {
  vector<lower=0>[M] k;     // rate constants
  real<lower=0> sigma;      // measurement noise
}
transformed parameters {
  array[N] vector[M] y_hat;
  if(solver == 1){
    y_hat = ode_bdf(ode_system, y0, t0 - 1E-9, ts, k, x_r, x_i);
  }
  if(solver == 2){
    y_hat = ode_rk45(ode_system, y0, t0 - 1E-9, ts, k, x_r, x_i);
  }
}
model {
  // Priors
  k ~ normal(0, 0.5);
  sigma ~ normal(0.1, 0.05);

  // Likelihood: each replicate at each timepoint for each analyte
  for (n in 1:N) {
    for (r in 1:R) {
      y_obs[n, r] ~ normal(y_hat[n], sigma);
    }
  }
}
generated quantities {
  array[N, R, M] real y_rep;
  for (n in 1:N) {
    for (r in 1:R) {
      for (m in 1:M) {
        y_rep[n, r, m] = normal_rng(y_hat[n, m], sigma);
      }
    }
  }
}
'

# ------------------------------------------------------------
# Fit the model using cmdstanr
# ------------------------------------------------------------
mod_chain <- cmdstan_model(write_stan_file(stan_model_code), 
                     cpp_options = list(stan_threads = TRUE))
dat <- list(
  solver = 1,
  N = n_time,
  R = n_reps,
  M = num_states,
  t0 = 0,
  ts = times,
  y_obs = obs_data,
  y0 = y0,
  x_r = numeric(0),      # empty real array
  x_i = c(num_states)    # int array holding M
)

fit <- mod$sample(
  refresh = 10,
  data=dat,
  seed=123,
  chains=4,
  parallel_chains=4,
  iter_warmup=1000,
  iter_sampling=1000,
  threads_per_chain = 1
)
fit_pf <- mod$pathfinder(
  refresh = 10,
  data=dat,
  seed=123,
  num_threads = 4
)

#extract results
summ <- fit$summary()
print(summ[order(summ$ess_bulk), c("variable", 
                                   "rhat", "ess_bulk", "ess_tail", 
                                   "mean", "sd")])
print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                               "rhat", "ess_bulk", "ess_tail", 
                                               "mean", "sd")])
samps <- as.data.frame(as_draws_df(fit$draws(variables = c("k", "sigma"))))
samps_pf <- as.data.frame(as_draws_df(fit_pf$draws(variables = c("k", "sigma"))))

# plot
par(mfrow = c(2,1))
hist(samps$sigma); abline(v=true_sigma, col = 2, lwd = 5)
k_samps <- samps[grepl("k\\[", colnames(samps))]
k_pmean <- apply(k_samps, 2, mean)
k_90CI <- apply(k_samps, 2, quantile, prob = c(0.05, 0.95))
k_ord <- order(k_true)
plot(k_pmean[k_ord], ylim = range(c(k_90CI, k_true)))
for(i in 1:num_params){
  lines(c(i, i), y = k_90CI[,k_ord[i]], lwd = 3)
}
points(k_true[k_ord], col = 2)

#compare to pathfinder
hist(samps_pf$sigma); abline(v=true_sigma, col = 2, lwd = 5)
k_samps_pf <- samps_pf[grepl("k\\[", colnames(samps_pf))]
k_pmean_pf <- apply(k_samps_pf, 2, mean)
k_90CI_pf <- apply(k_samps_pf, 2, quantile, prob = c(0.05, 0.95))
k_ord <- order(k_true)
plot(k_pmean_pf[k_ord], ylim = range(c(k_90CI_pf, k_true)))
for(i in 1:num_params){
  lines(c(i, i), y = k_90CI_pf[,k_ord[i]], lwd = 3)
}
points(k_true[k_ord], col = 2)
