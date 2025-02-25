library(deSolve)
library(cmdstanr)
library(posterior)
library(rstan)
library(igraph)
library(RColorBrewer)
library(bridgesampling)
rstan_options(auto_write = TRUE)   # Allow auto-writing of compiled models
options(mc.cores = parallel::detectCores())  # Use multiple cores

seeds <- c(123,124)
set.seed(seeds[1])

#### sample model(s) ####

# Number of states and parameters
num_states <- 10
num_extra_params <- 10
num_params <- num_states + num_extra_params # Beyond the base dependencies

# Define dependencies
# 1) Start with a baseline dependency: each state depends on one random state.
base_dependencies <- cbind(
  influenced = 1:num_states, 
  influencer = 1:num_states
)

# 2) Generate a pool of possible extra dependencies (all pairs)
all_pairs <- as.matrix(expand.grid(1:num_states, 1:num_states))
colnames(all_pairs) <- c("influenced", "influencer")
possible_dependencies <- all_pairs[apply(all_pairs, 1, function(x) x[1] != x[2]),]

# Now choose `num_extra_params` extra dependencies from the pool
extra_dependencies <- possible_dependencies[sample(1:nrow(possible_dependencies), num_extra_params),]

# Combine the base and extra dependencies
true_dependencies <- rbind(base_dependencies, extra_dependencies)

# Total number of dependencies
L <- nrow(true_dependencies)

# Generate true rate constants for all dependencies
k_true <- c(rnorm(num_states, -0.5, 0.25), rnorm(num_extra_params, 0, 0.5))

# Attach k to true_dependencies for clarity
true_dependencies <- cbind(true_dependencies, k = k_true)

#plot true dependencies graph
graph <- graph_from_data_frame(
  d = data.frame(
    from = true_dependencies[, "influencer"],
    to = true_dependencies[, "influenced"],
    weight = true_dependencies[, "k"]
  ),
  directed = TRUE
)
edge_weights <- E(graph)$weight
edge_thickness <- abs(edge_weights) / max(abs(edge_weights)) * 5
arrow_size <- abs(edge_weights) / max(abs(edge_weights)) * 0.5
colors <- adjustcolor(colorRampPalette(c("blue", "red"))(100), 0.9)
edge_colors <- colors[as.numeric(cut(edge_weights, breaks = 100))]
layout_wheel <- layout_in_circle(graph)

#### plot model(s) ####

# Adjust angles for self-loops to fan outward
self_loops <- 1:num_states  # Identify self-loops
loop_angles <- seq(0, -2 * pi, 
                   length.out = length(self_loops) + 1)  # Spread angles evenly
E(graph)$loop.angle <- rep(NA, ecount(graph))  # Initialize all angles to NA
E(graph)$loop.angle[self_loops] <- loop_angles[self_loops]  # Assign angles to self-loops

# plot the graph
plot.new()
plot.window(xlim = c(-1.5, 1.5), ylim = c(-1.5, 2), asp = 1)
par(mar = c(1,1,1,1))
plot(
  graph,
  layout = layout_wheel,
  edge.width = edge_thickness,
  edge.color = edge_colors,
  vertex.label = V(graph)$name,
  vertex.size = 20,
  vertex.color = "lightgrey",
  edge.arrow.size = arrow_size,
  add = T
)

# Draw the scale bar
library(sf)
source("~/repos/polylines/R/functions.R")
npts <- 50  # Number of points in the scale
k_range <- c(-1,1) * max(abs(range(true_dependencies[, "k"])))
k_seq <- seq(k_range[1], k_range[2], 
               length.out = npts)
scale_y <- seq(-0.5, 0.5, length.out = npts)  # Y positions for the gradient
scale_x <- rep(2, npts)  # Center x position for the hourglass
lwd <- abs(k_seq) / max(abs(k_seq)) * 0.1  # Scale the thickness based on |k|


# Use polylines to draw the hourglass
polylines(
  x = scale_x,  # Left and right edges
  y = scale_y,             # Mirror the y-coordinates
  lwd = lwd,                                # Horizontal thickness
  col = colorRampPalette(c("blue", "red"))(npts - 1),  # Gradient colors
  complex = FALSE,                          # Keep it simple
  border = NA                               # No border
)

# Add labels for the scale
text(x = scale_x, y = scale_y[1], round(k_range[1], 2), pos = 4)
text(x = scale_x, y = scale_y[npts], round(k_range[2], 2), pos = 4)
text(x = scale_x, y = scale_y[floor(npts/2)], 0, pos = 4)
text(x = scale_x, y = scale_y[npts] + diff(scale_y[1:2]) * 3, "dx/dy", pos = 3, cex = 1.1, font = 3)

#also a title
text(x = 0, y = 1.5, "Model #2", cex = 1.5, font = 2, col = 1, pos = 3, xpd = NA)

#### simulate data ####

# Initial conditions
y0 <- rnorm(num_states)

# ODE system function using the dependency matrix
ode_system <- function(t, y, parms) {
  # parms is a list with: 
  # parms$influenced, parms$influencer, parms$k
  dy <- numeric(num_states)
  # For each dependency: dy[influenced] += k * y[influencer]
  for (i in seq_along(parms$k)) {
    influenced <- parms$influenced[i]
    influencer <- parms$influencer[i]
    rate <- parms$k[i]
    dy[influenced] <- dy[influenced] + rate * y[influencer]
  }
  return(list(dy))
}

n_time <- 6
times <- c(0, seq(1, 10, length.out = n_time - 1))   # 10 timepoints
n_reps <- 3

# Simulate the system
params_list <- list(
  influenced = true_dependencies[,1],
  influencer = true_dependencies[,2],
  k = true_dependencies[,3]
)
out <- ode(y=y0, times=times, func=ode_system, parms=params_list)
true_states <- out[,-1]

# Determine a sigma for noise based on data spread
ts_sds <- apply(true_states, 2, sd)
true_sigma <- exp(mean(log(ts_sds)))

# Simulate noisy observations
obs_data <- array(NA, dim=c(n_time, n_reps, num_states))
for (i in seq_len(n_time)) {
  for (r in seq_len(n_reps)) {
    obs_data[i,r,] <- rnorm(num_states, mean=true_states[i,], sd=true_sigma)
  }
}

# Quick visualization (optional)
# quartz(width=5, height=5)
par(mfrow = c(4,3), mar = c(2,2,2,2))
for(k_i in 1:num_states){
  plot(true_states[,k_i], type = "l", 
       ylim = range(obs_data[,,k_i]), main = paste0("State ", k_i), lwd = 5)
  points(true_states[,k_i], cex = 2)
  for(i in 1:n_reps){
    lines(obs_data[,i,k_i], col = 2)
  }
}
plot.new()
plot.window(c(0,1), c(0,1))
legend("topleft", legend = c("True Values", "Observed Values"), 
       col = c(adjustcolor(1, 0.5), "red"), lty = c(1, 1), 
       pch = c(1, NA), lwd = c(5, 1), cex = 1.2, xpd = NA, bty = "n")

# Write Stan model
stan_model_code <- '
functions {
  vector ode_system(
    real t,
    vector y,
    vector theta,
    array[] real x_r,
    array[,] int dependencies
  ) {
    int M = size(y);
    array[2] int dep_dims = dims(dependencies);
    int L = dep_dims[1]; // number of rows
    // dep_dims[2] should be 2 (number of columns)
    vector[M] dy = rep_vector(0, M);

    for (l in 1:L) {
      int infl_state = dependencies[l, 1];
      int infl_by = dependencies[l, 2];
      dy[infl_state] += theta[l] * y[infl_by];
    }

    return dy;
  }
}
data {
  int<lower=1,upper=3> solver;
  int<lower=1> N;
  int<lower=1> R;
  int<lower=1> M;
  real t0;
  array[N] real ts;
  array[N,R,M] real y_obs;
  vector[M] y0;
  array[0] real x_r;
  int<lower=1> L;
  array[L,2] int dependencies; // You can declare them like this in data
}
parameters {
  vector[L] k;
  real k_self_mean;
  real<lower=0> k_self_sd;
  real<lower=0> sigma;
}
transformed parameters {
  array[N] vector[M] y_hat;
  if (solver == 1) {
    y_hat = ode_bdf(ode_system, y0, t0, ts, k, x_r, dependencies);
  }
  if (solver == 2) {
    y_hat = ode_rk45(ode_system, y0, t0, ts, k, x_r, dependencies);
  }
  if (solver == 3) {
    matrix[M, M] A = rep_matrix(0, M, M);
    for (l in 1:L) {
      int infl_state = dependencies[l, 1];
      int infl_by    = dependencies[l, 2];
      A[infl_state, infl_by] = k[l];
    }
    for (n in 1:N) {
      y_hat[n] = matrix_exp((ts[n] - t0) * A) * y0;
    }
  }
}
model {
  k_self_sd ~ normal(0,0.5);
  k_self_mean ~ normal(0,0.5);
  k[1:M] ~ normal(k_self_mean, k_self_sd);
  k[(M+1):L] ~ normal(0, 0.5);
  sigma ~ normal(0.1, 0.05);

  for (n in 1:N) {
    for (r in 1:R) {
      y_obs[n, r] ~ normal(y_hat[n], sigma);
    }
  }
}
'

#### cmdstanr fit ####
mod <- cmdstan_model(write_stan_file(stan_model_code), cpp_options = list(stan_threads = TRUE))
dependencies <- as.matrix(true_dependencies[, c("influenced", "influencer")])
dat <- list(
  solver = 3,                      # Use ode_rk45
  N = n_time,
  R = n_reps,
  M = num_states,
  t0 = -1E-6,                          # Initial time
  ts = times,
  y_obs = obs_data,
  y0 = y0,
  x_r = numeric(0),                # Empty real array
  dependencies = dependencies,     # 2D array of dependencies
  L = nrow(dependencies)           # Number of dependencies
)
fit <- mod$sample(
  refresh =10,
  data=dat,
  seed=123,
  chains=4,
  parallel_chains=4,
  iter_warmup=1000,
  iter_sampling=1000,
  threads_per_chain = 1,
  init = 0.05
)
fit_pf <- mod$pathfinder(
  refresh = 10,
  data=dat,
  seed=123,
  num_threads = 4,
  init = 0.05
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

#### rstan fit ####
rstan_model <- write_stan_file(stan_model_code)
fit_rstan <- stan(
  file = rstan_model,        # Path to your Stan model
  data = dat,              # Data list
  seed = 123,              # Seed for reproducibility
  chains = 4,              # Number of MCMC chains
  iter = 1000,             # Total iterations (warmup + sampling)
  warmup = 900,           # Warmup iterations
  init = 0.05,             # Initial values (scalar for all parameters)
  cores = 4                # Use 4 cores for parallel sampling
)
summ_rstan <- data.frame(summary(fit_rstan)$summary)
print(head(summ_rstan[order(summ_rstan$n_eff), c("Rhat", "n_eff",
                                   "mean", "sd")]))
print(head(summ_rstan[order(summ_rstan$Rhat, decreasing = T), c("Rhat", "n_eff",
                                               "mean", "sd")]))
fit.bridge <- bridge_sampler(fit_rstan, silent = TRUE)

#### false dependencies ####
set.seed(seeds[2])
false_extra_dependencies <- possible_dependencies[sample(1:nrow(possible_dependencies), num_extra_params),]
false_dependencies <- rbind(base_dependencies, false_extra_dependencies)
false_dat <- list(
  solver = 3,                      # Use ode_rk45
  N = n_time,
  R = n_reps,
  M = num_states,
  t0 = -1E-6,                          # Initial time
  ts = times,
  y_obs = obs_data,
  y0 = y0,
  x_r = numeric(0),                # Empty real array
  dependencies = false_dependencies,     # 2D array of dependencies
  L = nrow(false_dependencies)           # Number of dependencies
)
false_fit_rstan <- stan(
  file = rstan_model,        # Path to your Stan model
  data = false_dat,              # Data list
  chains = 4,              # Number of MCMC chains
  iter = 2000,             # Total iterations (warmup + sampling)
  warmup = 1000,           # Warmup iterations
  init = 0.05,             # Initial values (scalar for all parameters)
  cores = 4                # Use 4 cores for parallel sampling
)
false_summ_rstan <- data.frame(summary(false_fit_rstan)$summary)
print(head(false_summ_rstan[order(false_summ_rstan$n_eff), c("Rhat", "n_eff",
                                                 "mean", "sd")]))
print(head(false_summ_rstan[order(false_summ_rstan$Rhat, decreasing = T), c("Rhat", "n_eff",
                                                                "mean", "sd")]))
false_fit.bridge <- bridge_sampler(false_fit_rstan, silent = TRUE)

fit.bridge
false_fit.bridge

#### results for mcmc vs pathfinder ####
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
