#### packages ####
library(deSolve)
library(cmdstanr)
library(posterior)
library(rstan)
library(torch)
library(igraph)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### Stan Model ####
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
    int L = dep_dims[1];
    vector[M] dy = rep_vector(0, M);
    
    for (l in 1:L) {
      int infl_state = dependencies[l, 1];
      int infl_by    = dependencies[l, 2];
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
  array[L,2] int dependencies;
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
    // matrix_exp solver
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
  // Priors
  k_self_sd ~ normal(0,0.5);
  k_self_mean ~ normal(0,0.5);
  k[1:M] ~ normal(k_self_mean, k_self_sd);
  k[(M+1):L] ~ normal(0, 0.5);
  sigma ~ normal(0.1, 0.05);
  
  // Likelihood
  for (n in 1:N) {
    for (r in 1:R) {
      y_obs[n, r] ~ normal(y_hat[n], sigma);
    }
  }
}
'

#### Define ODE ####

#use random dependencies + self dependencies
num_states <- 10
num_extra_params <- 20
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

#evaluate dynamics by looking at adjacency matrix
adj <- matrix(0, num_states, num_states)
for (i in 1:nrow(true_dependencies)) {
  infl <- true_dependencies[i, "influenced"]
  by   <- true_dependencies[i, "influencer"]
  k    <- true_dependencies[i, "k"]
  adj[infl, by] <- k
}
real_components_evalues <- Re(eigen(adj)$values)

ode_system <- function(t, y, parms) {
  dy <- numeric(num_states)
  for (i in seq_along(parms$k)) {
    infl <- parms$influenced[i]
    by   <- parms$influencer[i]
    dy[infl] <- dy[infl] + parms$k[i]*y[by]
  }
  list(dy)
}

params_list <- list(
  influenced = true_dependencies[,1],
  influencer = true_dependencies[,2],
  k = true_dependencies[,3]
)

#### Plot ODE ####
library(igraph)
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
par(mar = c(1,1,1,1), mfrow = c(1,1))
plot.new()
plot.window(xlim = c(-1.5, 1.5), ylim = c(-1.5, 2), asp = 1)

# Adjust angles for self-loops to fan outward
self_loops <- 1:num_states
loop_angles <- seq(0, -2 * pi, length.out = length(self_loops) + 1)
E(graph)$loop.angle <- rep(NA, ecount(graph))
E(graph)$loop.angle[self_loops] <- loop_angles[self_loops]

plot(
  graph,
  layout = layout_wheel,
  edge.width = edge_thickness,
  edge.color = edge_colors,
  vertex.label = V(graph)$name,
  vertex.size = 20,
  vertex.color = "lightgrey",
  edge.arrow.size = arrow_size,
  add = TRUE
)

# Scale bar code (requires sf + polylines)
library(sf)
source("~/repos/polylines/R/functions.R")

npts <- 50
k_range <- c(-1,1) * max(abs(range(true_dependencies[, "k"])))
k_seq <- seq(k_range[1], k_range[2], length.out = npts)
scale_y <- seq(-0.5, 0.5, length.out = npts)
scale_x <- rep(1.75, npts)
lwd <- abs(k_seq) / max(abs(k_seq)) * 0.1

polylines(
  x = scale_x,
  y = scale_y,
  lwd = lwd,
  col = colorRampPalette(c("blue", "red"))(npts - 1),
  complex = FALSE,
  border = NA
)
text(x = scale_x, y = scale_y[1], round(k_range[1], 2), pos = 4, xpd = NA)
text(x = scale_x, y = scale_y[npts], round(k_range[2], 2), pos = 4, )
text(x = scale_x, y = scale_y[floor(npts/2)], 0, pos = 4)
text(x = scale_x, y = scale_y[npts] + diff(scale_y[1:2]) * 3, "dx/dy", pos = 3, cex = 1.1, font = 3)

# Add title
text(x = 0, y = 1.5, "Self + Random Dependencies", cex = 1.5, font = 2, col = 1,
     pos = 3, xpd = NA)

par(fig = c(0.0, 0.3, 0.0, 0.3), new = TRUE)
hist(real_components_evalues, main = "", xpd = NA, freq = F, cex.lab = 0.7, xlab = "Eigenvalues")
title(main = "Real Parts of\nEigenvalues of Jacobian", cex.main = 0.7, line = 0.5, xpd = NA)


#### Sample from ODE ####
library(deSolve)

# Choose how many timepoints & replicates
n_time <- 6
max_time <- 20
times <- c(0, seq(1, max_time, length.out = n_time - 1))
n_reps <- 3

# Solve the deterministic ODE for the "true" trajectory
y0 <- rnorm(num_states, 0, 1)
out <- ode(y = y0, times = times, func = ode_system, parms = params_list)
true_states <- out[, -1]  # drop the time column

# Add noisy observations for multiple replicates
ts_sds <- apply(true_states, 2, sd)
true_sigma <- exp(mean(log(ts_sds)))  # typical scale
obs_data <- array(NA, dim = c(n_time, n_reps, num_states))
for (i in seq_len(n_time)) {
  for (r in seq_len(n_reps)) {
    obs_data[i, r, ] <- rnorm(num_states, mean = true_states[i, ], sd = true_sigma)
  }
}

# Plot the "true vs. observed" timecourses
par(mfrow = c(4,3), mar = c(2,2,2,2))
for (k_i in 1:num_states) {
  plot(true_states[, k_i], type = "l", 
       ylim = range(obs_data[, , k_i]), 
       main = paste0("State ", k_i), lwd = 5)
  # points(true_states[, k_i], cex = 2)
  for (i in 1:n_reps) {
    lines(obs_data[, i, k_i], col = 2)
  }
  lines(true_states[, k_i], lwd = 5)
}
plot.new()
plot.window(c(0,1), c(0,1))
legend("topleft", legend = c("True Values", "Observed Values"), 
       col = c(adjustcolor(1, 0.5), "red"), lty = c(1, 1), 
       pch = c(1, NA), lwd = c(5, 1), cex = 1.2, xpd = NA, bty = "n")

#### Generate x-sectional data ####
# 2d) We'll generate many time-series with random initial conditions
N_series <- 500  # number of time series
Tmax <- 20        # maximum time
times_fine <- seq(0, Tmax, length.out=20)

# For cross-sectional data, we pick exactly 1 random timepoint from each series
X_cross <- matrix(NA, nrow=N_series, ncol=num_states)

for (i in seq_len(N_series)) {
  # random initial conditions
  y0_i <- rnorm(num_states)
  sol <- ode(y=y0_i, times=times_fine,
                      func=ode_system, parms=params_list)
  sol_states <- sol[,-1]  # drop times column -> 50 x num_states
  
  # pick one random timepoint
  idx <- sample(nrow(sol_states), 1)
  X_cross[i,] <- sol_states[idx,]
}


#### fit VAE ####

# R + torch VAE
VAE <- nn_module(
  "VAE",
  initialize = function(input_dim, hidden_dims = c(32, 16), latent_dim = 2) {
    self$input_dim <- input_dim
    self$latent_dim <- latent_dim
    
    # Encoder
    self$enc_fc1 <- nn_linear(input_dim, hidden_dims[1])
    self$enc_fc2 <- nn_linear(hidden_dims[1], hidden_dims[2])
    self$enc_mu <- nn_linear(hidden_dims[2], latent_dim)
    self$enc_logvar <- nn_linear(hidden_dims[2], latent_dim)
    
    # Decoder
    self$dec_fc1 <- nn_linear(latent_dim, hidden_dims[2])
    self$dec_fc2 <- nn_linear(hidden_dims[2], hidden_dims[1])
    self$dec_out <- nn_linear(hidden_dims[1], input_dim)
  },
  encode = function(x) {
    h1 <- torch_relu(self$enc_fc1(x))
    h2 <- torch_relu(self$enc_fc2(h1))
    mu <- self$enc_mu(h2)
    logvar <- self$enc_logvar(h2)
    list(mu, logvar)
  },
  reparameterize = function(mu, logvar) {
    std <- torch_exp(0.5 * logvar)
    eps <- torch_randn_like(std)
    mu + eps * std
  },
  decode = function(z) {
    h1 <- torch_relu(self$dec_fc1(z))
    h2 <- torch_relu(self$dec_fc2(h1))
    # For purely continuous data, no final activation or maybe just identity
    self$dec_out(h2)
  },
  forward = function(x) {
    encoded <- self$encode(x)
    z <- self$reparameterize(encoded[[1]], encoded[[2]])
    decoded <- self$decode(z)
    list(decoded, encoded[[1]], encoded[[2]])
  }
)

vae_loss <- function(recon_x, x, mu, logvar) {
  recon_loss <- nnf_mse_loss(recon_x, x, reduction = "sum")
  kl_loss <- -0.5 * torch_sum(1 + logvar - mu^2 - torch_exp(logvar))
  recon_loss + kl_loss
}

# Convert cross-sectional data to torch tensor
X_tensor <- torch_tensor(X_cross, dtype=torch_float())

# Create dataloader
batch_size <- 128
dataset <- tensor_dataset(X_tensor)
dloader <- dataloader(dataset, batch_size=batch_size, shuffle=TRUE)

# Instantiate VAE
vae <- VAE(
  input_dim = num_states,
  hidden_dims = c(32, 16, 8),
  latent_dim = 8
)
optimizer <- optim_adam(vae$parameters, lr=1e-3)

# Train VAE
epochs <- 2000
for (epoch in seq_len(epochs)) {
  total_loss <- 0
  coro::loop(for (b in dloader) {
    optimizer$zero_grad()
    out <- vae(b[[1]])
    loss <- vae_loss(out[[1]], b[[1]], out[[2]], out[[3]])
    loss$backward()
    optimizer$step()
    total_loss <- total_loss + loss$item()
  })
  if (epoch %% 50 == 0) {
    cat(sprintf("Epoch %d | Loss: %3.2f\n", epoch, total_loss))
  }
}

#### fit MVN ####
use_asinh <- F
return_mean <- F
nit_cmvn <- 100
if(use_asinh){
  sigma <- cov(asinh(X_cross))
  mu <- colMeans(asinh(X_cross))
  out_cmvn <- sinh(iterate_cmvn_traj(mu, sigma, asinh(y0), nit_cmvn, return_mean))
} else {
  sigma <- cov(X_cross)
  mu <- colMeans(X_cross)
  out_cmvn <- iterate_cmvn_traj(mu, sigma, y0, nit_cmvn, return_mean)
}


#### Iterate fitted VAE ####

start_idx <- 1
y0 <- rep(0, num_states)
y0 <- rnorm(num_states)
start_state <- torch_tensor(y0, dtype = torch_float())
n_iter <- 5000
fraction_to_next <- 1/50
out_vae <- array(NA, dim = c(n_iter+1, num_states))
out_vae[1,] <- as.numeric(start_state)

current <- start_state$clone()
curr <- as.numeric(current)
for (tstep in seq_len(n_iter)) {
  # Insert a dimension so shape is (1, 10)
  out <- vae(current$unsqueeze(dim = 1))
  # out[[1]] is shape (1, 10) (the reconstructed output)
  next_state <- out[[1]][1, ] 
  next_state <- current + (next_state - current) * fraction_to_next
  
  nxt <- as.numeric(next_state)
  out_vae[tstep+1,] <- nxt
  
  current <- next_state$clone()
  curr <- as.numeric(current)
}

#Compare to ODE simulation
nt_ode <- 200
max_t <- 6
times_ode <- seq(0, max_t, length.out = nt_ode+1)  
out_ode <- ode(y = y0, times = times_ode, func = ode_system, 
               parms = params_list)
out_ode <- out_ode[,-1]

#plot
par(mfrow = c(4, 3), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 1))  # ← added oma for outer margin
for(i in 1:9){
  plot(0:nt_ode/nt_ode, out_ode[,i], type = "l", col = "darkblue", lwd = 2,
       ylim = range(c(out_ode[,i], out_vae[,i], mu[i])), xlab = "", ylab = "")
  lines(0:n_iter/n_iter, out_vae[,i], col = "darkred", lwd = 2)
  abline(h = mu[i], lty = 3, lwd = 2, col = "darkviolet")
}
mtext("Relative Time", side = 1, outer = TRUE, line = -9, cex = 1.5, col = 1)
mtext("\t\tState Value", side = 2, outer = TRUE, line = 1, 
      cex = 1.5, col = 1)

#add legend
plot.new()
plot.window(c(0,1), c(0,1))
legend("topleft", 
       legend = c("ODE Simulation",  
                  "Iterative VAE Output",
                  "X-Sectional Sample Mean"), 
       col = c("darkblue", "darkred", "darkviolet"), lty = c(1,1,3), 
       lwd = c(2, 2, 2), cex = 1.2, xpd = NA, bty = "n")


#### iter cond-MVN? ####
conditional_mvn <- function(mu, sigma, observed_indices, observed_values) {
  # Check input dimensions
  p <- length(mu)
  if (!all(dim(sigma) == c(p, p))) {
    stop("sigma matrix must be p x p.")
  }
  if (length(observed_indices) != length(observed_values)) {
    stop("Length of observed_indices must match length of observed values.")
  }
  
  # Identify observed_indices for unobserved and observed variables
  all_indices <- 1:p
  unobserved_indices <- setdiff(all_indices, observed_indices)
  
  # Partition mu vector and sigma matrix
  mu1 <- mu[unobserved_indices]   # mu of unobserved variables
  mu2 <- mu[observed_indices]              # mu of observed variables
  
  Sigma11 <- sigma[unobserved_indices, unobserved_indices] # sigma of unobserved
  Sigma12 <- sigma[unobserved_indices, observed_indices]            # Cross-sigma
  Sigma21 <- sigma[observed_indices, unobserved_indices]
  Sigma22 <- sigma[observed_indices, observed_indices]                       # sigma of observed
  Sigma22_inv <- solve(Sigma22) # Inverse of Sigma22
  
  # Compute conditional mu and sigma
  conditional_mu <- mu1 + Sigma12 %*% Sigma22_inv %*% (observed_values - mu2)
  conditional_sigma <- Sigma11 - Sigma12 %*% Sigma22_inv %*% Sigma21
  
  # Return result as a list
  list(mu = as.vector(conditional_mu), sigma = conditional_sigma)
}

conditional_mvn_mean_loo <- function(mu, sigma, observed_values, return_mean = T){
  p <- length(mu)
  cond_dists <- lapply(1:p, function(idx){
    out <- conditional_mvn(mu, sigma, setdiff(1:p, idx), 
                    observed_values[setdiff(1:p, idx)])
    return(c(out$mu, out$sigma))
  })
  mus <- do.call(rbind, cond_dists)[,1]
  if(return_mean){
    return(mus)
  } else {
    return(observed_values + (mus - observed_values) / 10)  
  }
}

iterate_cmvn_traj <- function(mu, sigma, y0, nsteps, return_mean = T){
  p <- length(mu)
  trajmat <- matrix(0, nsteps+1, p)
  trajmat[1,] <- y0
  for(i in 2:(nsteps + 1)){
    trajmat[i,] <- conditional_mvn_mean_loo(mu, sigma, 
                                            trajmat[i-1,],
                                            return_mean)
  }
  return(trajmat)
}

#simulate again
use_asinh <- F
return_mean <- F
nit_cmvn <- 1000
if(use_asinh){
  sigma <- cov(asinh(X_cross))
  mu <- colMeans(asinh(X_cross))
  out_cmvn <- sinh(iterate_cmvn_traj(mu, sigma, asinh(y0), nit_cmvn, return_mean))
} else {
  sigma <- cov(X_cross)
  mu <- colMeans(X_cross)
  out_cmvn <- iterate_cmvn_traj(mu, sigma, y0, nit_cmvn, return_mean)
}

#plot
par(mfrow = c(4, 3), mar = c(2, 2, 2, 2), oma = c(4, 4, 2, 1))  # ← added oma for outer margin
for(i in 1:9){
  plot(0:nt_ode/nt_ode, out_ode[,i], type = "l", col = "darkblue", lwd = 2,
       ylim = range(c(out_ode[,i], out_vae[,i], out_cmvn[,i], mu[i])), xlab = "", ylab = "")
  lines(0:n_iter/n_iter, out_vae[,i], col = adjustcolor("darkred", 0.25), lwd = 2)
  lines(0:nit_cmvn/nit_cmvn, out_cmvn[,i], col = "darkorange", lwd = 2)
  abline(h = mu[i], lty = 3, lwd = 2, col = "darkviolet")
}
mtext("Relative Time", side = 1, outer = TRUE, line = -9, cex = 1.5, col = 1)
mtext("\t\tState Value", side = 2, outer = TRUE, line = 1, 
      cex = 1.5, col = 1)

#add legend
plot.new()
plot.window(c(0,1), c(0,1))
legend("topleft", 
       legend = c("ODE Simulation", 
                  "Iterative LOO Conditional MVN", 
                  "Iterative VAE Output",
                  "X-Sectional Sample Mean"), 
       col = c("darkblue", "darkorange", adjustcolor("darkred", 0.25), "darkviolet"), lty = c(1,1,1,3), 
       lwd = c(2, 2), cex = 1.2, xpd = NA, bty = "n")


#### fit ODE model to VAE output ####

# same dependencies as we used to generate the data
dependencies <- true_dependencies[,c("influenced","influencer")]

# We'll call the times 0..n_iter
tsteps <- 0:n_iter
N_tsteps <- length(tsteps)

# Build an array [time x replicate x state]
# replicate=1 here
obs_data_iter <- array(NA, dim=c(N_tsteps, 1, num_states))
obs_data_iter[,1,] <- out_vae

# Construct data list for cmdstanr using your original Stan code
dat_iter <- list(
  solver = 3,  # use matrix_exp solver
  N = N_tsteps,
  R = 1,
  M = num_states,
  t0 = -1e-6,
  ts = tsteps,
  y_obs = obs_data_iter,
  y0 = out_vae[1,],
  x_r = numeric(0),
  dependencies = dependencies,
  L = nrow(dependencies)
)

# Compile Stan model with cmdstanr
mod <- cmdstan_model(write_stan_file(stan_model_code), cpp_options = list(stan_threads = TRUE))

# Fit
fit_iter <- mod$sample(
  data = dat_iter,
  seed = 999,
  chains = 4,
  parallel_chains = 4,
  iter_warmup=500,
  iter_sampling=500,
  refresh=100
)

# Check summary
summ_iter <- fit_iter$summary()
print(summ_iter[order(summ_iter$ess_bulk), 
                c("variable","mean","sd","rhat","ess_bulk","ess_tail")])

# Or pathfinder, if you like
fit_iter_pf <- mod$pathfinder(
  data=dat_iter,
  seed=999,
  num_threads=4,
  init=0.05
)
fit_iter_pf$draws() |> as_draws_df() |> summary()

#### Inspect Fit ####
# Plot the VAE-iterated path in one or two states vs. time, for instance:

par(mfrow=c(2,2))
plot(0:n_iter, out_vae[,1], type="o", main="Iterated State 1")
plot(0:n_iter, out_vae[,2], type="o", main="Iterated State 2")
