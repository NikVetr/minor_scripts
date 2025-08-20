# --- 0. Load Libraries ---
library(cmdstanr)
library(posterior)
library(dplyr)

# --- 1. Simulation Parameters ---
K <- 20         # Number of groups (distributions)
Ns <- c(2,2)
N_vec <- Ns[rbinom(K, 1, 0.5) + 1] # Vector of sample sizes per group
N_total <- sum(N_vec) # Total number of observations

# True parameter values for simulation
true_tau <- 1
true_mu <- rnorm(K, 0, true_tau) # True group means
sigma_y <- 10                   # Fixed observation SD

cat("Simulating data with:\n")
cat("K =", K, "\n")
cat("N per group:", N_vec, "\n")
cat("Total N =", N_total, "\n")
cat("True tau =", round(true_tau, 3), "\n")
# cat("True mu:", round(true_mu, 3), "\n") # Optional: print true means

# --- 2. Simulate Data ---
y_list <- vector("list", K)
group_idx_list <- vector("list", K)

for (k in 1:K) {
  y_list[[k]] <- rnorm(N_vec[k], true_mu[k], sigma_y)
  group_idx_list[[k]] <- rep(k, N_vec[k])
}

# Combine into single vectors for Stan
y_obs <- unlist(y_list)
group_idx <- unlist(group_idx_list)

# Check lengths
stopifnot(length(y_obs) == N_total)
stopifnot(length(group_idx) == N_total)

# --- 3. Prepare Data for Stan ---
stan_data <- list(
  K = K,
  N = N_vec,
  N_total = N_total,
  y = y_obs,
  group_idx = group_idx,
  sigma_y = sigma_y # Pass fixed sigma_y to Stan
)

# --- 4. Define Stan Models ---

# Centered Parameterization (CP)
stan_code_cp <- "
data {
  int<lower=1> K;         // Number of groups
  array[K] int<lower=1> N; // Sample sizes per group
  int<lower=1> N_total;   // Total number of observations
  vector[N_total] y;      // Observations
  array[N_total] int<lower=1, upper=K> group_idx; // Group index for each obs
  real<lower=0> sigma_y; // Fixed observation SD
}
parameters {
  real<lower=0> tau;      // Group mean standard deviation (hyperparameter)
  vector[K] mu;           // Group means (parameters)
}
model {
  // Priors
  tau ~ normal(0, 1);     // Half-normal prior for tau (implicit lower=0 bound)
  mu ~ normal(0, tau);    // Prior for group means (CENTERED)

  // Likelihood (vectorized)
  y ~ normal(mu[group_idx], sigma_y);
}
"

# Non-Centered Parameterization (NCP)
stan_code_ncp <- "
data {
  int<lower=1> K;         // Number of groups
  array[K] int<lower=1> N; // Sample sizes per group
  int<lower=1> N_total;   // Total number of observations
  vector[N_total] y;      // Observations
  array[N_total] int<lower=1, upper=K> group_idx; // Group index for each obs
  real<lower=0> sigma_y; // Fixed observation SD
}
parameters {
  real<lower=0> tau;      // Group mean standard deviation (hyperparameter)
  vector[K] mu_raw;       // Raw group means (sampled from standard normal)
}
transformed parameters {
  vector[K] mu = tau * mu_raw; // Actual group means (NON-CENTERED calculation)
                            // Assumes global mean = 0
}
model {
  // Priors
  tau ~ normal(0, 1);       // Half-normal prior for tau
  mu_raw ~ std_normal();    // Prior for raw parameters (independent of tau)

  // Likelihood (vectorized)
  // Uses the transformed 'mu'
  y ~ normal(mu[group_idx], sigma_y);
}
"

stan_code_sqrt <- "
data {
  int<lower=1> K;         // Number of groups
  array[K] int<lower=1> N; // Sample sizes per group
  int<lower=1> N_total;   // Total number of observations
  vector[N_total] y;      // Observations
  array[N_total] int<lower=1, upper=K> group_idx; // Group index for each obs
  real<lower=0> sigma_y; // Fixed observation SD
}
parameters {
  real<lower=0> tau;      // Group mean standard deviation (hyperparameter)
  vector[K] z;           // Transformed parameters
}
transformed parameters {

  vector[K] z_sign;  
  for (k in 1:K) {
    z_sign[k] = z[k] < 0 ? -1 : 1;
  }
  
  vector[K] mu;           // Actual group means
  for (k in 1:K) {
    // Reconstruct mu from z: mu = sign(z) * z^2
    // Using pow() might be slightly more stable than z[k]*abs(z[k])
    mu[k] = z_sign[k] * pow(z[k], 2);
  }
}
model {
  // Priors
  tau ~ normal(0, 1);     // Half-normal prior for tau

  // Prior on mu (applied implicitly via z) + Jacobian adjustment
  // log p(z) = log p(mu(z)) + log | jacobian |
  // Target increment is log p(mu(z)) + log | jacobian |

  // Prior term: log p(mu(z))
  target += normal_lpdf(mu | 0, tau); // Prior is conceptually on mu

  // Jacobian adjustment term: log | d(mu)/d(z) |
  // Need to add log(2*|z[k]|) for each k
  // Use vectorization where possible
  // Note: log(fabs(z)) can be problematic if z hits exactly 0,
  // though HMC should avoid this. Stan handles log(0) safely (-inf).
  target += K * log(2.0); // log(2) for each of the K components
  target += sum(log(abs(z))); // Sum of log(|z[k]|) for k=1..K

  // Likelihood (vectorized)
  y ~ normal(mu[group_idx], sigma_y);
}"


#maybe try squaring
stan_code_sqrt <- "// Stan Model: Inverse Sqrt-like Transformation (INV-SQRT-P)
data {
  int<lower=1> K;         // Number of groups
  array[K] int<lower=1> N; // Sample sizes per group
  int<lower=1> N_total;   // Total number of observations
  vector[N_total] y;      // Observations
  array[N_total] int<lower=1, upper=K> group_idx; // Group index for each obs
  real<lower=0> sigma_y; // Fixed observation SD
}
parameters {
  real<lower=0> tau;           // Group mean standard deviation (hyperparameter)
  vector<lower=0>[K] mu_sq;  // Sample the squared magnitude
  vector[K] s_raw;           // Sample auxiliary variable for sign
}
transformed parameters {
  vector[K] mu;                // Actual group means
  for (k in 1:K) {
    // Reconstruct mu = sign(s_raw) * sqrt(mu_sq)
    // sign() function is 1 if >0, -1 if <0, 0 if 0.
    // Multiplying by 0 if s_raw is 0 is fine as mu_sq would also be 0.
    // Add small epsilon to sqrt for numerical stability near 0 if needed,
    // but typically okay.
    real sign_k = s_raw[k] > 0 ? 1.0 : (s_raw[k] < 0 ? -1.0 : 0.0);
    mu[k] = sign_k * sqrt(mu_sq[k]);
  }
}
model {
  // Prior for tau
  tau ~ normal(0, 1);     // Half-normal prior for tau

  // Prior for sign (implicitly handles the Bernoulli(0.5) for sign(mu))
  s_raw ~ std_normal();   // Equivalent to sign being +/- 1 with prob 0.5

  // Prior for mu_sq derived from mu ~ normal(0, tau)
  // mu_sq ~ Gamma(0.5, 1.0 / (2.0 * tau^2));
  // Need to handle tau=0 case? The rate would go to infinity.
  // Stan's gamma_lpdf should handle rate=infinity gracefully (density is 0 if mu_sq>0).
  // Avoid division by zero explicitly if worried:
  if (tau > 1e-9) { // Add a small tolerance check
      target += gamma_lpdf(mu_sq | 0.5, 1.0 / (2.0 * tau^2));
  } else {
      // If tau is effectively zero, mu must be zero, so mu_sq must be zero.
      // gamma_lpdf handles mu_sq=0 correctly (log-density is 0 if shape > 0)
      // but the rate -> infinity implies density is 0 for mu_sq > 0.
      // Check if mu_sq is near zero. A simple approach is to let Stan handle it,
      // or add a penalty if mu_sq is not zero when tau is zero.
      // For simplicity here, let Stan's gamma_lpdf handle it.
      // Alternatively, implement check: if any(mu_sq > 1e-9) target += -infinity;
      target += gamma_lpdf(mu_sq | 0.5, 1.0 / (2.0 * tau^2 + 1e-9)); // Add epsilon to denom
  }


  // Likelihood (vectorized)
  y ~ normal(mu[group_idx], sigma_y);
}"
# --- 5. Compile Stan Models ---
# Create temporary files for Stan code
file_cp <- write_stan_file(stan_code_cp)
file_ncp <- write_stan_file(stan_code_ncp)
file_sqrt <- write_stan_file(stan_code_sqrt)

# Compile using cmdstanr
# Consider adding `cpp_options = list(stan_threads = TRUE)` and `threads_per_chain = N`
# if you have multiple cores available for parallel execution within chains.
model_cp <- cmdstan_model(file_cp, quiet = FALSE)
model_ncp <- cmdstan_model(file_ncp, quiet = FALSE)
model_sqrt <- cmdstan_model(file_sqrt, quiet = FALSE)

# --- 6. Run Stan Models ---
# MCMC settings
n_iter <- 1000
n_warmup <- 500
n_chains <- 4

fit_cp <- model_cp$sample(
  data = stan_data,
  chains = n_chains,
  parallel_chains = min(n_chains, parallel::detectCores()), # Use available cores
  iter_warmup = n_warmup,
  iter_sampling = n_iter - n_warmup,
  refresh = 100, # How often to print progress
  # adapt_delta = 0.8 # Default is 0.8, might need increasing for CP if divergences occur
  max_treedepth = 10 # Default is 10
)

fit_ncp <- model_ncp$sample(
  data = stan_data,
  chains = n_chains,
  parallel_chains = min(n_chains, parallel::detectCores()),
  iter_warmup = n_warmup,
  iter_sampling = n_iter - n_warmup,
  refresh = 100
  # No need to typically adjust adapt_delta or max_treedepth for NCP
)

fit_sqrt <- model_sqrt$sample(
  data = stan_data,
  chains = n_chains,
  parallel_chains = min(n_chains, parallel::detectCores()),
  iter_warmup = n_warmup,
  iter_sampling = n_iter - n_warmup,
  refresh = 100
)

summary_cp <- fit_cp$summary("mu")
summary_ncp <- fit_ncp$summary("mu")
summary_sqrt <- fit_sqrt$summary("mu")

# --- 8. Compare ESS Specifically ---
n_head <- 3
print(head(summary_cp[order(summary_cp$ess_bulk), c("variable", 
                                               "rhat", "ess_bulk", "ess_tail", 
                                               "mean", "sd")], n_head))
print(head(summary_ncp[order(summary_cp$ess_bulk), c("variable", 
                                                     "rhat", "ess_bulk", "ess_tail", 
                                                     "mean", "sd")], n_head))
print(head(summary_sqrt[order(summary_cp$ess_bulk), c("variable", 
                                                     "rhat", "ess_bulk", "ess_tail", 
                                                     "mean", "sd")], n_head))
