# Load necessary libraries
library(cmdstanr)
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(posterior)  # For summarizing posterior draws

# --- Stan Model Code as String ---
# Using R's raw string literal r"(...)" to avoid escaping characters
stan_model_code <- r"(
data {
  int<lower=0> n;
  int<lower=0> n_gene;
  int<lower=0> n_indiv;
  array[n] int<lower=1, upper=n_gene> gene_i;
  array[n] int<lower=1, upper=n_indiv> indiv_j;
  array[n] int<lower=1, upper=2> eQTL_het; //2 if they have a cis-eqtl, 1 otherwise
  vector[n] ase_coef;
}

parameters {
  
  //logit mixture param for allelic balance
  real mu_base;
  vector[n_indiv] mu_indiv;
  real<lower=0> sd_mu_indiv;
  vector[n_gene] mu_gene;
  real<lower=0> sd_mu_gene;
  real mu_eQTL_het;
  
  //sigma parameter
  real sigma_base;
  vector[n_indiv] sigma_indiv;
  real<lower=0> sd_sigma_indiv;
  vector[n_gene] sigma_gene;
  real<lower=0> sd_sigma_gene;
  real sigma_eQTL_het;
  
}

transformed parameters {
  
}

model {
  
  //priors and hyperpriors
  // Priors for variance of ASE coef
  sigma_base ~ std_normal();
  sigma_indiv ~ std_normal();
  sd_sigma_indiv ~ std_normal();
  sigma_gene ~ std_normal();
  sd_sigma_gene ~ std_normal();
  sigma_eQTL_het ~ std_normal();
  
  // Priors for mean of ASE coef
  mu_base ~ std_normal();
  mu_indiv ~ std_normal();
  sd_mu_indiv ~ std_normal();
  mu_gene ~ std_normal();
  sd_mu_gene ~ std_normal();
  mu_eQTL_het ~ std_normal();
  
  //model
  
  //////////////////////////
  // sigma scale PARAMS //
  //////////////////////////
  
  //construct mean eQTL_het effect
  vector[2] sigma_eQTL_het_split;
  sigma_eQTL_het_split[1] = sigma_eQTL_het / 2;
  sigma_eQTL_het_split[2] = -sigma_eQTL_het / 2;
  
  //construct positive sigma scale vector
  vector[n] sigma = exp(sigma_base + 0 +
    sigma_gene[gene_i] * sd_sigma_gene + 
    sigma_indiv[indiv_j] * sd_sigma_indiv +
    sigma_eQTL_het_split[eQTL_het]
  );

  // mean of ase coefs
  
  //construct mean eQTL_het effect
  vector[2] mu_eQTL_het_split;
  mu_eQTL_het_split[1] = mu_eQTL_het / 2;
  mu_eQTL_het_split[2] = -mu_eQTL_het / 2;
  
  //construct positive mu scale vector
  vector[n] mu = (mu_base + 0 +
    mu_gene[gene_i] * sd_mu_gene + 
    mu_indiv[indiv_j] * sd_mu_indiv +
    mu_eQTL_het_split[eQTL_het]
  );

  ////////////////
  // Likelihood //
  ////////////////
  
  ase_coef ~ normal(mu, sigma);
  
}
)" # End of raw string literal


# --- Simulation Setup ---

# Set seed for reproducibility
set.seed(123)

# Define simulation dimensions
n_gene_sim <- 100
n_indiv_sim <- 5
# Assuming full crossing for simplicity (each indiv has data for each gene)
n_sim <- n_gene_sim * n_indiv_sim

# Generate index variables
gene_i_sim <- rep(1:n_gene_sim, times = n_indiv_sim)
indiv_j_sim <- rep(1:n_indiv_sim, each = n_gene_sim)

# Generate eQTL status (e.g., ~30% have eQTL status 2)
eQTL_het_sim <- sample(1:2, size = n_sim, replace = TRUE, prob = c(0.7, 0.3))

# --- Define TRUE Parameter Values ---
# Choose values that are somewhat reasonable under std_normal priors,
# but specific enough to test recovery.

true_params <- list()

# Mu parameters
true_params$mu_base <- 0.2
true_params$sd_mu_indiv <- 0.5 # True SD for individual effects on mu
true_params$sd_mu_gene <- 1   # True SD for gene effects on mu
true_params$mu_eQTL_het <- 0.6  # True difference effect for eQTL on mu

# Draw standardized random effects from N(0,1) as in the model's priors
true_params$mu_indiv <- rnorm(n_indiv_sim, 0, 1)
true_params$mu_gene <- rnorm(n_gene_sim, 0, 1)

# Sigma parameters (remember these are on log scale before exp())
true_params$sigma_base <- log(0.8) # Base log-sd -> sd = 0.8
true_params$sd_sigma_indiv <- 0.3 # True SD for individual effects on log_sigma
true_params$sd_sigma_gene <- 0.4  # True SD for gene effects on log_sigma
true_params$sigma_eQTL_het <- -0.5 # True difference effect for eQTL on log_sigma

# Draw standardized random effects from N(0,1)
true_params$sigma_indiv <- rnorm(n_indiv_sim, 0, 1)
true_params$sigma_gene <- rnorm(n_gene_sim, 0, 1)


# --- Generate Data Based on TRUE Parameters ---

# Calculate mu and sigma for each observation using the TRUE parameters
# Mimic the structure in the Stan model's transformed parameters block

# Split eQTL effects
mu_eQTL_het_split_true <- c(true_params$mu_eQTL_het / 2, -true_params$mu_eQTL_het / 2)
sigma_eQTL_het_split_true <- c(true_params$sigma_eQTL_het / 2, -true_params$sigma_eQTL_het / 2)

# Calculate linear predictors (log scale for sigma, also log scale for mu based on model)
log_mu_linpred_true <- numeric(n_sim)
log_sigma_linpred_true <- numeric(n_sim)

for (i in 1:n_sim) {
  log_mu_linpred_true[i] <- true_params$mu_base +
    true_params$mu_gene[gene_i_sim[i]] * true_params$sd_mu_gene +
    true_params$mu_indiv[indiv_j_sim[i]] * true_params$sd_mu_indiv +
    mu_eQTL_het_split_true[eQTL_het_sim[i]]
  
  log_sigma_linpred_true[i] <- true_params$sigma_base +
    true_params$sigma_gene[gene_i_sim[i]] * true_params$sd_sigma_gene +
    true_params$sigma_indiv[indiv_j_sim[i]] * true_params$sd_sigma_indiv +
    sigma_eQTL_het_split_true[eQTL_het_sim[i]]
}

# Calculate final mu and sigma on the outcome scale
mu_true <- (log_mu_linpred_true)
sigma_true <- exp(log_sigma_linpred_true)

# Generate the response variable 'ase_coef'
ase_coef_sim <- rnorm(n = n_sim, mean = mu_true, sd = sigma_true)

# --- Prepare Data List for Stan ---
stan_data <- list(
  n = n_sim,
  n_gene = n_gene_sim,
  n_indiv = n_indiv_sim,
  gene_i = gene_i_sim,
  indiv_j = indiv_j_sim,
  eQTL_het = eQTL_het_sim,
  ase_coef = ase_coef_sim
)

# --- Compile and Fit the Stan Model ---

# Compile the model from the string variable 'stan_model_code'
# This might take a minute the first time
mod <- cmdstan_model(stan_file = write_stan_file(stan_model_code))

# Fit the model using MCMC
# Use multiple chains and cores for better diagnostics
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 10,      # Print progress every 500 iterations
  adapt_delta = 0.85  # Adjust if divergences occur
)

# --- Check MCMC Diagnostics ---
# Check for divergences, Rhat, effective sample size
print(fit$summary())
# Check specific diagnostics
fit$cmdstan_diagnose()

# --- Evaluate Parameter Recovery ---

# Extract posterior summaries
posterior_summary <- fit$summary()

# Function to compare a single scalar parameter
compare_scalar <- function(param_name, true_value, summary_df) {
  est_row <- summary_df %>% filter(variable == param_name)
  if (nrow(est_row) == 0) {
    warning("Parameter ", param_name, " not found in summary.")
    return(NULL)
  }
  data.frame(
    parameter = param_name,
    true = true_value,
    mean = est_row$mean,
    median = est_row$median,
    sd = est_row$sd,
    q5 = est_row$q5,
    q95 = est_row$q95,
    rhat = est_row$rhat,
    ess_bulk = est_row$ess_bulk,
    ess_tail = est_row$ess_tail
  )
}

# Function to compare vector parameters (like random effects)
compare_vector <- function(param_name_base, true_vector, summary_df) {
  # Match names like "mu_gene[1]", "mu_gene[2]", ...
  pattern <- paste0("^", param_name_base, "\\[")
  est_rows <- summary_df %>% filter(grepl(pattern, variable))
  if (nrow(est_rows) == 0) {
    warning("Parameter vector ", param_name_base, " not found in summary.")
    return(NULL)
  }
  # Ensure order matches
  indices <- as.numeric(sub(paste0(param_name_base, "\\[([0-9]+)\\]"), "\\1", est_rows$variable))
  est_rows <- est_rows[order(indices), ]
  
  if(length(true_vector) != nrow(est_rows)) {
    warning("Length mismatch for ", param_name_base, ": True=", length(true_vector), ", Est=", nrow(est_rows))
    return(NULL)
  }
  
  data.frame(
    parameter = est_rows$variable,
    index = 1:length(true_vector),
    true = true_vector,
    mean = est_rows$mean,
    median = est_rows$median,
    sd = est_rows$sd,
    q5 = est_rows$q5,
    q95 = est_rows$q95,
    rhat = est_rows$rhat,
    ess_bulk = est_rows$ess_bulk,
    ess_tail = est_rows$ess_tail
  )
}

# Compare scalar parameters
scalar_params <- c("mu_base", "sd_mu_indiv", "sd_mu_gene", "mu_eQTL_het",
                   "sigma_base", "sd_sigma_indiv", "sd_sigma_gene", "sigma_eQTL_het")

scalar_comparison <- bind_rows(lapply(scalar_params, function(p) {
  compare_scalar(p, true_params[[p]], posterior_summary)
}))

# Compare vector parameters
vector_comparison <- list()
vector_comparison$mu_indiv <- compare_vector("mu_indiv", true_params$mu_indiv, posterior_summary)
vector_comparison$mu_gene <- compare_vector("mu_gene", true_params$mu_gene, posterior_summary)
vector_comparison$sigma_indiv <- compare_vector("sigma_indiv", true_params$sigma_indiv, posterior_summary)
vector_comparison$sigma_gene <- compare_vector("sigma_gene", true_params$sigma_gene, posterior_summary)

# Print comparisons
cat("\n--- Scalar Parameter Recovery ---\n")
print(scalar_comparison %>% select(parameter, true, mean, q5, q95, rhat, ess_bulk))

cat("\n--- Vector Parameter Recovery (Example: mu_gene) ---\n")
# Quick look at one vector parameter
if (!is.null(vector_comparison$mu_gene)) {
  print(head(vector_comparison$mu_gene %>% select(parameter, true, mean, q5, q95, rhat)))
} else {
  cat("mu_gene comparison data frame is NULL.\n")
}

# --- Visualization of Recovery (Optional) ---

# Plot true vs. estimated for scalar parameters
if (!is.null(scalar_comparison) && nrow(scalar_comparison) > 0) {
  plot_scalar <- ggplot(scalar_comparison, aes(x = true, y = mean)) +
    geom_pointrange(aes(ymin = q5, ymax = q95)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~parameter, scales = "free") +
    labs(title = "Scalar Parameter Recovery", x = "True Value", y = "Posterior Mean (with 90% CI)") +
    theme_bw()
  print(plot_scalar)
}


# Plot true vs. estimated for a vector parameter (e.g., mu_gene)
if (!is.null(vector_comparison$mu_gene)) {
  plot_mu_gene <- ggplot(vector_comparison$mu_gene, aes(x = true, y = mean)) +
    geom_pointrange(aes(ymin = q5, ymax = q95), alpha=0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Recovery of mu_gene Effects (Standardized)", x = "True Value", y = "Posterior Mean (with 90% CI)") +
    theme_bw()
  print(plot_mu_gene)
}
# Plot true vs. estimated for a SD parameter random effects (e.g., sigma_gene)
if (!is.null(vector_comparison$sigma_gene)) {
  plot_sigma_gene <- ggplot(vector_comparison$sigma_gene, aes(x = true, y = mean)) +
    geom_pointrange(aes(ymin = q5, ymax = q95), alpha=0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Recovery of sigma_gene Effects (Standardized)", x = "True Value", y = "Posterior Mean (with 90% CI)") +
    theme_bw()
  print(plot_sigma_gene)
}


# Further checks: Density plots of posterior vs true value
# Example for sd_mu_indiv
draws_df <- as_draws_df(fit$draws())

if ("sd_mu_indiv" %in% names(draws_df)) {
  plot_sd_mu_indiv_density <- ggplot(draws_df, aes(x = sd_mu_indiv)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = true_params$sd_mu_indiv, color = "red", linetype = "dashed", size=1) +
    labs(title = "Posterior Density for sd_mu_indiv", x = "sd_mu_indiv",
         subtitle = paste("True value =", round(true_params$sd_mu_indiv, 3))) +
    theme_bw()
  print(plot_sd_mu_indiv_density)
}


cat("\nSimulation and analysis complete.\n")