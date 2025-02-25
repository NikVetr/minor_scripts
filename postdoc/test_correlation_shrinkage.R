# Load required libraries
library(cmdstanr)
library(mvtnorm)
library(ggplot2)
library(dplyr)

# Parameters for simulation
n <- 10
p <- 10
r_ij <- 0.5

# Create covariance matrix: 1s on diagonal, 0.5 on off-diagonals
Sigma <- matrix(r_ij, nrow = p, ncol = p)
diag(Sigma) <- 1.0

# Simulate multivariate normal data
x_full <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
samp_Sigma <- cor(x_full)
samp_Sigma[1,2]

#### models ####
# Define Stan model as a string
stan_model_lkj <- "
data {
  int<lower=1> n;          // number of observations
  int<lower=1> p;          // dimension
  real<lower=-1, upper=1> r;          // other correlations
  matrix[n, p] x;          // data
}
parameters {
  cholesky_factor_corr[p] Lcorr;       // Cholesky factor of correlation
}
model {
  // priors
  Lcorr ~ lkj_corr_cholesky(1);
  
  // likelihood
  for(i in 1:n) {
    x[i] ~ multi_normal_cholesky(rep_vector(0, p), Lcorr);
  }
}
generated quantities {
  // Reconstruct correlation matrix from Lcorr
  matrix[p,p] cor_mat;
  cor_mat = multiply_lower_tri_self_transpose(Lcorr);
  real cor_12 = cor_mat[1,2];
}
"

stan_model_luo <- "
functions {
  //--- Utility: fill a p x p symmetric matrix G from gamma (off-diagonal)
  matrix build_G(vector gamma_raw, int p) {
    matrix[p,p] G;
    // zero-initialize
    G = rep_matrix(0.0, p, p);

    // fill strictly-lower triangle from gamma_raw
    // We'll do a simple index scheme:
    int idx = 1;
    for(i in 1:(p-1)) {
      for(j in (i+1):p) {
        G[j,i] = gamma_raw[idx]; 
        G[i,j] = gamma_raw[idx]; 
        idx += 1;
      }
    }
    return G;
  }

  //--- Utility: matrix exponential via built-in Stan function
  //    \"matrix_exp\" is available in Stan 2.32+ (or \"matrix_exp_multiply\" if older).
  //    We'll just use matrix_exp (an alternative is \"expm\" from Math library).
  
  //--- Key iterative function:
  //    Adjust diagonal of G so that exp(G) has diag ~ 1
  matrix fix_diag_logC(matrix G, int max_iter, real tol) {
    int p = rows(G);
    matrix[p,p] Gtemp = G;
    for(iter in 1:max_iter) {
      matrix[p,p] Ctemp = matrix_exp(Gtemp);
      vector[p] d = diagonal(Ctemp);
      // measure how far from 1
      real err = max(abs(d - rep_vector(1.0, p)));
      if(err < tol) break;

      // update Gii <- Gii - log(Ctemp[ii])
      for(i in 1:p) {
        Gtemp[i,i] = Gtemp[i,i] - log(d[i]);
      }
    }
    return Gtemp;
  }
}

data {
  int<lower=1> n;          // number of observations
  int<lower=1> p;          // dimension
  matrix[n, p] x;          // data
}

parameters {
  // gamma_raw: length = p*(p-1)/2
  // These are the off-diagonal entries of log C
  vector[ ( (p*(p-1)) %/% 2 ) ] gamma_raw;
}

transformed parameters {
  matrix[p,p] G;            // logC before diagonal fix
  matrix[p,p] G_fixed;      // logC after diagonal fix
  matrix[p,p] C_samp;       // final correlation matrix

  // 1) Build the initial G from gamma_raw (diagonal=0)
  G = build_G(gamma_raw, p);

  // 2) Perform iterative fix to get diag(exp(G))=1
  //    We'll do, say, 8 iterations. That often suffices for moderate p.
  G_fixed = fix_diag_logC(G, 8, 1e-8);

  // 3) Now exponentiate to get final correlation matrix
  C_samp = matrix_exp(G_fixed);
}

model {
  // 1) Put a prior on gamma_raw:
  //    For demonstration, use i.i.d normal(0, 1) for each off-diagonal
  gamma_raw ~ normal(0, 1);

  // 2) Likelihood: x[i] ~ MVN(0, C_samp)
  //    We can use multi_normal_lpdf with C_samp
  for(i in 1:n) {
    x[i] ~ multi_normal(rep_vector(0.0, p), C_samp);
  }
}

generated quantities {
  real cor_12;
  cor_12 = C_samp[1,2];
}
"

stan_model_trace <- "
data {
  int<lower=1> n;          // number of observations
  int<lower=1> p;          // dimension
  real<lower=-1, upper=1> r;          // other correlations
  matrix[n, p] x;          // data
  int<lower=0> incl_ll;          // toggle to include the log likelihood
}
parameters {
  cholesky_factor_corr[p] Lcorr;       // Cholesky factor of correlation
}
transformed parameters {
  vector[p] diag_Lcorr = diagonal(Lcorr);
  real scaled_trace = (sum(diag_Lcorr) - 1) / p;
}
model {
  // priors
  scaled_trace ~ beta(p,1);  //in contrast, beta(1,p) would weigh towards the identity

  // likelihood
  if(incl_ll == 1){
    for(i in 1:n) {
      x[i] ~ multi_normal_cholesky(rep_vector(0, p), Lcorr);
    }  
  }

}
generated quantities {
  // Reconstruct correlation matrix from Lcorr
  matrix[p,p] cor_mat;
  cor_mat = multiply_lower_tri_self_transpose(Lcorr);
  real cor_12 = cor_mat[1,2];
}
"

# Define Stan model as a string
stan_model_1r <- "
data {
  int<lower=1> n;          // number of observations
  int<lower=1> p;          // dimension
  real<lower=-1, upper=1> r;          // other correlations
  matrix[n, p] x;          // data
}
parameters {
  real<lower=-1, upper=1> cor_12;       // correlation in the 1,2 position
}
model {
  // priors
  // implicit uniform prior on cor_12 over support
  
  matrix[p,p] R;
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        R[i,j] = 1;
      } else if(i==1 && j==2){
        R[i,j] = cor_12;
      } else if(i==2 && j==1){
        R[i,j] = cor_12;
      } else{
        R[i,j] = r;
      }
    }
  }
  
  // likelihood
  for(i in 1:n) {
    x[i] ~ multi_normal(rep_vector(0, p), R);
  }
}
"


# Define Stan model as a string
stan_model_xr <- "
data {
  int<lower=1> n;          // number of observations
  int<lower=1> p;          // dimension
  real<lower=-1, upper=1> r;          // other correlations
  matrix[n, p] x;          // data
}
parameters {
  real<lower=-1, upper=1> cor_12;       // correlation in the 1,2 position
  vector<lower=-1, upper=1>[p-2] cor_1x;       // correlation in the 1,2 position
}
model {
  // priors
  // implicit uniform prior on cor_12 over support

  matrix[p,p] R;
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        R[i,j] = 1;
      } else if(i==1){
        if(j==2){
          R[i,j] = cor_12;
        } else {
          R[i,j] = cor_1x[j-2];
        }
      } else if(j==1){
        if(i==2){
          R[i,j] = cor_12;
        } else {
          R[i,j] = cor_1x[i-2];
        }
      } else{
        R[i,j] = r;
      }
    }
  }

  // likelihood
  for(i in 1:n) {
    x[i] ~ multi_normal(rep_vector(0, p), R);
  }
}
"

# Compile the Stan model
model_to_use <- 5
stan_models <- list(stan_model_lkj,
                    stan_model_1r,
                    stan_model_xr,
                    stan_model_trace,
                    stan_model_luo
)
stan_model <- stan_models[[model_to_use]]
mod <- cmdstan_model(write_stan_file(stan_model))

#### fit ####
# Fit the model for dimensions 2 through 10
posterior_samples <- list()

for(d in 2:p) {
  print(d)
  
  # Subset the data to the first d columns
  x_sub <- x_full[, 1:d, drop = FALSE]
  
  # Prepare data list for Stan
  data_list <- list(
    n = n,
    p = d,
    r = r_ij,
    x = x_sub
    # , incl_ll = 0
  )
  
  # Fit the model
  fit <- mod$sample(
    data = data_list,
    seed = 123,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 2000,
    refresh = 0
  )
  
  # Extract posterior draws of cor_12
  cor12_draws <- fit$draws("cor_12", inc_warmup = FALSE)
  cor12_vec <- as.vector(posterior::as_draws_matrix(cor12_draws))
  
  posterior_samples[[d]] <- cor12_vec
  
  # samps <- fit$draws("scaled_trace", inc_warmup = FALSE)
  # samps <- as.data.frame(posterior::as_draws_matrix(samps))
  # hist(samps$scaled_trace)
  #hmmm clearly not the specified distribution over the trace
  
}

# Combine all posterior draws into one data frame for plotting
plot_data <- do.call(rbind, lapply(seq_along(posterior_samples), function(d) {
  if(!is.null(posterior_samples[[d]])) {
    data.frame(
      cor_12 = posterior_samples[[d]],
      dimension = factor(d, levels = 2:p)
    )
  }
})) 

# Compute the sample correlation for the first two dimensions
sample_cor <- cor(x_full[, 1], x_full[, 2])

#### visualize ####
# Plot overlapping histograms with vertical lines for true and sample correlations
ggplot(plot_data, aes(x = cor_12, fill = dimension)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 100) +
  xlim(-1, 1) +
  geom_vline(xintercept = r_ij, color = "red", linetype = "dashed", size = 1, 
             alpha = 0.8, show.legend = FALSE) +  # True correlation
  geom_vline(xintercept = sample_cor, color = "blue", linetype = "dotted", size = 1, 
             alpha = 0.8, show.legend = FALSE) + # Sample correlation
  labs(
    title = "",
    x = "Posterior for Correlation Parameter",
    y = "Density"
  ) +
  annotate("text", x = r_ij+0.075, y = 0.3, label = "True r_ij", color = "red", 
           vjust = -27, size = 4) +
  annotate("text", x = sample_cor-0.125, y = 0.3, label = "Sample Cor", color = "blue", 
           vjust = -27, size = 4) +
  coord_cartesian(clip = "off") +  # Ensures text can appear outside the plot area
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 20, 10))
