library(cmdstanr)
library(posterior)
library(mvtnorm)

# Save the Stan code to a file
stan_code <- "
functions {
  matrix cholesky_corr_constrain_tfp_super_lp (vector y, int K) {
    matrix[K, K] L = identity_matrix(K);
    int counter = 1;
    
    for (i in 2 : K) {
        row_vector[i - 1] y_star = y[counter:counter + i - 2]';
          real dsy = dot_self(y_star);
          real alpha_r = 1 / (dsy  + 1);
          real gamma = sqrt( (1 - alpha_r^2) / dsy );
        L[i, : i] = append_col(gamma * y_star, alpha_r);
        target += -i * log1p(dsy);
        counter += i - 1;
      }
    return L;
  }
  
  matrix cholesky_corr_constrain_tfp_lp (vector y, int K) {
    matrix[K, K] L = identity_matrix(K);
    int counter = 1;
    real s;

    for (i in 2 : K) {
        L[i, 1:i - 1] = y[counter:counter + i - 2]';
        counter += i - 1;
        s = norm2(L[i,  : i]);
        L[i,  : i] = L[i,  : i] / s;
        target += - (i + 1) * log(s);
      }
    return L;
  }
}
data {
  int<lower=0> N;
  int bool_tfp;
  int bool_stan;
  real<lower=0> eta;
  int nobs;
  matrix[nobs, N] x;
  int bool_ll;

}
parameters {
 vector[choose(N, 2)] y;
 cholesky_factor_corr[N] L_stan;
}
transformed parameters {
 matrix[N, N] L;
 if(bool_tfp == 1){
  L = cholesky_corr_constrain_tfp_lp(y, N);
 } else {
  L = cholesky_corr_constrain_tfp_super_lp(y, N);
 }
}
model {

  L_stan ~ lkj_corr_cholesky(eta);
  if(bool_ll == 1){
    for(i in 1:nobs) {
      x[i] ~ multi_normal_cholesky(rep_vector(0, N), L_stan);
    }
  }

  L ~ lkj_corr_cholesky(eta);
  if(bool_ll == 1){
    for(i in 1:nobs) {
      x[i] ~ multi_normal_cholesky(rep_vector(0, N), L);
    }
  }

}
generated quantities {
  // Reconstruct correlation matrix from lower cholesky factor
  matrix[N,N] C;
  C = multiply_lower_tri_self_transpose(L);
  matrix[N,N] C_stan;
  C_stan = multiply_lower_tri_self_transpose(L_stan);
}
"

# Compile the mod
mod <- cmdstan_model(write_stan_file(stan_code))

# define dimension of correlation matrix + lkj prior
N <- 15
eta <- 2.0  # LKJ prior shape parameter

# Create covariance matrix: 1s on diagonal, r_ij on off-diagonals
r_ij <- 0.5
Sigma <- matrix(r_ij, nrow = N, ncol = N)
diag(Sigma) <- 1.0

# Simulate multivariate normal data
nobs <- 1E2
x <- rmvnorm(nobs, mean = rep(0, N), sigma = Sigma)

# specify data list to pass to Stan
data_list <- list(
  N = N,
  bool_tfp = 1,
  bool_stan = 1,
  eta = eta,
  nobs = nobs,
  x = x, 
  bool_ll = 1
)

# Fit the mod
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 500,
  refresh = 100
)

# Print summary of results
summ_tfp <- fit$summary("C")
print(summ[order(summ_tfp$ess_bulk),])
print(summ[order(summ_tfp$rhat, decreasing = T),])

summ_stan <- fit$summary("C_stan")
print(summ[order(summ_stan$ess_bulk),])
print(summ[order(summ_stan$rhat, decreasing = T),])

# Extract samples
samps <- data.frame(as_draws_df(fit$draws()))

indices <- c(3,4)
samp1 <- samps[, paste0("C.", indices[1], ".", indices[2], ".")]
samp2 <- samps[, paste0("C_stan.", indices[1], ".", indices[2], ".")]

hist(samp1, breaks = -10:10/10)
hist(samp2, add = T, 
     col = adjustcolor(2,0.5), breaks = -10:10/10)

# Q-Q plot
quantiles <- 1:99/100
qs1 <- quantile(samp1, quantiles)
qs2 <- quantile(samp2, quantiles)
plot(qs1, qs2, type = "l")
abline(0, 1, col = "red", lwd = 2)

