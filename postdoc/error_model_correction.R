# heteroskedastic error-in-variables correlation estimator
hetero_cor <- function(x, y,
                       var_x,           # vector of measurement-error variances σ²_{x,i}
                       var_y,           # vector of measurement-error variances σ²_{y,i}
                       cov_xy = NULL,   # optional vector of error covariances τ_i (default 0)
                       ci = FALSE,      # logical: return 95 % bootstrap CI?
                       n_boot = 999,    # bootstrap replicates
                       seed = NULL) {
  
  # basic checks
  if (!(length(x) == length(y) &&
        length(x) == length(var_x) &&
        length(x) == length(var_y))) {
    stop("x, y, var_x, and var_y must have the same length")
  }
  n <- length(x)
  if (n < 3) stop("need at least 3 observations")
  
  # allow for optional correlated errors
  if (is.null(cov_xy)) cov_xy <- rep(0, n)
  if (length(cov_xy) != n) stop("cov_xy must have length n or be NULL")
  
  # sample moments (unbiased, denominator n-1)
  Sxx <- stats::var(x)
  Syy <- stats::var(y)
  Sxy <- stats::cov(x, y)
  
  # average measurement-error terms
  bar_sig2_x <- mean(var_x)
  bar_sig2_y <- mean(var_y)
  bar_tau     <- mean(cov_xy)
  
  # debias
  tilde_Sxx <- Sxx - bar_sig2_x
  tilde_Syy <- Syy - bar_sig2_y
  tilde_Sxy <- Sxy - bar_tau
  
  # guard against negative debiased variances
  if (tilde_Sxx <= 0 || tilde_Syy <= 0)
    stop("debiased variance ≤ 0; measurement noise overwhelms the signal")
  
  rho_hat <- tilde_Sxy / sqrt(tilde_Sxx * tilde_Syy)
  
  # optional bootstrap CI
  if (!ci) return(rho_hat)
  
  if (!is.null(seed)) set.seed(seed)
  boot_vals <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    Sxx_b <- stats::var(x[idx])
    Syy_b <- stats::var(y[idx])
    Sxy_b <- stats::cov(x[idx], y[idx])
    rho_b <- (Sxy_b - mean(cov_xy[idx])) /
      sqrt((Sxx_b - mean(var_x[idx])) *
             (Syy_b - mean(var_y[idx])))
    rho_b
  })
  
  ci_bounds <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
  
  list(rho = rho_hat,
       ci95 = as.numeric(ci_bounds),
       boot_values = boot_vals)
}

# latent signals with true correlation 0.8
n  <- 1E3
X0 <- rnorm(n, 0, 1)
r <- 0.7
Y0 <- r * X0 + sqrt(1 - r^2) * rnorm(n)
sd_e_x <- 1
sd_e_y <- 1

# heteroskedastic measurement noise
se_x <- abs(rnorm(n)) * sd_e_x          # each obs has its own s.e.
se_y <- abs(rnorm(n)) * sd_e_y 

x_obs <- X0 + rnorm(n, 0, se_x)
y_obs <- Y0 + rnorm(n, 0, se_y)

# naive and corrected correlations
cor_naive      <- cor(x_obs, y_obs)
cor_corrected  <- hetero_cor(x_obs, y_obs, se_x^2, se_y^2, ci = T)

cor(X0, Y0)
cor_naive
cor_corrected$ci95
