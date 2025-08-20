# run_student_t_mixture_prior.R
# Samples x and y from their priors and compares them via a Q‑Q plot (base R).

library(cmdstanr)
library(statmod)      # Gauss‑Laguerre quadrature
library(posterior)

# ---- inputs ----
df    <- 4
loc   <- 0
scale <- 1
K     <- 8          # mixture components
ndraw <- 20000       # Stan posterior draws = prior samples
disc_method <- c("laguerre", "quantile")[2]
# ----------------

# 1.  Build mixture parameters ------------------------------------------------

shape <- df / 2
rate  <- df / 2   # Gamma(shape, rate)

# 1.  Build mixture parameters --------------------------------------------
if (K == 1) {
  # Moment‑matched single Normal:  Var = σ²·df/(df−2)
  lambda <- (df - 2) / df   # precision = 1 / (df/(df−2))
  pi     <- 1
} else {
  if (identical(disc_method, "laguerre")) {
    # Gauss–Laguerre quadrature (statmod uses *scale*, not rate)
    scale_gamma <- 1 / rate
    quad   <- statmod::gauss.quad.prob(K, dist = "gamma",
                                       alpha = shape, beta = scale_gamma)
    lambda <- quad$nodes
    pi     <- quad$weights / sum(quad$weights)
    
  } else if (identical(disc_method, "quantile")) {
    # Equal‑probability (quantile) bins -----------------------------------
    edges  <- qgamma(seq(0, 1, length.out = K + 1), shape = shape, rate = rate)
    pi     <- rep(1 / K, K)
    lambda <- numeric(K)
    for (k in seq_len(K)) {
      lambda[k] <- integrate(function(l) l * dgamma(l, shape = shape, rate = rate),
                             lower = edges[k], upper = edges[k + 1])$value / pi[k]
    }
  } else {
    stop("Unknown discretisation method: ", disc_method)
  }
}

data_list <- list(K = K, df = df, loc = loc, scale = scale,
                  lambda = lambda, pi = pi)

# 2.  Compile & sample --------------------------------------------------------
mod  <- cmdstan_model("~/scripts/minor_scripts/postdoc/student_t_mixture.stan")
fit <- mod$sample(data = data_list,
                  chains = 4,
                  iter_warmup = 500,
                  iter_sampling = ndraw / 4,
                  refresh = 0)

draws <- as_draws_df(fit$draws())
x <- draws$x
y <- draws$y

# 3.  Q‑Q plot (base R) -------------------------------------------------------
qp <- 1:999/1000
qx <- quantile(x, qp)
qy <- quantile(y, qp)
plot(qx, qy, 
       xlab = "Student‑t quantiles", ylab = "Mixture quantiles",
       main = sprintf("Q‑Q plot, df = %g, K = %d", df, K), type = "l", col = 2)
abline(0, 1, lty = 2)
abline(v=-10:10, lty = 3, lwd = 0.5)

