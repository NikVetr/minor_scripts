# Load necessary libraries
library(cmdstanr)
library(posterior)  # For working with draws and computing ESS

#### untruncated models ####

# Define the Stan models as strings
centered_funnel <- "
data {
  real<lower=0> sd_y;  // standard deviation of y
  int<lower=1> n_x;  // length of x
}
parameters {
  real y;
  vector[n_x] x;
}
model {
  y ~ normal(0, sd_y);        // Prior on y
  x ~ normal(0, exp(y));     // Prior on x given y (centered)
}
"

noncentered_funnel <- "
data {
  real<lower=0> sd_y;  // standard deviation of y
  int<lower=1> n_x;  // length of x
}
parameters {
  real y_raw;
  vector[n_x] x_raw;
}
transformed parameters {
  real y = y_raw * sd_y;
  vector[n_x] x = x_raw * exp(y);  // Reparameterized x
}
model {
  y_raw ~ std_normal();   // Prior on y
  x_raw ~ std_normal();         // Standard normal for z (non-centered)
}
"

quantile_funnel <- "
data {
  real<lower=0> sd_y;  // standard deviation of y
  int<lower=1> n_x;  // length of x
}
parameters {
  real<lower=0, upper=1> y_q;  // Uniform(0, 1) sample for y
  vector<lower=0, upper=1>[n_x] x_q;  // Uniform(0, 1) samples for x
}
transformed parameters {
  real y = sd_y * inv_Phi(y_q);  // Transform y_q to normal(0, sd_y)
  vector[n_x] x;  // Transform x_q to normal(0, exp(y))
  for (i in 1:n_x) {
    x[i] = inv_Phi(x_q[i]) * exp(y);  // Apply inverse CDF transformation
  }
}
model {
  // Implicitly using Uniform(0, 1) priors for y_q and x_q
}
"


# Compile the Stan models
centered_model <- cmdstan_model(write_stan_file(centered_funnel))
noncentered_model <- cmdstan_model(write_stan_file(noncentered_funnel))
quantile_model <- cmdstan_model(write_stan_file(quantile_funnel))

# Define data to pass to models
funnel_data <- list(sd_y = 1,
                    n_x = 100)

# Fit models
niter <- 5E3
nsamp <- 1E3
nref <- 20
fit_centered <- centered_model$sample(
  data = funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
  
)
fit_noncentered <- noncentered_model$sample(
  data = funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)
fit_quantile <- quantile_model$sample(
  data = funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter/2, iter_sampling = niter, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

# Extract draws
draws_centered <- data.frame(as_draws_df(fit_centered$draws()))
draws_noncentered <- data.frame(as_draws_df(fit_noncentered$draws()))
draws_quantile <- data.frame(as_draws_df(fit_quantile$draws()))

# Sample from target distributions directly
n <- niter * 1E2
y <- rnorm(n, 0, funnel_data$sd_y)
x <- rnorm(n, 0, sd = exp(y))

# Inspect q-q plots
par(mfrow = c(2,3), mar = c(5,5,2,2))

#for y
qqplot(
  y, draws_centered$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Centered Parameterization Posterior Samples",
  main = "Q-Q Plot: Y, Centered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  y, draws_noncentered$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Noncentered Parameterization Posterior Samples",
  main = "Q-Q Plot: Y, Noncentered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  y, draws_quantile$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Quantile Parameterization Posterior Samples",
  main = "Q-Q Plot: Y, Quantile",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

#for x
qqplot(
  asinh(x), asinh(draws_centered$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Non-Centered Parameterization Posterior Samples",
  main = "Q-Q Plot: asinh(X), Centered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  asinh(x), asinh(draws_noncentered$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Centered Parameterization Posterior Samples",
  main = "Q-Q Plot: asinh(x), Noncentered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  asinh(x), asinh(draws_quantile$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Quantile Parameterization Posterior Samples",
  main = "Q-Q Plot: asinh(x), Noncentered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

# Examine MCMC diagnostics
summ_centered <- fit_centered$summary()
summ_noncentered <- fit_noncentered$summary()
summ_quantile <- fit_quantile$summary()

mean(fit_centered$time()$chains$total)
print(head(summ_centered[order(summ_centered$ess_bulk),], 3))

mean(fit_noncentered$time()$chains$total)
print(head(summ_noncentered[order(summ_noncentered$ess_bulk),], 3))

mean(fit_quantile$time()$chains$total)
print(head(summ_quantile[order(summ_quantile$ess_bulk),], 3))



#### truncated models ####

# Define the Stan models as strings
trunc_centered_funnel <- "
data {
  real<lower=0> sd_y;  // Standard deviation of y
  int<lower=1> n_x;  // Length of x
  vector[2] y_b;  // Bounds for truncation of y (lower and upper)
  vector[2] x_b;  // Bounds for truncation of x (lower and upper)
}
parameters {
  real<lower=y_b[1], upper=y_b[2]> y;  // Truncated y
  vector<lower=x_b[1], upper=x_b[2]>[n_x] x;  // Truncated x
}
model {
  y ~ normal(0, sd_y);  // Prior on y
  x ~ normal(0, exp(y));  // Prior on x given y (centered)
}
"

trunc_noncentered_funnel <- "
data {
  real<lower=0> sd_y;  // Standard deviation of y
  int<lower=1> n_x;  // Length of x
  vector[2] y_b;  // Bounds for truncation of y (lower and upper)
  vector[2] x_b;  // Bounds for truncation of x (lower and upper)
}
transformed data {
  real lb_y = (y_b[1] / sd_y);
  real ub_y = (y_b[2] / sd_y);
}
parameters {
  real<lower=lb_y, upper=ub_y> y_raw;  // Truncated raw y
  vector<lower=x_b[1] / exp(y_raw * sd_y), 
         upper=x_b[2] / exp(y_raw * sd_y)>[n_x] x_raw;  // Raw x values
}
transformed parameters {
  real y = y_raw * sd_y;
  vector[n_x] x = x_raw * exp(y);
}
model {
  y_raw ~ std_normal();  // Truncated normal for y
  x_raw ~ std_normal();  // Standard normal for x_raw
}
"

trunc_quantile_funnel <- "
data {
  real<lower=0> sd_y;  // Standard deviation of y
  int<lower=1> n_x;  // Length of x
  vector[2] y_b;  // Bounds for truncation of y (lower and upper)
  vector[2] x_b;  // Bounds for truncation of x (lower and upper)
}
transformed data {
  real lb_y = Phi(y_b[1] / sd_y);
  real ub_y = Phi(y_b[2] / sd_y);
}
parameters {
  real<lower=0, upper=1> y_q;  // Uniform(0, 1) sample for y
  vector<lower=0, upper=1>[n_x] x_q;  // Uniform(0, 1) samples for x
}
transformed parameters {
  real y = sd_y * inv_Phi(y_q * (ub_y - lb_y) + lb_y);
  real lb_x = Phi(x_b[1] / exp(y));
  real ub_x = Phi(x_b[2] / exp(y));
  vector[n_x] x = inv_Phi(x_q * (ub_x - lb_x) + lb_x) * exp(y);
}
model {
  // Implicitly using Uniform(0, 1) priors for y_q and x_q
}
"


trunc_zquantile_funnel <- "
data {
  real<lower=0> sd_y;  // Standard deviation of y
  int<lower=1> n_x;  // Length of x
  vector[2] y_b;  // Bounds for truncation of y (lower and upper)
  vector[2] x_b;  // Bounds for truncation of x (lower and upper)
}
transformed data {
  real lb_y = Phi(y_b[1] / sd_y);
  real ub_y = Phi(y_b[2] / sd_y);
}
parameters {
  real y_z;
  vector[n_x] x_z;
}
transformed parameters {
  real y = sd_y * inv_Phi(Phi(y_z) * (ub_y - lb_y) + lb_y);
  real lb_x = Phi(x_b[1] / exp(y));
  real ub_x = Phi(x_b[2] / exp(y));
  vector[n_x] x = inv_Phi(Phi(x_z) * (ub_x - lb_x) + lb_x) * exp(y);
}
model {
  y_z ~ std_normal();
  x_z ~ std_normal();
}
"

# Compile the Stan models
trunc_centered_model <- cmdstan_model(write_stan_file(trunc_centered_funnel))
trunc_noncentered_model <- cmdstan_model(write_stan_file(trunc_noncentered_funnel))
trunc_quantile_model <- cmdstan_model(write_stan_file(trunc_quantile_funnel))
trunc_zquantile_model <- cmdstan_model(write_stan_file(trunc_zquantile_funnel))

# Define data to pass to models
trunc_funnel_data <- list(sd_y = 1,
                          n_x = 100,
                          y_b = c(-2,3),
                          x_b = c(-10, 40)
)

# Fit models
niter <- 5E3
nsamp <- 5E3
nref <- 20

fit_trunc_centered <- trunc_centered_model$sample(
  data = trunc_funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter, iter_sampling = niter/2, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

fit_trunc_noncentered <- trunc_noncentered_model$sample(
  data = trunc_funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter, iter_sampling = niter/2, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

fit_trunc_quantile <- trunc_quantile_model$sample(
  data = trunc_funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter, iter_sampling = niter/2, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

fit_trunc_zquantile <- trunc_zquantile_model$sample(
  data = trunc_funnel_data, chains = 4, parallel_chains = 4, 
  iter_warmup = niter, iter_sampling = niter/2, thin = ceiling(niter / nsamp),
  refresh = ceiling(niter / nref)
)

# Extract draws
draws_trunc_centered <- data.frame(as_draws_df(fit_trunc_centered$draws()))
draws_trunc_noncentered <- data.frame(as_draws_df(fit_trunc_noncentered$draws()))
draws_trunc_quantile <- data.frame(as_draws_df(fit_trunc_quantile$draws()))
draws_trunc_zquantile <- data.frame(as_draws_df(fit_trunc_zquantile$draws()))

# Sample from target distributions directly
n <- niter * 1E2
yf <- rnorm(n, 0, trunc_funnel_data$sd_y)
y <- yf[yf > trunc_funnel_data$y_b[1] & yf < trunc_funnel_data$y_b[2]]
xf <- rnorm(length(y), 0, sd = exp(y))
x <- xf[xf > trunc_funnel_data$x_b[1] & xf < trunc_funnel_data$x_b[2]]

# Examine MCMC diagnostics
summ_trunc_centered <- fit_trunc_centered$summary(c("x", "y"))
summ_trunc_noncentered <- fit_trunc_noncentered$summary(c("x_raw", "y_raw"))
summ_trunc_quantile <- fit_trunc_quantile$summary(c("x_q", "y_q"))
summ_trunc_zquantile <- fit_trunc_zquantile$summary(c("x_z", "y_z"))

mean_chain_times <- list(
  centered = mean(fit_trunc_centered$time()$chains$total),
  noncentered = mean(fit_trunc_noncentered$time()$chains$total),
  quantile = mean(fit_trunc_quantile$time()$chains$total),
  zquantile = mean(fit_trunc_zquantile$time()$chains$total)
)

summ_trunc_centered$ess_bulk_per_s <- summ_trunc_centered$ess_bulk / 
  mean_chain_times$centered
summ_trunc_centered$ess_tail_per_s <- summ_trunc_centered$ess_tail / 
  mean_chain_times$centered

summ_trunc_noncentered$ess_bulk_per_s <- summ_trunc_noncentered$ess_bulk / 
  mean_chain_times$noncentered
summ_trunc_noncentered$ess_tail_per_s <- summ_trunc_noncentered$ess_tail / 
  mean_chain_times$noncentered

summ_trunc_quantile$ess_bulk_per_s <- summ_trunc_quantile$ess_bulk / 
  mean_chain_times$quantile
summ_trunc_quantile$ess_tail_per_s <- summ_trunc_quantile$ess_tail / 
  mean_chain_times$quantile

summ_trunc_zquantile$ess_bulk_per_s <- summ_trunc_zquantile$ess_bulk / 
  mean_chain_times$zquantile
summ_trunc_zquantile$ess_tail_per_s <- summ_trunc_zquantile$ess_tail / 
  mean_chain_times$zquantile

#get the ranges to make the axes the same
range_ess_bulk <- range(
  c(summ_trunc_centered$ess_bulk_per_s,
    summ_trunc_noncentered$ess_bulk_per_s,
    summ_trunc_quantile$ess_bulk_per_s,
    summ_trunc_zquantile$ess_bulk_per_s)
)

range_ess_tail <- range(
  c(summ_trunc_centered$ess_tail_per_s,
    summ_trunc_noncentered$ess_tail_per_s,
    summ_trunc_quantile$ess_tail_per_s,
    summ_trunc_zquantile$ess_tail_per_s)
)

#### truncated plotting ####
par(mfrow = c(4,4), mar = c(5,5,2,2))

# Inspect q-q plots

#for y
qqplot(
  y, draws_trunc_centered$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Centered Parameterization",
  main = "Q-Q Plot (Truncated): Y, Centered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  y, draws_trunc_noncentered$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Noncentered Parameterization",
  main = "Q-Q Plot (Truncated): Y, Noncentered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  y, draws_trunc_quantile$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Quantile Parameterization",
  main = "Q-Q Plot (Truncated): Y, Quantile",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  y, draws_trunc_zquantile$y, 
  xlab = "Target Distribution Samples", 
  ylab = "Z-Quantile Parameterization",
  main = "Q-Q Plot (Truncated): Y, Quantile",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

#for x
qqplot(
  asinh(x), asinh(draws_trunc_centered$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Centered Parameterization",
  main = "Q-Q Plot (Truncated): asinh(X), Centered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  asinh(x), asinh(draws_trunc_noncentered$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Noncentered Parameterization",
  main = "Q-Q Plot (Truncated): asinh(X), Noncentered",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  asinh(x), asinh(draws_trunc_quantile$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Quantile Parameterization",
  main = "Q-Q Plot (Truncated): asinh(X), Quantile",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

qqplot(
  asinh(x), asinh(draws_trunc_zquantile$x.1.), 
  xlab = "Target Distribution Samples", 
  ylab = "Quantile Parameterization",
  main = "Q-Q Plot (Truncated): asinh(X), Z-Quantile",
  type = "l", lwd = 2
)
abline(0, 1, col = "red", lty = 2)  # Add a 45-degree reference line

#for actual sampled distributions
nscatter <- 1E3
subsamp_inds <- round(seq(1, (niter*2), length.out = nscatter))

plot(
  draws_trunc_centered$x.1.[subsamp_inds],
  draws_trunc_centered$y[subsamp_inds],
  xlab = "Sampled X Variable", 
  ylab = "Sampled Y Variable",
  main = "Centered Parameterization",
  col = adjustcolor(1,0.1)
)

plot(
  draws_trunc_noncentered$x_raw.1.[subsamp_inds],
  draws_trunc_noncentered$y_raw[subsamp_inds],
  xlab = "Sampled X Variable", 
  ylab = "Sampled Y Variable",
  main = "Noncentered Parameterization",
  col = adjustcolor(1,0.1)
)

plot(
  draws_trunc_quantile$x_q.1.[subsamp_inds],
  draws_trunc_quantile$y_q[subsamp_inds],
  xlab = "Sampled X Variable", 
  ylab = "Sampled Y Variable",
  main = "Quantile Parameterization",
  col = adjustcolor(1,0.1)
)

plot(
  draws_trunc_zquantile$x_z.1.[subsamp_inds],
  draws_trunc_zquantile$y[subsamp_inds],
  xlab = "Sampled X Variable", 
  ylab = "Sampled Y Variable",
  main = "Z-Quantile Parameterization",
  col = adjustcolor(1,0.1)
)

#also scatters of ESS / second 
plot(
  summ_trunc_centered$ess_bulk_per_s,
  summ_trunc_centered$ess_tail_per_s,
  xlab = "ESS Bulk / second", 
  ylab = "ESS Tail / second",
  main = "Centered Parameterization",
  col = adjustcolor(1,0.1), log = "xy",
  xlim = range_ess_bulk, ylim = range_ess_tail
)

plot(
  summ_trunc_noncentered$ess_bulk_per_s,
  summ_trunc_noncentered$ess_tail_per_s,
  xlab = "ESS Bulk / second", 
  ylab = "ESS Tail / second",
  main = "Noncentered Parameterization",
  col = adjustcolor(1,0.1), log = "xy",
  xlim = range_ess_bulk, ylim = range_ess_tail
)


plot(
  summ_trunc_quantile$ess_bulk_per_s,
  summ_trunc_quantile$ess_tail_per_s,
  xlab = "ESS Bulk / second", 
  ylab = "ESS Tail / second",
  main = "Quantile Parameterization",
  col = adjustcolor(1,0.1), log = "xy",
  xlim = range_ess_bulk, ylim = range_ess_tail
)

plot(
  summ_trunc_zquantile$ess_bulk_per_s,
  summ_trunc_zquantile$ess_tail_per_s,
  xlab = "ESS Bulk / second", 
  ylab = "ESS Tail / second",
  main = "Z-Quantile Parameterization",
  col = adjustcolor(1,0.1), log = "xy",
  xlim = range_ess_bulk, ylim = range_ess_tail
)

#compare just z-quantile and quantile
par(mfrow = c(1,2))
plot(summ_trunc_quantile$ess_bulk_per_s, 
     summ_trunc_zquantile$ess_bulk_per_s,
     xlim = range(summ_trunc_quantile$ess_bulk_per_s,
                  summ_trunc_zquantile$ess_bulk_per_s),
     ylim = range(summ_trunc_quantile$ess_bulk_per_s,
                  summ_trunc_zquantile$ess_bulk_per_s),
     xlab = "ESS Bulk / s (Quantile Parameterization)",
     ylab = "ESS Bulk / s (Z-Quantile Parameterization)") 
abline(0,1,lty=2,col=2,lwd=2)

plot(summ_trunc_quantile$ess_tail_per_s, 
     summ_trunc_zquantile$ess_tail_per_s,
     xlim = range(summ_trunc_quantile$ess_tail_per_s,
                  summ_trunc_zquantile$ess_tail_per_s),
     ylim = range(summ_trunc_quantile$ess_tail_per_s,
                  summ_trunc_zquantile$ess_tail_per_s),
     xlab = "ESS Tail / s (Quantile Parameterization)",
     ylab = "ESS Tail / s (Z-Quantile Parameterization)")
abline(0,1,lty=2,col=2,lwd=2)

