library(cmdstanr)
library(posterior)

# specify stan model
stan_model_code <- '
data {
  int<lower=0> count;
  int<lower=count> total;
  real<lower=0> sd_logconc;
}

parameters {
  real log_concentration;
}

model {
  log_concentration ~ normal(0, sd_logconc);
}

generated quantities {
  // Convert to alpha, beta
  real shape = exp(log_concentration) / 2;
  
  // Compare binomial PMF at p = 0.5 vs. beta-binomial
  real log_binomial_05      = binomial_lpmf(count | total, 0.5);
  real log_beta_binomial    = beta_binomial_lpmf(count | total, shape, shape);
  real diff_lprob    = log_beta_binomial - log_binomial_05;
}
'
mod_path <- "~/beta-binomial_error.stan"
writeLines(stan_model_code, mod_path)
mod <- cmdstan_model(mod_path)

# specify data
stan_data <- list(count = 57, total = 117, sd_logconc = 50)

# fit model
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 1,   # for demonstration
  iter_warmup = 500,
  iter_sampling = 2000
)
fit$summary()

# extract and inspect samples
samps <- data.frame(as_draws_df(fit$draws()))
asinh_diff_logp <- asinh(samps$diff_lprob)
par(mar = c(6,6,2,2))
plot(samps$log_concentration, asinh_diff_logp,
     xlab = "log_concentration", 
     ylab = "diff_lprob", yaxt = "n")
axis(2, at = asinh(round(sinh(pretty(asinh_diff_logp)), 0)), 
     round(sinh(pretty(asinh_diff_logp)), 0))
abline(h=0, col = 2, lty = 2)

#mark where it goes wild
wild_lprob_at <- min(samps$log_concentration[asinh_diff_logp > 0.01])
abline(v=wild_lprob_at, lty = 3)
text(x = wild_lprob_at, y = par("usr")[4], labels = round(wild_lprob_at, 2), pos = 3, xpd = NA)
