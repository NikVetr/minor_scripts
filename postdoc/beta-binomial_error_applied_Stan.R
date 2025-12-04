# dependencies: cmdstanr
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::check_cmdstan_toolchain()
# cmdstanr::install_cmdstan() # if needed

library(cmdstanr)

set.seed(42)

# ---- simulate beta-binomial data ----
rbetabinom <- function(n, size, alpha, beta) {
  # draw theta ~ Beta(alpha, beta), then y ~ Binomial(size, theta)
  theta <- rbeta(n, alpha, beta)
  stats::rbinom(n, size, theta)
}

# true parameters
N <- 10
size <- rep(30L, N)
p_true <- 0.35
tau_true <- 20.0
alpha_true <- p_true * tau_true
beta_true  <- (1 - p_true) * tau_true
y <- rbetabinom(N, size, alpha_true, beta_true)

# sanity check
cat(sprintf("debug: mean(y/size)=%.4f (true p=%.4f)\n", mean(y/size), p_true))

# ---- write Stan model to a temp file ----
stan_code <- '
functions {
  // Rising Pochhammer in log-scale with log-parameter:
  real log_rising_factorial_log(real log_x, int m) {
    if (m <= 0) return 0.0;
    real acc = log_x;  // base case, ie j == 0
    for (j in 1:m - 1) {
      acc += log_sum_exp(log_x, log(j));
    }
    return acc;
  }

  real beta_binomial_log_lpmf(int y, int n, real log_alpha, real log_beta) {
    // lchoose is already on the log scale and exact for integers
    real lp = lchoose(n, y);

    lp += log_rising_factorial_log(log_alpha, y);           // lgamma(α+y)    - lgamma(α)
    lp += log_rising_factorial_log(log_beta,  n - y);       // lgamma(β+n−y)  - lgamma(β)

    real log_alpha_plus_beta = log_sum_exp(log_alpha, log_beta);
    lp -= log_rising_factorial_log(log_alpha_plus_beta, n);

    return lp;
  }
  
  real beta_binomial_stable_lpmf(int y, int n, real alpha, real beta) {
    // lchoose is already on the log scale and exact for integers
    real lp = lchoose(n, y);
    real log_alpha = log(alpha);
    real log_beta = log(beta);
    lp +=  beta_binomial_log_lpmf(y | n, log_alpha, log_beta);

    return lp;
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  array[N] int<lower=0> n;
  // prior switch: 1 = uniform(-200,200) on log_tau; 2 = half-normal(s_tau) on tau
  int<lower=1,upper=2> prior_type;
  real<lower=0> s_tau; // scale for half-normal on tau (used only if prior_type==2)
}

parameters {
  real logit_p;    // mean on logit scale
  real log_tau;    // log concentration
}

transformed parameters {
  real p = inv_logit(logit_p);
  real tau = exp(log_tau);
}

model {
  // priors on p (weakly informative, symmetric)
  logit_p ~ normal(0, 2.5);

  if (prior_type == 1) {
    // uniform(-200, 200) on log_tau
    target += uniform_lpdf(log_tau | -200, 200);
  } else {
    // half-normal(s_tau) prior on tau with Jacobian under log_tau parameterization
    // tau ~ normal(0, s_tau) T[0,]
    target += normal_lpdf(tau | 0, s_tau) - normal_lccdf(0 | 0, s_tau) + log_tau;
  }

  // likelihood
  {
    real log_a;
    real log_b;
    for (i in 1:N) {
      log_a = log(p) + log_tau;
      log_b = log1m(p) + log_tau;
      target += beta_binomial_stable_lpmf(y[i] | n[i], log_a, log_b);
    }
  }
}

generated quantities {
  // expose for quick checks
  real p_post = p;
  real tau_post = tau;
}
'

stan_file <- write_stan_file(stan_code)

# ---- compile ----
mod <- cmdstan_model(stan_file)

# ---- common data list ----
data_base <- list(
  N = N,
  y = y,
  n = as.array(size),
  s_tau = 20.0  # half-normal scale (roughly puts mass over tau in [0, ~60])
)

# ---- fit with uniform(-200, 200) on log_tau ----
data_unif <- data_base
data_unif$prior_type <- 1L

fit_unif <- mod$sample(
  data = data_unif,
  seed = 123,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  max_treedepth = 12, adapt_delta = 0.9
)

summ <- fit_unif$summary()
summ$ess_bulk
