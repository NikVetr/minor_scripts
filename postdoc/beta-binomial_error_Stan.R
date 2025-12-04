library(cmdstanr)
library(posterior)

# specify stan model
stan_model_code <- '
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
  int<lower=0> count;
  int<lower=count> total;
  real<lower=0> sd_logconc;
  int<lower=1> version;
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
  real log_beta_binomial;
  if(version == 1){
    log_beta_binomial    = beta_binomial_lpmf(count | total, shape, shape);
  } else if(version == 2){
    log_beta_binomial    = beta_binomial_stable_lpmf(count | total, shape, shape);
  } else if(version == 3){
    log_beta_binomial    = beta_binomial_lpmf(count | total, shape, shape);
  } else if(version == 4){
    log_beta_binomial    = beta_binomial_lpmf(count | total, shape, shape);
  }
  real diff_lprob    = log_beta_binomial - log_binomial_05;
}
'
mod_path <- "~/beta-binomial_error.stan"
writeLines(stan_model_code, mod_path)
mod <- cmdstan_model(mod_path)

# specify data
stan_data <- list(count = 57, 
                  total = 117, 
                  sd_logconc = 50, 
                  version = 2)

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


# plot(samps$log_concentration, samps$log_beta_binomial,
#      xlab = "log_concentration", 
#      ylab = "log_beta_binomial")
# abline(v = 0)
# hist(samps$log_concentration[samps$log_beta_binomial == max(samps$log_beta_binomial)], breaks = 100)
# min(samps$log_concentration[samps$log_beta_binomial == max(samps$log_beta_binomial)])

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

#compare to r-versions
dbetabinom_log <- function(k, size, alpha, beta) {
  lchoose(size, k) + lbeta(alpha + k, beta + size - k) - lbeta(alpha, beta)
}

# numerically stable beta-binomial log pmf via log rising factorials
# computes S(x, m) = log Γ(x+m) - log Γ(x) using m*log(x) + sum log1p(i/x)
# vectorized over k, with O(n) precomputation and O(1) per k lookups

dbetabinom_log_stable <- function(k, size, alpha, beta, check = TRUE, debug = FALSE) {
  n <- as.integer(size)
  k <- as.integer(k)
  if (length(k) == 0) return(numeric(0))
  if (any(k < 0 | k > n)) stop("all k must be in 0..size")
  
  # helper: prefix sums of log1p(i/x) for i=0..m-1, with prefix[0]=0
  prefix_log1p <- function(x, m) {
    if (m == 0) return(0)
    # i/x can be very small when x is huge; log1p keeps accuracy
    v <- seq_len(m) - 1L
    c(0, cumsum(log1p(v / x)))
  }
  
  # precompute prefixes up to n for the three S(., .) terms
  # each S(x, m) = m*log(x) + prefix[x][m]
  pref_a   <- prefix_log1p(alpha, n)
  pref_b   <- prefix_log1p(beta, n)
  pref_abn <- prefix_log1p(alpha + beta, n)
  
  # fast accessors using the prefixes
  S <- function(x, m, pref) {
    # m can be a vector; pref is length n+1 with pref[1]==0 for m=0
    m * log(x) + pref[m + 1L]
  }
  
  # core formula: lchoose + S(a,k) + S(b,n-k) - S(a+b, n)
  # lchoose is already stable for large n
  lpmf <- lchoose(n, k) +
    S(alpha, k, pref_a) +
    S(beta,  n - k, pref_b) -
    S(alpha + beta, n, pref_abn)
  
  if (debug) {
    # quick normalization check via log-sum-exp across full support
    kk <- 0:n
    l_all <- lchoose(n, kk) +
      S(alpha, kk, pref_a) +
      S(beta,  n - kk, pref_b) -
      S(alpha + beta, n, pref_abn)
    m0 <- max(l_all)
    lse <- m0 + log(sum(exp(l_all - m0)))
    message(sprintf("debug: log-sum-exp over k=0..%d = %.12g (should be 0)", n, lse))
  }
  lpmf
}


bb_log <- dbetabinom_log(k = stan_data$count, size = stan_data$total, alpha = samps$shape, beta = samps$shape)
bb_log <- extraDistr::dbbinom(stan_data$count, size = stan_data$total, alpha = samps$shape, beta = samps$shape, log = TRUE)
bb_log <- VGAM::dbetabinom.ab(stan_data$count, size = stan_data$total, 
                              shape1 = samps$shape, shape2 = samps$shape, log = TRUE)
bb_log <- sapply(samps$shape, function(si) 
  dbetabinom_log_stable(k = stan_data$count, size = stan_data$total, alpha = si, beta = si)
)

plot(bb_log, samps$log_beta_binomial)
abline(h = max(samps$log_beta_binomial))
abline(0,1)

plot(samps$log_concentration, bb_log, xlim =)
min(samps$log_concentration[bb_log = which.max(bb_log)])


dbetabinom_log_stable(k = 10, size = 21, alpha = 1E8, beta = 1E8)
dbetabinom_log_stable(k = 10, size = 21, alpha = 1E20, beta = 1E20)
