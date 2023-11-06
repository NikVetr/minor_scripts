n <- 1E4
s <- c(0.5,2)
p <- rbeta(n, s[1], s[2])
k <- 2
x <- rbinom(n, k, p)
y <- k - x
ep <- (x + 0.01) / (k + 0.02)

tep <- table(ep) / n
tep
bdens <- dbeta(as.numeric(names(tep)), s[1], s[2])
plot(c(tep), bdens / sum(bdens))

# fitdistrplus::fitdist(ep, "beta", method = "mle")
library(cmdstanr)

#fit bayesian model
d <- list(n = n,
         k = k,
         x = x)

stan_program <- "
data {
  int<lower=0> n;
  int<lower=0> k;
  int<lower=0,upper=k> x[n];
}
parameters {
  //vector<lower=0,upper=1>[n] p;
  vector<lower=0>[2] s;
}
model {
  //x ~ binomial(k, p);
  //p ~ beta(s[1], s[2]);
  x ~ beta_binomial(k, s[1], s[2]);
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, adapt_delta = 0.95)
samps <- data.frame(posterior::as_draws_df(out$draws()))
par(mfrow = c(2,1))
hist(samps$s.1.); abline(v = s[1], col = 2, lwd = 2)
hist(samps$s.2.); abline(v = s[2], col = 2, lwd = 2)

#why doesn't Stan auto-parse unique sites?
#maybe look for lp statements that are identical during warmup?
