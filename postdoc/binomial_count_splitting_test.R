#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#specify functions
logit <- function(p) log(p/(1-p))
invlogit <- function(x) {exp(x) / (1+exp(x))}

#simulate data
n <- 1000
n0 <- round(n / 2)
x <- rnorm(n)
xb <- rnorm(n)
xf <- x + xb / 2
xc <- x - xb / 2
pf <- invlogit(xf)
pc <- invlogit(xc)
tfs <- c(10, 5)
tf <- c(rep(tfs[1], n - n0), rep(tfs[2], n0))
tcs <- 2000
tc <- rep(tcs, n)
cf <- rbinom(n, tf, pf)
cc <- rbinom(n, tc, pc)

#fit bayesian model
d <- list(n = n,
          tf = tf,
          tc = tc,
          cf = cf,
          cc = cc
)

stan_program <- "
data {
  int<lower=1> n;
  int<lower=0> tf[n];
  int<lower=0> tc[n];
  int<lower=0> cf[n];
  int<lower=0> cc[n];
}
parameters {
  vector[n] x;
  vector[n] xb;
  real<lower=0> x_sd;
  real<lower=0> xb_sd;
}
transformed parameters {
  vector[n] xf = x + xb / 2;
  vector[n] xc = x - xb / 2;
}
model {
  //priors
  x ~ normal(0,x_sd);
  xb ~ normal(0,xb_sd);
  x_sd ~ std_normal();
  xb_sd ~ std_normal();

  //likelihood
  cf ~ binomial_logit(tf, xf);
  cc ~ binomial_logit(tc, xc);
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, adapt_delta = 0.85)
# summ <- out$summary()
# summ[order(summ$ess_bulk),]
# summ[order(summ$rhat, decreasing = T),]

samps <- data.frame(as_draws_df(out$draws()))

subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

#do some plotting
par(mar = c(4,4,4,4))
layout(rbind(c(1,2,3,4), c(5,6,6,7)))

plot(xf, apply(subset_samps("xf", "_", samps), 2, mean), pch = 19, col = c(rep(1, n - n0), rep(3, n0)), ylab = "posterior mean for xf"); abline(0, 1, lwd = 3, col = 2)
legend("topleft", pch = 19, col = c(1,3), legend = c(paste0("tf = ", tfs)), cex = 0.75)
legend("bottomright", lwd = 3, col = 2, legend = "1-to-1 line", cex = 0.75)

plot(xc, apply(subset_samps("xc", "_", samps), 2, mean), pch = 19, col = c(rep(1, n - n0), rep(3, n0)), ylab = "posterior mean for xc"); abline(0, 1, lwd = 3, col = 2)
legend("topleft", pch = 19, col = c(1,3), legend = c(paste0("tc = ", rep(tcs, 2))), cex = 0.75)
legend("bottomright", lwd = 3, col = 2, legend = "1-to-1 line", cex = 0.75)

plot(xb, apply(subset_samps("xb", "_", samps), 2, mean), pch = 19, col = c(rep(1, n - n0), rep(3, n0)), ylab = "posterior mean for xb"); abline(0, 1, lwd = 3, col = 2)
legend("topleft", pch = 19, col = c(1,3), legend = c(paste0("tf = ", tfs)), cex = 0.75)
legend("bottomright", lwd = 3, col = 2, legend = "1-to-1 line", cex = 0.75)

hist(samps$x_sd, probability = T, xlab = "posterior draws for x_sd"); abline(v = 1, lwd = 3, col = 2)
legend("topleft", lwd = 3, col = 2, legend = "true value", cex = 0.75)

hist(samps$xb_sd, probability = T, xlab = "posterior draws for xb_sd"); abline(v = 1, lwd = 3, col = 2)
legend("topleft", lwd = 3, col = 2, legend = "true value", cex = 0.75)

prop_greater_than_0 <- function(x) mean(x>0)
barplot(apply(subset_samps("xb", "_", samps = samps), 2, prop_greater_than_0), 
        col = c(rep(1, n - n0), rep(3, n0)), border = c(rep(1, n - n0), rep(3, n0)),
        ylab = "proportion of posterior mass above 0 (xb)", width = 0.1,
        xlab = "parameter index for xb")

hist(apply(subset_samps("xb", "_", samps = samps), 2, prop_greater_than_0)[(n-n0+1):n],
     xlab = "proportion of posterior mass above 0 (xb)", breaks = 100, col = 3, border = 3, main = "")

