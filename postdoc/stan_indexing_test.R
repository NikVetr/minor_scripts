#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#specify functions


#fit bayesian model
n <- 100
d <- list(nA = 5,
          nB = 6,
          nC = 3,
          n = n,
          A = sample(1:5, n, T),
          B = sample(1:6, n, T),
          C = sample(1:3, n, T)
)

stan_program <- "
data {
  int<lower=1> nA;
  int<lower=1> nB;
  int<lower=1> nC;
  int<lower=1> n;
  
  int<lower=1, upper=nA> A[n];
  int<lower=1, upper=nB> B[n];
  int<lower=1, upper=nC> C[n];
}
parameters {
  real x[nA,nB,nC];
}
transformed parameters {
  vector[n] y = x[A,B,C];
}
model {
  x ~ normal(0,1);
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
