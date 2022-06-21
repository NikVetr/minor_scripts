#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#specify functions
logit <- function(p) log(p/(1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

#### example 1 ####
n <- 1E3
x <- rnorm(n)
pr_x_m <- invlogit(x - 0.5)
xmi <- which(as.logical(rbinom(n, 1, pr_x_m)))
xm <- x; xm[xmi] <- NA
xc <- xm[!is.na(xm)]

breaks <- seq(min(x), max(x), length.out = 20)
hist(x, col = 2, breaks = breaks)
hist(xc, add = T, col = 3, breaks = breaks)
legend("topleft", legend = c("observed", "missing"), col = c(2,3), pch = 15, pt.cex = 2)
mean(x)
mean(xc)

#fit bayesian model
n <- 100
d <- list(
  nc = length(xc),
  nm = length(xmi),
  xc = xc
)

stan_program <- "
data {
  int<lower=1> nc;
  int<lower=0> nm;
  real xc[nc];
}
parameters {
  real mu;
  real<lower=0> sigma;
  real xm[nm];
  real xm_pr_b;
}
model {
  //missingness model
  target += bernoulli_logit_lpmf(0 | to_vector(xc) + rep_vector(xm_pr_b, nc));
  target += bernoulli_logit_lpmf(1 | to_vector(xm) + rep_vector(xm_pr_b, nm));

  //likelihood
  xc ~ normal(mu,sigma);
  xm ~ normal(mu,sigma);
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
mean(x)
mean(xc)
mean(x[xmi])
mean(samps$mu)
mean(samps$sigma)
mean(samps$xm_pr_b)


#####
n <- 1E3
x <- rnorm(n)
b1 <- 0.5
y <- x * b1 + rnorm(n)
a <- 3
b2 <- 0.75
b3 <- 1.5
z <- a + x * b2 + y * b3 + rnorm(n)

d <- data.frame(x = x, y = y, z = z)
lm(z ~ x + y)

#missingness model
pr_x_m <- abs(x + y)
pr_x_m <- invlogit((pr_x_m - quantile(pr_x_m, 0.9)) / sd(pr_x_m))

pr_y_m <- abs(x + z)
pr_y_m <- invlogit((pr_y_m - quantile(pr_y_m, 0.9)) / sd(pr_y_m))

xmi <- which(as.logical(rbinom(n, 1, pr_x_m)))
ymi <- which(as.logical(rbinom(n, 1, pr_y_m)))

xm <- x; xm[xmi] <- NA
ym <- y; ym[ymi] <- NA

dm <- data.frame(x = xm, y = ym, z = z)
dc <- dm[complete.cases(dm),]

if(n <= 1E3){
  pairs(dc)
  pairs(d)
}


lm(z ~ x + y, data = d)
lm(z ~ x + y, data = dc)
