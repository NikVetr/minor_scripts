#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (exp(x) + 1)
softmax <- function(x) exp(x) / sum(exp(x))
subset_samps <- function(include = NULL, exclude = NULL, samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

#params
n <- 1000
p <- 5
n_groups <- sample(2:10, p, T)
group_freqs <- lapply(lapply(n_groups, rnorm), softmax)
indivs <- t(sapply(1:n, function(i) sapply(1:p, function(j) sample(x = 1:n_groups[j], 1, prob = group_freqs[[j]]))))

sds <- abs(rnorm(p) * 2)
effects <- lapply(lapply(n_groups, rnorm), function(x) (x - mean(x)) / sd(x))
mean_lodds <- 0

success_probs <- invlogit(apply(indivs, 1, function(x) mean_lodds + sum(sapply(1:p, function(i) effects[[i]][x[i]]))))
n_rounds <- 50
pass <- replicate(n_rounds, rbinom(n, 1, success_probs))
ko <- apply(pass, 1, function(x) min(which(c(x, 0) == 0)))
rev(cumsum(rev(table(ko))))
hist(ko) 
#those who make it to the highest round have the features likely to let them keep going -- need to adjust strength of baseline filter
#whereas those who are penalized get knocked out super early
#need to just do an iterated bernoulli model

#fit bayesian model
d <- list(n = n,
          p = p,
          n_groups = n_groups,
          n_rounds = n_rounds,
          indivs = indivs,
          ko = ko,
          n_failure = sum(ko <= n_rounds),
          inds_failure = which(ko <= n_rounds),
          inds_success = which(ko > n_rounds)
)

#fundamentally a right-censored negative binomial model
stan_program <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> n_groups[p];
  int<lower=1> n_rounds;
  int indivs[n, p];
  int<lower=1> ko[n];
  int<lower=0> n_failure;
  int inds_failure[n_failure];
  int inds_success[n - n_failure];
}
transformed data{
  int<lower=1> b_size = sum(n_groups);
  int<lower=1> b_cut[p];
  b_cut[1] = 1;
  for(i in 2:p){
    b_cut[i] = b_cut[i-1] + n_groups[i-1];
  }
}
parameters {
  real a;
  vector[b_size] b;
  vector<lower=0>[p] sigma;
}
transformed parameters{
  vector[n] b_indiv;
  for(i in 1:n){
    real b_temp = a;
    for(j in 1:p){
      b_temp = b_temp + b[b_cut[j] + indivs[i,j] - 1];
    }
    b_indiv[i] = b_temp;
  }
}
model {
  //priors  
  a ~ normal(0,2);
  sigma ~ normal(0,1);
  for(i in 1:p){
    segment(b, b_cut[i], n_groups[i]) ~ normal(0,sigma[i]);  
  }

  //likelihood
  target += neg_binomial_lpmf(ko[inds_failure] | 1 , inv_logit(b_indiv[inds_failure]));
  target += log(1 - exp(neg_binomial_lcdf(n_rounds + 1 | 1, inv_logit(b_indiv[inds_success]))));
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
# summ <- out$summary()
# summ[order(summ$ess_bulk),]
# summ[order(summ$rhat, decreasing = T),]

samps <- data.frame(as_draws_df(out$draws()))

b_post <- subset_samps("b\\.", exclude = "b_indiv", samps = samps)

hist(samps$a); abline(v = -mean_lodds, lwd = 2, col = 2)
plot(unlist(effects), apply(b_post, 2, mean))
abline(0,-1, lwd = 2, col = 2)


#### alternate parameterization w/ binomials only ####
n <- 500
p <- 3
n_groups <- sample(2:10, p, T)
group_freqs <- lapply(lapply(n_groups, rnorm), softmax)
indivs <- t(sapply(1:n, function(i) sapply(1:p, function(j) sample(x = 1:n_groups[j], 1, prob = group_freqs[[j]]))))

sds <- abs(rnorm(p) * 2)
effects <- lapply(lapply(n_groups, rnorm), function(x) (x - mean(x)) / sd(x))
n_rounds <- 5
base_lodds <- 0.5 - 0:(n_rounds-1) * 1

success_prob_deviation <- apply(indivs, 1, function(x) sum(sapply(1:p, function(i) effects[[i]][x[i]])))
pass <- matrix(0, nrow = n, ncol = n_rounds)
pass[,1] <- rbinom(n, 1, invlogit(base_lodds[1] + success_prob_deviation))
for(i in 2:n_rounds){
  prev_passers <- pass[,i-1] == 1
  pass[prev_passers,i] <- rbinom(sum(prev_passers), 1, invlogit(base_lodds[i] + success_prob_deviation[prev_passers]))
  # mean(invlogit(success_prob_deviation[which(pass[,i-1] == 1)]))
}
ko <- apply(pass, 1, function(x) min(which(c(x, 0) == 0)))
n_applicants_per_round <- rev(cumsum(rev(table(ko))))
par(mfrow = c(2,1))
plot(n_applicants_per_round); lines(n_applicants_per_round)
hist(ko)

#k-way interactions?
sapply(1:p, function(i){sum(apply(combn(p, i), 2, function(x) prod(n_groups[x])))})
interactions <- apply(combn(p, 2), 2, function(x) apply(indivs[,x], 1, paste0, collapse = ":"))
int_tab <- apply(interactions, 2, table)
length(unlist(int_tab))

#fit bayesian model
d <- list(n = n,
          p = p,
          n_groups = n_groups,
          n_rounds = n_rounds,
          indivs = indivs,
          ko = ko
)

#fundamentally an iterated bernoulli model
stan_program <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> n_groups[p];
  int<lower=1> n_rounds;
  int indivs[n, p];
  int<lower=1> ko[n];
}
transformed data{
  int<lower=1> b_size = sum(n_groups);
  int<lower=1> b_cut[p];
  b_cut[1] = 1;
  for(i in 2:p){
    b_cut[i] = b_cut[i-1] + n_groups[i-1];
  }
}
parameters {
  real a[n_rounds];
  vector[b_size] b;
  vector<lower=0>[p] sigma;
}
model {
  //priors  
  a ~ normal(0,5);
  sigma ~ normal(0,1);
  for(i in 1:p){
    segment(b, b_cut[i], n_groups[i]) ~ normal(0,sigma[i]);  
  }

  //model for group bias
  vector[n] b_indiv;
  for(i in 1:n){
    real b_temp = 0;
    for(j in 1:p){
      b_temp = b_temp + b[b_cut[j] + indivs[i,j] - 1];
    }
    b_indiv[i] = b_temp;
  }
  
  //likelihood
  for(i in 1:n){
    if(ko[i] <= n_rounds){
      target += binomial_logit_lpmf(0 | 1 , a[ko[i]] + b_indiv[i]);
    }
    if(ko[i] > 1){
      for(j in 1:(ko[i]-1)){
        target += binomial_logit_lpmf(1 | 1 , a[j] + b_indiv[i]);
      }
    }
  }
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E2, iter_warmup = 2E2, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10)
summ <- out$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]

samps <- data.frame(as_draws_df(out$draws()))
a_post <- subset_samps("a\\.", exclude = "sigma", samps = samps)
a_90QI <- apply(a_post, 2, quantile, prob = c(0.05, 0.95))
b_post <- subset_samps("b\\.", exclude = "b_indiv", samps = samps)
b_90QI <- apply(b_post, 2, quantile, prob = c(0.05, 0.95))
true_b <- unlist(effects)

plot(base_lodds, apply(a_post, 2, mean), ylim = range(a_90QI)); abline(0,1, lwd = 2, col = 2)
for(i in seq_along(base_lodds)){segments(x0 = base_lodds[i], x1 = base_lodds[i], y0 = a_90QI[1, i], y1 = a_90QI[2, i])}

plot(true_b, apply(b_post, 2, mean), ylim = range(b_90QI)); abline(0,1, lwd = 2, col = 2)
for(i in seq_along(true_b)){segments(x0 = true_b[i], x1 = true_b[i], y0 = b_90QI[1, i], y1 = b_90QI[2, i])}

#### trying to figure out relation between pass rate and idiosyncratic effects ####
# logit <- function(p) log(p / (1-p))
# invlogit <- function(x) exp(x) / (exp(x) + 1)
# z <- (rnorm(50, mean = 2, sd = 2))
# x <- -200:200/10
# y <- sapply(x, function(i) (mean(invlogit(z + i)))) - mean(invlogit(z))
# plot(x, y, type = "l", ylab = "q - p")
# for(i in 1:20){
#   z <- (rnorm(50, mean = 2, sd = 2))
#   y <- sapply(x, function(i) (mean(invlogit(z + i)))) - mean(invlogit(z))
#   lines(x, y)
# }
# 
# 
# logistic <- function(L, k, x, x0, u){
#   L / (1 + exp(-k*(x-x0))) + u
# }
# 
# sslogistic <- function(par, data){
#   parl <- as.list(par)
#   preds <- logistic(parl$L, parl$k, data$x, parl$x0, parl$u)
#   return(sum((data$y-preds)^2))
# }
# 
# out <- optimx::optimx(par = c(L = 0, k = 0, x0 = 0, u = 0), 
#       fn = sslogistic,
#       data = list(x = x, y = y),
#       method = "BFGS"
# )
# 
# out
# plot(x, y, type = "l", lwd = 3, col = adjustcolor(1,0.5),)
# lines(x, logistic(L = out$L, k = out$k, x0 = out$x0, x = x, u = out$u), col = adjustcolor(2,0.5), lty = 2, lwd = 3)

