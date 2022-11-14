#libraries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(mvtnorm)

#functions
invlogit <- function(x) exp(x)/(1+exp(x))
prop_greater_than_0 <- function(x) mean(x>0)
subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

#### basic univariate test ####
n <- 20
p <- 500
mu_sd <- 0
mu <- rnorm(p, sd = mu_sd)
x <- sapply(mu, function(mu_i) rnorm(n, mu_i))

d <- list(n = n,
          p = p,
          i = rep(1:p, each = n),
          x = c(x))

base = "basic_multilevel_model"
stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> p;
    int<lower=1, upper=p> i[n*p];
    real x[n*p];
}
parameters {
    vector[p] raw_mu;
    real<lower=0> mu_sd;
    real mu_mu;
    real<lower=0> sd;
}
transformed parameters {
  vector[p] mu = raw_mu * mu_sd + mu_mu;
}
model {
    //priors
    mu_sd ~ normal(0, 0.1);
    sd ~ std_normal();
    mu_mu ~ normal(0,5);
    raw_mu ~ std_normal();
    
    //uncensored obs likelihood
    x ~ normal(mu[i], sd);
}
'


if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
write_stan_file(stan_program, dir = "~/Desktop/", basename = paste0(base))
write_stan_json(d, paste0("~/Desktop/", paste0(base, ".json")))
fit_model <- T
if(fit_model){
  out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 100, init = 0.1, max_treedepth = 15, thin = 2)
  summ <- out$summary()
  print(summ[order(summ$ess_bulk),])
  print(summ[order(summ$rhat, decreasing = T),])
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
}


samps <- data.frame(as_draws_df(out$draws()))

mu_samps <- subset_samps("mu", c("mu_mu", "mu_sd", "raw_mu"), samps)

#####

par(mfrow = c(1,2))
q.mu_i <- sapply(1:p, function(i) mean(mu_samps[,i] > mu[i]))
hist(q.mu_i, breaks = 0:20/20, probability = T, xlab = latex2exp::TeX("quantile of true $\\mu_i$ in marginal posterior"), main = latex2exp::TeX("$\\mu_i$"))
CI_range <- 0:100/100
CI_coverage.mu_i <- sapply(CI_range, function(CIr) sum((sapply(q.mu_i, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.mu_i))
plot(CI_range, CI_coverage.mu_i, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value", main = latex2exp::TeX("$\\mu_i$"))
abline(0,1,lty=2,lwd=2,col=2)
legend("topleft", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")


#####


#####
# n <- 20
# p <- 500
# x_mu_1 <- rnorm(p)
# x_mu_2 <- rnorm(p)
# x_mu_diff <- x_mu_1 - x_mu_2
# x_1 <- sapply(x_mu_1, function(x_mu) rnorm(n, x_mu))
# x_2 <- sapply(x_mu_2, function(x_mu) rnorm(n, x_mu))
# 
# d <- list(n = n,
#           p = p,
#           i = rep(1:p, each = n),
#           x_1 = c(x_1),
#           x_2 = c(x_2))
# 
# 
# #load libraries
# library(MASS)
# library(mvtnorm)
# library(cmdstanr)
# library(posterior)
# 
# base = "basic_multilevel_difference_1"
# stan_program <- '
# data {
#     int<lower=1> n;
#     int<lower=1> p;
#     int<lower=1, upper=p> i[n*p];
#     real x_1[n*p];
#     real x_2[n*p];
# }
# parameters {
#     vector[p] x_mu_1;
#     vector[p] x_mu_2;
# }
# transformed parameters {
# }
# model {
#     //priors
#     x_mu_1 ~ std_normal();
#     x_mu_2 ~ std_normal();
#     
#     //uncensored obs likelihood
#     x_1 ~ normal(x_mu_1[i], 1);
#     x_2 ~ normal(x_mu_2[i], 1);
# }
# generated quantities {
#     vector[p] x_mu_diff = x_mu_1 - x_mu_2;
# }
# '
# 
# 
# if(!exists("curr_stan_program") || stan_program != curr_stan_program){
#   curr_stan_program <- stan_program
#   f <- write_stan_file(stan_program)
# }
# mod <- cmdstan_model(f)
# 
# #fit model
# write_stan_file(stan_program, dir = "~/Desktop/", basename = paste0(base))
# write_stan_json(d, paste0("~/Desktop/", paste0(base, ".json")))
# fit_model <- F
# if(fit_model){
#   out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
#                     adapt_delta = 0.85, refresh = 100, init = 0.1, max_treedepth = 15, thin = 2)
#   summ <- out$summary()
#   print(summ[order(summ$ess_bulk),])
#   print(summ[order(summ$rhat, decreasing = T),])
#   save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
# } else {
#   load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
# }
# 
# #####