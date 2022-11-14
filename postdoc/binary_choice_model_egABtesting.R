#TODO add covariance? user effects w/ repeated measures?

#libaries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

k <- 5
n <- 2000
mus <- rnorm(k, 0, 2)
sds <- exp(rnorm(k, log(0.5), 0.5))
inds <- t(replicate(n, sample(1:k, 2)))
liab <- t(apply(inds, 1, function(x){
  c(rnorm(1, mus[x[1]], sds[x[1]]),
    rnorm(1, mus[x[2]], sds[x[2]]))
  }))

choice <- apply(liab, 1, which.max)
winner <- sapply(1:n, function(i) inds[i,choice[i]])
cbind(inds, winner)
inds_ord <- t(sapply(1:n, function(i) c(loser = inds[i,-choice[i]], winner = inds[i,choice[i]])))

#naive ranking
n_win <- n_try <- setNames(rep(0, k), 1:k)
n_win[names(table(winner))] <- table(winner)
n_try[names(table(inds))] <- table(inds)
plot(order(n_win / n_try), order(mus))

#bayesian ranking
d <- list(k = k,
          n = n,
          win = inds_ord[,"winner"],
          lose = inds_ord[,"loser"])

# STAN model
stan_program <- '
data {
    int<lower=1> k; // number of discrete choices
    int<lower=1> n; // number of comparisons
    int<lower=1,upper=k> win[n]; // index of winning choice
    int<lower=1,upper=k> lose[n]; // index of losing choice
}
parameters {
    vector[k] mu;
    vector[k] log_sd;
    vector[n] lose_liab;      
    vector<lower=0>[n] diff_liab;      
}
transformed parameters {
    vector[n] win_liab = lose_liab + diff_liab;      
}
model {
    // priors
    log_sd ~ std_normal(); // constrain log_sd to mean = 0 for identifiability
    mu ~ normal(0,10); // weakly regularizing prior for mu in case an item always wins
    
    // the liability model itself
    lose_liab ~ normal(mu[lose], exp(log_sd[lose]));
    win_liab ~ normal(mu[win], exp(log_sd[win]));
    
    // jacobian adjustment ot needed bc it is a linear transformation
    // target += log(fabs(1));

}
generated quantities {
}
'

# or use an ordered vector
stan_program <- '
data {
    int<lower=1> k; // number of discrete choices
    int<lower=1> n; // number of comparisons
    int<lower=1,upper=k> win[n]; // index of winning choice
    int<lower=1,upper=k> lose[n]; // index of losing choice
}
parameters {
    vector[k] mu;
    vector[k] log_sd;
    ordered[2] liab[n];
}
model {
    // priors
    log_sd ~ std_normal(); // constrain log_sd to mean = 0 for identifiability
    mu ~ normal(0,10); // weakly regularizing prior for mu in case an item always wins
    
    // the liability model itself
    liab[,1] ~ normal(mu[lose], exp(log_sd[lose]));
    liab[,2] ~ normal(mu[win], exp(log_sd[win]));
    
    // jacobian adjustment ot needed bc it is a linear transformation
    // target += log(fabs(1));

}
generated quantities {
}
'
# stan_program <- stan_program_direct
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = 0.1, max_treedepth = 20)
# out <- mod$variational(data = d)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  if(exclude != ""){
    excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
    return_inds <- setdiff(incl_inds, excl_inds)  
    return(samps[,return_inds])
  } else {
    return(samps[,incl_inds])
  }
}

pmean_mu <- apply(subset_samps("mu", samps = samps), 2, mean)
pmean_sd <- apply(exp(subset_samps("log_sd", samps = samps)), 2, mean)
plot(order(n_win / n_try), order(mus))
plot(order(pmean_mu), order(mus))

plot(mus, pmean_mu)
plot(sds, pmean_sd)

