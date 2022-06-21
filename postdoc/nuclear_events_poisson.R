#packages
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(bayesplot)

#write data
count <- c(1, 5, 1, 2, 1, 0, 0)


d <- list(count = count, 
          decade = seq_along(count), 
          n = length(count)
          )

#basic multilevel model
stan_program <- '
data {
    int<lower=1> n;
    int<lower=0> count[n];
    int<lower=1,upper=n> decade[n];
}
parameters {
    vector[n] raw_lograte;
    real lograte_mean;
    real<lower=0> lograte_sd;
}
transformed parameters {
    vector[n] lograte = raw_lograte * lograte_sd + lograte_mean;
}
model {
    raw_lograte ~ std_normal();
    lograte_sd ~ std_normal();
    lograte_mean ~ std_normal();
    count ~ poisson_log(lograte);
}
generated quantities {

}
'

#autoregressive model with trend
stan_program <- '
data {
    int<lower=1> n;
    int<lower=0> count[n];
    int<lower=1,upper=n> decade[n];
}
transformed data {
    vector[n] std_decade = to_vector(decade) - mean(to_vector(decade));
}
parameters {
    vector[n] raw_lograte_error;
    real lograte_intercept;
    real<lower=0> lograte_sd;
    real lograte_slope;
    real<lower=-1, upper=1> rho;
}
transformed parameters {
    vector[n] lograte_error;
    lograte_error[1] = raw_lograte_error[1] * 2;
    for(t in 2:n){
      lograte_error[t] = raw_lograte_error[t] * lograte_sd + lograte_error[t-1] * rho;
    }
}
model {
    raw_lograte_error ~ std_normal();
    rho ~ std_normal();
    lograte_sd ~ std_normal();
    lograte_intercept ~ std_normal();
    lograte_slope ~ std_normal();
    
    vector[n] lograte = lograte_intercept + lograte_slope * std_decade + lograte_error;
    count ~ poisson_log(lograte);
}
generated quantities {
    real lograte_next = normal_rng(lograte_intercept + lograte_slope * (std_decade[n] + 1) + lograte_error[n] * rho, lograte_sd);
    int count_next = poisson_log_rng(lograte_next);
}
'


if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)


#fit model
base = "nuclear_poisson"
write_stan_file(stan_program, dir = "~/Desktop/", basename = base)
write_stan_json(d, paste0("~/Desktop/", base,".json"))
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                  adapt_delta = 0.9, refresh = 100, init = 0.1, max_treedepth = 20, thin = 1, step_size = 0.1)

#check mcmc diagnostics
summ <- out$summary()
head(summ[order(summ$ess_bulk),])
head(summ[order(summ$rhat, decreasing = T),])

#query posterior & posterior predictive
samps <- data.frame(as_draws_df(out$draws()))
mean(samps$lograte_slope < 0)
table(samps$count_next) / length(samps$count_next)
