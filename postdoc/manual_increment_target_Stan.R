
#### ok, what about a Bayesian DE model ####
# dsub <- d[d$time %in% c(0,2),]
dsub <- d
dat <- list(n = 1E4,
            count = rpois(1E4, 1500)
)

stan_program_1 <- '
data {
    int<lower=1> n;
    int<lower=0> count[n];
}
parameters {
    real<lower=0> rate;
}
model {
    rate ~ normal(0, 1000);
    count ~ poisson(rate);
}
'

stan_program_2 <- '
data {
    int<lower=1> n;
    vector<lower=0>[n] count;
}
parameters {
    real<lower=0> rate;
}
model {
    rate ~ normal(0, 1000);
    //count ~ poisson(rate);
    //target += count * log(rate) - rate - log(tgamma(count + 1.0));
    target += count * log(rate) - rate - (count + 0.5) .* log(count) - count + log(2.0 * pi()) / 2;
}
'

if(!exists("curr_stan_program_1") || stan_program_1 != curr_stan_program_1){
  curr_stan_program_1 <- stan_program_1
  f_1 <- write_stan_file(stan_program_1)
}
mod_1 <- cmdstan_model(f_1)

if(!exists("curr_stan_program_2") || stan_program_2 != curr_stan_program_2){
  curr_stan_program_2 <- stan_program_2
  f_2 <- write_stan_file(stan_program_2)
}
mod_2 <- cmdstan_model(f_2)


out_1 <- mod_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 1000, init = 0.1, max_treedepth = 15, thin = 1)
out_2 <- mod_2$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 1000, init = 2000, max_treedepth = 15, thin = 1)

out_1$summary()
out_2$summary()

samps_1 <- data.frame(as_draws_df(out_1$draws()))
samps_2 <- data.frame(as_draws_df(out_2$draws()))

plot(quantile(samps_1$rate, probs = 1:99/100), quantile(samps_2$rate[1:3000], probs = 1:99/100))
abline(0,1,col=2,lty=2)
