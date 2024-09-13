library(cmdstanr)
library(posterior)

#### attempt 1 ####

stan_model <- "
parameters {
  real<lower=0> x;
  real<lower=0> y;
}
transformed parameters {
  real z1 = log(y)^2 + log(x)^2;
  //real z2 = (log(y) + log(x)) / sqrt(2);
}

model {
  // centered distributions on latent variables
  z1 ~ chi_square(2);
  //z2 ~ std_normal();
  x ~ lognormal(0,1);

  // Jacobian adjustment?
  target += log(abs(2 * (1/y * log(y))));
}
"

# D(expression(log(x)^2 + log(y)^2), "x")
# D(expression((log(x) + log(y)) / sqrt(2)), "y")
# D(expression(log(y)^2), "y")
# D(expression(x + y), "y")

mod <- cmdstan_model(write_stan_file(stan_model))
fit <- mod$sample(chains = 4, iter_sampling = 1E4, iter_warmup = 1E4,
                  adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 2, init = 0.1)

#test continuation
fit_metadata <- fit$metadata()
step_size <- fit_metadata$step_size_adaptation
last_draws <- fit$draws(variables = NULL) # Extract all draws
init_values <- lapply(1:dim(last_draws)[2], function(chain) {
  as.list(last_draws[1, chain, , drop = TRUE])
})
new_fit <- mod$sample(chains = 4, 
                      iter_sampling = 1E4,  # New number of iterations
                      iter_warmup = 0,      # Skip warmup as you've already done it
                      step_size = step_size, 
                      init = init_values, 
                      adapt_delta = 0.9, 
                      parallel_chains = 4, 
                      refresh = 100, 
                      max_treedepth = 10, 
                      thin = 2,
                      adapt_engaged = FALSE)

summ <- fit$summary()
print(summ[order(summ$rhat, decreasing = T),])
samps <- data.frame(as_draws_df(fit$draws()))

varns <- c("x", "y", "z1")
for(varn in varns){
  qs <- 1:999/1000
  fsamps <- samps[,varn]
  if(varn == "z1"){
    dirqs <- qchisq(qs, df = 2)  
  } else {
    dirqs <- qlnorm(qs)  
  }
  
  plot(dirqs, quantile(fsamps, probs = qs), type = "l", 
       xlab = paste0("true ", varn, " quantiles"), 
       ylab = paste0("sampled ", varn, " quantiles"))
  abline(0, 1, col = "red")
  
}

#### test 2 ####

stan_model <- "
parameters {
  real<lower=0> x;
}
transformed parameters {
  real z = log(x);
}

model {
  // centered distributions on latent variables
  z ~ normal(0, 1);
  
  // Jacobian adjustment?
  target += -log(x);
}
"

mod <- cmdstan_model(write_stan_file(stan_model))
fit <- mod$sample(chains = 4, iter_sampling = 1E4, iter_warmup = 1E4,
                  adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 2, init = 0.1)
summ <- fit$summary()
print(summ[order(summ$rhat, decreasing = T),])
samps <- data.frame(as_draws_df(fit$draws()))

varn <- "x"
# dir_samp <- rchisq(n = 1E4, df = 2)
qs <- 1:999/1000
fsamps <- samps[,varn]
dirqs <- qlnorm(qs)  

plot(dirqs, quantile(fsamps, probs = qs), type = "l")
abline(0, 1, col = "red")

