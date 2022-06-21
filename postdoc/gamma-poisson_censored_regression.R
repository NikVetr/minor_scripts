library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(mvtnorm)

n_indiv <- 100
n_prot <- 1000

t0_mean <- log(6E5)
t0_mean <- 2

t0_prot_sd <- 1
t0_prot <- rnorm(n_prot, 0, t0_prot_sd)

t0_indiv_sd <- 0.5
t0_indiv <- rnorm(n_indiv, 0, t0_indiv_sd)

t01_mean <- 1
t01_prot_sd <- 2
t01_prot <- rnorm(n_prot, 0, t01_prot_sd)

t01_indiv_sd <- 0.5
t01_indiv <- rnorm(n_indiv, 0, t01_indiv_sd)

gamma_dispersion <- 0.01

gamma_alpha <- 1 / gamma_dispersion
gamma_beta_t0 <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) gamma_alpha / exp(t0_mean + t0_prot[prot_i] + t0_indiv[indiv_i])))
t0_rates <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) 
  rgamma(1, shape = gamma_alpha, rate = gamma_beta_t0[[prot_i]][indiv_i])))

gamma_beta_t1 <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) gamma_alpha / 
                                                            exp(t0_mean + t0_prot[prot_i] + t0_indiv[indiv_i] + t01_mean + t01_prot[prot_i] + t01_indiv[indiv_i])))
t1_rates <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) 
  rgamma(1, shape = gamma_alpha, rate = gamma_beta_t1[[prot_i]][indiv_i])))

t0_counts <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) rpois(1, lambda = t0_rates[[prot_i]][indiv_i])))
t1_counts <- lapply(1:n_prot, function(prot_i) sapply(1:n_indiv, function(indiv_i) rpois(1, lambda = t1_rates[[prot_i]][indiv_i])))

data <- as.data.frame(rbind(cbind(count = unlist(t0_counts), timepoint = 1, 
                                  indiv = rep(1:n_indiv, n_prot), 
                                  prot = unlist(lapply(1:n_prot, function(pi) rep(pi, n_indiv)))),
                            cbind(count = unlist(t1_counts), timepoint = 2, 
                                  indiv = rep(1:n_indiv, n_prot), 
                                  prot = unlist(lapply(1:n_prot, function(pi) rep(pi, n_indiv))))
                            ))

censor_threshold <- round(quantile(data$count, 0.3))
# censor_threshold <- 0

par(mfrow = c(2,1))
hist(log(data$count), breaks = quantile(log(data$count), 0:10/10)); abline(v = log(censor_threshold), col = "red", lwd = 2)
plot(log(data$count), col = as.integer(log(data$count) < log(censor_threshold)) + 1); abline(h = log(censor_threshold), col = "red", lwd = 2)

d_comp <- data[data$count > censor_threshold,]
d_cens <- data[data$count <= censor_threshold,]
unique(d_cens)

d <- list(n_indiv = n_indiv,
          n_prot = n_prot,
          n_comp = nrow(d_comp),
          n_cens = nrow(d_cens),
          n_timepoint = 2,
          count = d_comp$count,
          timepoint_comp = d_comp$timepoint,
          indiv_comp = d_comp$indiv,
          prot_comp = d_comp$prot,
          timepoint_cens = d_cens$timepoint,
          indiv_cens = d_cens$indiv,
          prot_cens = d_cens$prot,
          censor_threshold = censor_threshold
)


d <- list(n_indiv = n_indiv,
          n_prot = n_prot,
          
          n_comp_0 = sum(d_comp$timepoint == 1),
          n_cens_0 = sum(d_cens$timepoint == 1),
          count_0 = d_comp$count[d_comp$timepoint == 1],
          indiv_comp_0 = d_comp$indiv[d_comp$timepoint == 1],
          prot_comp_0 = d_comp$prot[d_comp$timepoint == 1],
          indiv_cens_0 = d_cens$indiv[d_cens$timepoint == 1],
          prot_cens_0 = d_cens$prot[d_cens$timepoint == 1],
          
          n_comp_1 = sum(d_comp$timepoint == 2),
          n_cens_1 = sum(d_cens$timepoint == 2),
          count_1 = d_comp$count[d_comp$timepoint == 2],
          indiv_comp_1 = d_comp$indiv[d_comp$timepoint == 2],
          prot_comp_1 = d_comp$prot[d_comp$timepoint == 2],
          indiv_cens_1 = d_cens$indiv[d_cens$timepoint == 2],
          prot_cens_1 = d_cens$prot[d_cens$timepoint == 2],
          
          censor_threshold = censor_threshold
)

# data {
#   int<lower=1> n_indiv;
#   int<lower=1> n_prot;
#   int<lower=1> n_comp;
#   int<lower=1> n_cens;
#   int<lower=1> n_timepoint;
#   
#   int<lower=0> count[n_comp];
#   
#   int<lower=1,upper=n_timepoint> timepoint_comp[n_comp];
#   int<lower=1,upper=n_indiv> indiv_comp[n_comp];
#   int<lower=1,upper=n_prot> prot_comp[n_comp];
#   
#   int<lower=1,upper=n_timepoint> timepoint_cens[n_cens];
#   int<lower=1,upper=n_indiv> indiv_cens[n_cens];
#   int<lower=1,upper=n_prot> prot_cens[n_cens];
#   
#   int<lower=0> censor_threshold;
# }
# //vector[n_timepoint] mu_lr;
# //matrix[n_prot, n_timepoint] raw_prot_lr;
# //matrix[n_indiv, n_timepoint] raw_indiv_lr;
# //vector<lower=0>[n_timepoint] prot_lr_sd;
# //vector<lower=0>[n_timepoint] indiv_lr_sd;

base = "basic_2timepoint_censored_poisson"
# STAN model
stan_program <- '
data {
    int<lower=1> n_indiv;
    int<lower=1> n_prot;
    int<lower=1> n_comp_0;
    int<lower=1> n_comp_1;
    int<lower=0> n_cens_0;
    int<lower=0> n_cens_1;

    int<lower=0> count_0[n_comp_0];
    int<lower=0> count_1[n_comp_1];
    
    int<lower=1,upper=n_indiv> indiv_comp_0[n_comp_0];
    int<lower=1,upper=n_prot> prot_comp_0[n_comp_0];
    int<lower=1,upper=n_indiv> indiv_comp_1[n_comp_1];
    int<lower=1,upper=n_prot> prot_comp_1[n_comp_1];
    
    int<lower=1,upper=n_indiv> indiv_cens_0[n_cens_0];
    int<lower=1,upper=n_prot> prot_cens_0[n_cens_0];
    int<lower=1,upper=n_indiv> indiv_cens_1[n_cens_1];
    int<lower=1,upper=n_prot> prot_cens_1[n_cens_1];
    
    int<lower=0> censor_threshold;
}
parameters {
    real t0_mean;
    real t01_mean;
    
    vector[n_prot] raw_t0_prot;
    real<lower=0> t0_prot_sd;
    vector[n_prot] raw_t01_prot;
    real<lower=0> t01_prot_sd;
    
    vector[n_indiv] raw_t0_indiv;
    real<lower=0> t0_indiv_sd;
    vector[n_indiv] raw_t01_indiv;
    real<lower=0> t01_indiv_sd;
    
}
transformed parameters {
    vector[n_prot] t0_prot = raw_t0_prot * t0_prot_sd;
    vector[n_prot] t01_prot = raw_t01_prot * t01_prot_sd;
    vector[n_prot] t1_prot = t0_prot + t01_prot;
    
    vector[n_indiv] t0_indiv = raw_t0_indiv * t0_indiv_sd;
    vector[n_indiv] t01_indiv = raw_t01_indiv * t01_indiv_sd;
    vector[n_indiv] t1_indiv = t0_indiv + t01_indiv;
    
    real t1_mean = t0_mean + t01_mean;
}
model {
    //priors
    t0_mean ~ normal(0,10);
    t01_mean ~ normal(0,5);
    
    raw_t0_prot ~ std_normal();
    t0_prot_sd ~ std_normal();
    raw_t01_prot ~ std_normal();
    t01_prot_sd ~ std_normal();
    
    raw_t0_indiv ~ std_normal();
    t0_indiv_sd ~ std_normal();
    raw_t01_indiv ~ std_normal();
    t01_indiv_sd ~ std_normal();
    
    //poisson model for noncensored data
    count_0 ~ poisson(exp(t0_mean + t0_prot[prot_comp_0] + t0_indiv[indiv_comp_0]));
    count_1 ~ poisson(exp(t1_mean + t1_prot[prot_comp_1] + t1_indiv[indiv_comp_1]));
    
    //censored obs
    target += poisson_lcdf(censor_threshold | exp(t0_mean + t0_prot[prot_cens_0] + t0_indiv[indiv_cens_0]));
    target += poisson_lcdf(censor_threshold | exp(t1_mean + t1_prot[prot_cens_1] + t1_indiv[indiv_cens_1]));
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
  out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 50, init = 0.1, max_treedepth = 20, thin = 2)
  summ <- out$summary()
  print(summ[order(summ$ess_bulk),])
  print(summ[order(summ$rhat, decreasing = T),])
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
}

samps <- data.frame(as_draws_df(out$draws()))


plot(samps$lp__, type = "l")



hist(samps$t0_mean)
abline(v = t0_mean, col = 2, lwd = 2)

hist(samps$t0_prot_sd)
abline(v = t0_prot_sd, col = 2, lwd = 2)

hist(samps$t0_indiv_sd)
abline(v = t0_indiv_sd, col = 2, lwd = 2)

hist(samps$t01_mean)
abline(v = t01_mean, col = 2, lwd = 2)

hist(samps$t01_prot_sd)
abline(v = t01_prot_sd, col = 2, lwd = 2)

hist(samps$t01_indiv_sd)
abline(v = t01_indiv_sd, col = 2, lwd = 2)


subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

plot(t0_prot, apply(subset_samps("t0_prot", c("sd", "raw"), samps), 2, mean)); abline(0,1,col=2)
plot(t01_prot, apply(subset_samps("t01_prot", c("sd", "raw"), samps), 2, mean)); abline(0,1,col=2)

plot(t0_indiv, apply(subset_samps("t0_indiv", c("sd", "raw"), samps), 2, mean)); abline(0,1,col=2)
plot(t01_indiv, apply(subset_samps("t01_indiv", c("sd", "raw"), samps), 2, mean)); abline(0,1,col=2)



