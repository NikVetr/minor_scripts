#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#simulate data
n_rep <- 1E2
n <- 20
sd_a <- 0.5
sd_b <- 0.1
sd_e <- 1
a <- rnorm(n_rep, sd = 0.5)
b <- rnorm(n_rep, sd = 0.1)
# a <- rep(0, n_rep)
# b <- rep(0, n_rep)
x <- replicate(n_rep, rnorm(n))
y <- sapply(1:n_rep, function(i) a[i] + b[i] * x[,i] + rnorm(n, sd = sd_e))
fits <- data.frame(t(sapply(1:n_rep, function(i){
  fit <- lm(y[,i] ~ x[,i])
  c(summary(fit)$coefficients, summary(fit)$r.squared)
})))
colnames(fits) <- c("est_a", "est_b", "se_a", "se_b", "t_a", "t_b", "p_a", "p_b", "r2")
sd(fits$est_a)
sd(fits$est_b)
# plot(abs(b), fits$r2)

#fit bayesian model
d <- list(n = n,
          n_rep = n_rep,
          x = t(x),
          y = t(y)
)

stan_program <- "
data {
  int<lower=1> n;
  int<lower=1> n_rep;
  matrix[n_rep, n] x;
  matrix[n_rep, n] y;
}
parameters {
  vector[n_rep] a;
  vector[n_rep] b;
  vector<lower=0>[n_rep] sigma_e;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  
}
model {
  a ~ normal(0,sigma_a);
  sigma_a ~ normal(0,1);
  b ~ normal(0,sigma_b);
  sigma_b ~ normal(0,1);
  sigma_e ~ normal(0,1);
  for(i in 1:n_rep){
    y[i,] ~ normal(a[i] + x[i,] * b[i], sigma_e[i]);
  }
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

subset_samps <- function(include = NULL, exclude = NULL, samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}

plot(a, apply(subset_samps("a\\.", "sigma", samps), 2, mean)); abline(0,1,col=2)
plot(a, fits$est_a); abline(0,1,col=2)
plot(fits$est_a, apply(subset_samps("a\\.", "sigma", samps), 2, mean)); abline(0,1,col=2)

plot(b, apply(subset_samps("b\\.", "sigma", samps), 2, mean)); abline(0,1,col=2)
plot(b, fits$est_b); abline(0,1,col=2)
plot(fits$est_b, apply(subset_samps("b\\.", "sigma", samps), 2, mean)); abline(0,1,col=2)

mass_to_one_side_of_0 <- function(x) max(mean(x>0), mean(x<0))
hist(apply(subset_samps("a\\.", "sigma", samps), 2, mass_to_one_side_of_0), breaks = 20)
hist(apply(subset_samps("b\\.", "sigma", samps), 2, mass_to_one_side_of_0), breaks = 20)


#fit post-hoc bayesian model
d_ph <- list(est_a = fits$est_a,
          est_b = fits$est_b,
          se_a = fits$se_a,
          se_b = fits$se_b,
          n_rep = n_rep
)

stan_program_ph <- "
data {
  int<lower=1> n_rep;
  vector[n_rep] est_a;
  vector<lower=0>[n_rep] se_a;
  vector[n_rep] est_b;
  vector<lower=0>[n_rep] se_b;
}
parameters {
  vector[n_rep] a;
  vector[n_rep] b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
}
model {
  a ~ normal(0,sigma_a);
  sigma_a ~ normal(0,1);
  b ~ normal(0,sigma_b);
  sigma_b ~ normal(0,1);
  est_a ~ normal(a,se_a);
  est_b ~ normal(b,se_b);
}
"

if(!exists("curr_stan_program_ph") || stan_program_ph!= curr_stan_program_ph){
  curr_stan_program_ph <- stan_program_ph
  f_ph <- write_stan_file(stan_program_ph)
}
mod_ph <- cmdstan_model(f_ph)

out_ph <- mod_ph$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d_ph, parallel_chains = 4, adapt_delta = 0.95)
# summ_ph <- out_ph$summary()
# summ_ph[order(summ_ph$ess_bulk),]
# summ_ph[order(summ_ph$rhat, decreasing = T),]

samps_ph <- data.frame(as_draws_df(out_ph$draws()))

mean(samps_ph$sigma_a)
mean(samps_ph$sigma_b)

mass_to_one_side_of_0 <- function(x) max(mean(x>0), mean(x<0))


hist(fits$p_a, breaks = 20)
hist(p.adjust(fits$p_a, method = "BH"), breaks = 20, xlim = c(0,1))
hist(p.adjust(fits$p_a, method = "bonf"), breaks = 20, xlim = c(0,1))
hist(apply(subset_samps("a\\.", "sigma", samps_ph), 2, mass_to_one_side_of_0), breaks = 20)

hist(fits$p_b, breaks = 20)
hist(p.adjust(fits$p_b, method = "BH"), breaks = 20, xlim = c(0,1))
hist(p.adjust(fits$p_b, method = "bonf"), breaks = 20, xlim = c(0,1))
hist(apply(subset_samps("b\\.", "sigma", samps_ph), 2, mass_to_one_side_of_0), breaks = 20)

plot(apply(subset_samps("a\\.", "sigma", samps), 2, mean), apply(subset_samps("a\\.", "sigma", samps_ph), 2, mean)); abline(0,1,col=2)
plot(apply(subset_samps("b\\.", "sigma", samps), 2, mean), apply(subset_samps("b\\.", "sigma", samps_ph), 2, mean)); abline(0,1,col=2)


plot(apply(subset_samps("a\\.", "sigma", samps), 2, mass_to_one_side_of_0), 
     apply(subset_samps("a\\.", "sigma", samps_ph), 2, mass_to_one_side_of_0)); abline(0,1,col=2)
plot(apply(subset_samps("b\\.", "sigma", samps), 2, mass_to_one_side_of_0), 
     apply(subset_samps("b\\.", "sigma", samps_ph), 2, mass_to_one_side_of_0)); abline(0,1,col=2)

#evaluate calibration
eval_cal <- function(subsamps, truevals){
  qis <- 0:100/100
  qs <- apply(subsamps, 2, quantile, probs = qis)
  cis <- cbind(50:0 + 1, 50 + 1:51)
  calib <- sapply(1:nrow(cis), function(i) mean(qs[cis[i,1],] < truevals & qs[cis[i,2],] > truevals))
  return(cbind(apply(cis, 1, diff) / 100, calib))
}
plot(eval_cal(subset_samps("b\\.", "sigma", samps), b), type = "l", xlab = "interval width", ylab = "prop true vals in interval"); abline(0,1,col=2)
plot(eval_cal(subset_samps("a\\.", "sigma", samps), a), type = "l", xlab = "interval width", ylab = "prop true vals in interval"); abline(0,1,col=2)
plot(eval_cal(subset_samps("b\\.", "sigma", samps_ph), b), type = "l", xlab = "interval width", ylab = "prop true vals in interval"); abline(0,1,col=2)
plot(eval_cal(subset_samps("a\\.", "sigma", samps_ph), a), type = "l", xlab = "interval width", ylab = "prop true vals in interval"); abline(0,1,col=2)

#fit unpooled model one at a time, then try to empirical bayes it?
stan_program_single <- "
data {
  int<lower=1> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;
  real b;
  real<lower=0> sigma_e;
}
model {
  a ~ normal(0,100);
  b ~ normal(0,100);
  sigma_e ~ normal(0,100);
  y ~ normal(a + b * x, sigma_e);
}
"


if(!exists("curr_stan_program_single") || stan_program_single!= curr_stan_program_single){
  curr_stan_program_single <- stan_program_single
  f_single <- write_stan_file(stan_program_single)
}
mod_single <- cmdstan_model(f_single)

out_single <- data.frame(do.call(rbind, parallel::mclapply(1:n_rep, function(i){
  system(sprintf('echo "%s"', paste0(i, " ")))
  d_single <- list(n = n,
                   x = x[,i],
                   y = y[,i]
  )
  out_single <- (mod_single$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                                  data = d_single, parallel_chains = 2, adapt_delta = 0.8))
  samps_single <- data.frame(as_draws_df(out_single$draws()))
  c(est_a = mean(samps_single$a), est_b = mean(samps_single$b), sd_a = sd(samps_single$a), sd_b= sd(samps_single$b))
}, mc.cores = 16)))

par(mfrow = c(2,2))
plot(out_single$est_a, fits$est_a); abline(0,1)
plot(out_single$sd_a, fits$se_a); abline(0,1)
plot(out_single$est_b, fits$est_b); abline(0,1)
plot(out_single$sd_b, fits$se_b); abline(0,1)

#### ok, now let's check to see if we can do some sort of multilevel empirical bayes maximum marginal likelihood thing ####
stan_program_multilevel <- "
data {
  int<lower=1> n;
  int<lower=1> n_rep;
  matrix[n_rep, n] x;
  matrix[n_rep, n] y;
}
parameters {
  vector[n_rep] a;
  vector[n_rep] b;
  vector<lower=0>[n_rep] sigma_e;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  
}
model {
  a ~ normal(0,sigma_a);
  sigma_a ~ normal(0,1);
  b ~ normal(0,sigma_b);
  sigma_b ~ normal(0,1);
  sigma_e ~ normal(0,1);
  for(i in 1:n_rep){
    y[i,] ~ normal(a[i] + x[i,] * b[i], sigma_e[i]);
  }
}
"
mod_multilevel <- cmdstan_model(write_stan_file(stan_program_multilevel))


#fit model
out_multilevel <- mod_multilevel$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
samps_multilevel <- data.frame(as_draws_df(out_multilevel$draws()))

#test if sd(a) is equal to sd_a
sd(apply(subset_samps("a\\.", samps = samps_multilevel), 2, mean))
sd(unlist(subset_samps("a\\.", samps = samps_multilevel)))
mean(samps_multilevel$sigma_a)
sd_a

sd(apply(subset_samps("b\\.", samps = samps_multilevel), 2, mean))
sd(unlist(subset_samps("b\\.", samps = samps_multilevel)))
mean(samps_multilevel$sigma_b)
sd_b

hist(samps_multilevel$sigma_a)
abline(v = sd(apply(subset_samps("a\\.", samps = samps_multilevel), 2, mean)), col = 2, lwd = 2)

hist(samps_multilevel$sigma_b)
abline(v = sd(apply(subset_samps("b\\.", samps = samps_multilevel), 2, mean)), col = 2, lwd = 2)

stan_program_singlelevel <- "
data {
  int<lower=1> n;
  int<lower=1> n_rep;
  matrix[n_rep, n] x;
  matrix[n_rep, n] y;
}
parameters {
  vector[n_rep] a;
  vector[n_rep] b;
  vector<lower=0>[n_rep] sigma_e;
}
model {
  a ~ normal(0,100);
  b ~ normal(0,100);
  sigma_e ~ normal(0,100);
  for(i in 1:n_rep){
    y[i,] ~ normal(a[i] + x[i,] * b[i], sigma_e[i]);
  }
}
"

mod_singlelevel <- cmdstan_model(write_stan_file(stan_program_singlelevel))

#fit model
out_singlelevel <- mod_singlelevel$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
samps_singlelevel <- data.frame(as_draws_df(out_singlelevel$draws()))

sd(apply(subset_samps("a\\.", samps = samps_singlelevel), 2, mean))
sd(unlist(subset_samps("a\\.", samps = samps_singlelevel)))

sd(apply(subset_samps("a\\.", samps = samps_multilevel), 2, mean))
sd(unlist(subset_samps("a\\.", samps = samps_multilevel)))

sd(apply(subset_samps("b\\.", samps = samps_singlelevel), 2, mean))
sd(unlist(subset_samps("b\\.", samps = samps_singlelevel)))

sd(apply(subset_samps("b\\.", samps = samps_multilevel), 2, mean))
sd(unlist(subset_samps("b\\.", samps = samps_multilevel)))
sd_b

mod_singlelevel <- cmdstan_model(write_stan_file(stan_program_singlelevel))

#fit model
out_singlelevel_vb <- mod_singlelevel$variational(data = d, adapt_engaged = T)
samps_singlelevel_vb <- data.frame(as_draws_df(out_singlelevel_vb$draws()))


sd(apply(subset_samps(include = "a\\.", samps = samps_singlelevel_vb), 2, mean))
sd(unlist(subset_samps("a\\.", samps = samps_singlelevel_vb)))
sd_a

sd(apply(subset_samps(include = "b\\.", samps = samps_singlelevel_vb), 2, mean))
sd(unlist(subset_samps("b\\.", samps = samps_singlelevel_vb)))
sd(apply(subset_samps(include = "b\\.", samps = samps_multilevel), 2, mean))
sd(unlist(subset_samps("b\\.", samps = samps_multilevel)))
mean(subset_samps("sigma_b", samps = samps_multilevel))
sd_b

sd(fits$est_b)


stan_program_vb <- "
data {
  int<lower=1> n;
  int<lower=1> n_rep;
  matrix[n_rep, n] x;
  matrix[n_rep, n] y;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
}
parameters {
  vector[n_rep] a;
  vector[n_rep] b;
  vector<lower=0>[n_rep] sigma_e;
}
model {
  a ~ normal(0, sigma_a);
  b ~ normal(0, sigma_b);
  sigma_e ~ normal(0, 100);
  for(i in 1:n_rep){
    y[i,] ~ normal(a[i] + x[i,] * b[i], sigma_e[i]);
  }
}
"

mod_vb <- cmdstan_model(write_stan_file(stan_program_vb))

d_vb <- list(n = n,
          n_rep = n_rep,
          x = t(x),
          y = t(y),
          sigma_a = 100,
          sigma_b = 100
)

out_vb <- mod_vb$variational(data = d_vb, adapt_engaged = T)

est_sigmas <- list(
  sigma_a = sd(unlist(subset_samps("a\\.", samps = data.frame(as_draws_df(out_vb$draws()))))),
  sigma_b = sd(unlist(subset_samps("b\\.", samps = data.frame(as_draws_df(out_vb$draws())))))
)
d_vb <- list(n = n,
             n_rep = n_rep,
             x = t(x),
             y = t(y),
             sigma_a = est_sigmas$sigma_a,
             sigma_b = est_sigmas$sigma_b
)
out_vb <- mod_vb$variational(data = d_vb, adapt_engaged = T)
est_sigmas

mean(samps_multilevel$sigma_a)
mean(samps_multilevel$sigma_b)
est_sigmas


#need to subtract out posterior uncertainty -- whoops!
#see file '~/scripts/minor_scripts/postdoc/scale_mixture_variances.R' for derivation + demonstration

#