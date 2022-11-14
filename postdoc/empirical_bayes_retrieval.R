#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

#define functions
subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  if(all(exclude == "")){
    return(samps[,incl_inds])
  }else{
    excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
    return_inds <- setdiff(incl_inds, excl_inds)
    return(samps[,return_inds])  
  }
}
my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  sr <- cor(x, y)
  r <- abs(sr)
  txt <- format(c(sr, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#simulate data
n = 15
p = 5E3
b_var <- 0.5
b <- rnorm(p, 0, sqrt(b_var))
a <- rnorm(p)

x <- rnorm(n)
sigma2 <- 1

y <- t(matrix(rep(a, n), p, n)) + t(t(x)) %*% t(b) + matrix(rnorm(n*p, 0, sqrt(sigma2)), n, p)

#ols estimate
fits <- as.data.frame(do.call(rbind, lapply(1:p, function(i){ 
  fit <- lm(y[,i] ~ x)
  fit <- summary(fit)$coefficients[2,1:2]
  c(est = as.numeric(fit[1]), se = as.numeric(fit[2]))
})))

plot(b, fits$est); abline(0,1)
lm(fits$est ~ b)

sd(fits$est)

#indiv bayesian models
stan_program_0 <- "
data {
  int<lower=0> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;
  real b;
  real<lower=0> sd;
}
model {
  y ~ normal(a + b * x,sd);
}"
if(!exists("curr_stan_program_0") || stan_program_0!= curr_stan_program_0){
  curr_stan_program_0 <- stan_program_0
  f_0 <- write_stan_file(stan_program_0)  
}

foo <- function(i){
  system(sprintf('echo "%s"', paste0(i, collapse="")))
  d_0 <- list(n = n,
            x = x,
            y = y[,i])
  sink(tempfile())
  out_0 <- mod_0$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 1E3, data = d_0, parallel_chains = 4, adapt_delta = 0.8, refresh = 0)
  sink()
  samps_0 <- data.frame(as_draws_df(out_0$draws()))$b
  c(pmean = mean(samps_0), sd = sd(samps_0))
}


bayes_marfits <- parallel::mclapply(1:p, function(i) foo(i), mc.cores = 12)
bayes_marfits <- data.frame(do.call(rbind, bayes_marfits))

plot(fits$est, bayes_marfits$pmean); abline(0,1,lty=2,col=2)
plot(log(fits$se), log(bayes_marfits$sd)); abline(0,1,lty=2,col=2)

#fit bayesian model
d <- list(p = p,
          b = fits$est,
          se = fits$se)

stan_program <- "
data {
  int<lower=0> p;
  vector[p] b;
  vector[p] se;
}
parameters {
  real mu;
  vector[p] b_aug;
  real<lower=0> sd;
}
model {
  mu ~ normal(0,10);
  sd ~ exponential(0.1);
  b ~ normal(b_aug, se);
  b_aug ~ normal(mu, sd);
}
"

if(!exists("curr_stan_program") || stan_program!= curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

d_2 <- list(p = p,
          n = n,
          x = x,
          y = y)

stan_program_2 <- "
data {
  int<lower=0> p;
  int<lower=0> n;
  vector[n] x;
  real y[n,p];
}
parameters {
  real mu_b;
  real<lower=0> sd_b;
  vector[p] b;
  vector[p] a;
  real<lower=0> sd[p];
}
model {
  mu_b ~ normal(0,10);
  sd_b ~ exponential(0.1);
  sd ~ exponential(0.5);
  for(i in 1:p){
    y[,i] ~ normal(a[i] + b[i] * x, sd[i]);
  }
  b ~ normal(mu_b, sd_b);
}
"

if(!exists("curr_stan_program_2") || stan_program_2!= curr_stan_program_2){
  curr_stan_program_2 <- stan_program_2
  f_2 <- write_stan_file(stan_program_2)  
}
mod_2 <- cmdstan_model(f_2)


#fit model
out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, adapt_delta = 0.8)
check_fit <- F
if(check_fit){
  summ <- out$summary()
  summ[order(summ$ess_bulk),]
  summ[order(summ$rhat, decreasing = T),]
}
samps <- data.frame(as_draws_df(out$draws()))

b_samp <- subset_samps("b_aug", samps = samps)
b_means <- apply(b_samp, 2, mean)
pairs(cbind(b, b_est = fits$est, b_means), lower.panel = my_line, upper.panel = panel.cor)

qs <- sapply(1:ncol(b_samp), function(i) mean(b_samp[,i] > b[i]))
hist(qs, breaks = 0:20/20)

#compare to analytic fit
out_2 <- mod_2$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d_2, parallel_chains = 4, adapt_delta = 0.8)
if(check_fit){
  summ_2 <- out_2$summary()
  summ_2[order(summ_2$ess_bulk),]
  summ_2[order(summ_2$rhat, decreasing = T),]
}
samps_2 <- data.frame(as_draws_df(out_2$draws()))

b_samp_2 <- subset_samps("b", exclude = c("mu", "sd"), samps = samps_2)

probs <- 1:99/100
plot(quantile(b_samp_2[,1], probs), quantile(b_samp[,1], probs), type = "l",
     xlim = range(b_samp_2), ylim = range(b_samp))
for(i in 2:p){
  lines(quantile(b_samp_2[,i], probs), quantile(b_samp[,i], probs), type = "l")
}
abline(0,1, col = 2, lty = 2)

#check if tail calibration is good?
prop_greater_than_0 <- function(x) mean(x > 0)
pg0 <- apply(b_samp_2, 2, prop_greater_than_0)
alpha <- 0.05
hits <- which(pg0 < alpha | pg0 > (1-alpha))

CI_range <- 0:100/100
par(mfrow = c(2,2))

#everything
q.b <- sapply(1:p, function(i) mean(b_samp_2[,i] > b[i]))
hist(q.b, xlim = c(0,1), breaks = 0:20/20, freq = F, xlab = "quantile of true value in marginal posterior", main = "b_S")
CI_coverage.b <- sapply(CI_range, function(CIr) sum((sapply(q.b, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.b))
plot(CI_range, CI_coverage.b, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value (0)", main = "b_S")
abline(0,1,lty=2,lwd=2,col=2)
legend("topleft", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")

#just the hits
q.b <- sapply(hits, function(i) mean(b_samp_2[,i] > b[i]))
hist(q.b, xlim = c(0,1), breaks = 0:20/20, freq = F, xlab = "quantile of true value in marginal posterior", main = "b_S")
CI_coverage.b <- sapply(CI_range, function(CIr) sum((sapply(q.b, function(qi) sum((qi * (1-1E-6) + 0.5 * 1E-6) > c(0.5 - CIr / 2, 0.5 + CIr / 2))) == 1)) / length(q.b))
plot(CI_range, CI_coverage.b, type = "l", xlab = "breadth of middle credibility interval", ylab = "coverage of true parameter value (0)", main = "b_S")
abline(0,1,lty=2,lwd=2,col=2)
legend("topleft", col=2,lty=2,lwd=2, legend = "1-to-1 line", bty = "n")
