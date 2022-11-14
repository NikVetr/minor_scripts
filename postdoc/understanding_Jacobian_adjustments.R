#TODO add covariance? user effects w/ repeated measures?

#libaries
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#functions
meanhist <- function(x, ...){hist(x, freq = F, main = paste0("mean = ", round(mean(x), 3), ", var = ", round(var(x), 3)), ...)}
panel.pts <- function(x, y, ...){
  points(x,y,pch = 19, col = adjustcolor(1, 0.25))
  rij <- round(cov(x, y), 3)
  lab <- paste0("cov = ", rij)
  if(rij >= 0) {
    text(x = par("usr")[1] + strwidth(lab, cex = 2) / 1.9, y = par("usr")[4], pos = 1,
         labels = lab, cex = 2, xpd = NA, col = "darkred")
  } else {
    text(x = par("usr")[2] - strwidth(lab, cex = 2) / 1.9, y = par("usr")[4], pos = 1,
         labels = lab, cex = 2, xpd = NA, col = "darkblue")  
  }
}
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = F)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(xleft = breaks[-nB], ybottom = 0, xright = breaks[-1], ytop = y, col = "grey50", ...)
  text(x = mean(par("usr")[1:2]), y = 1.5, pos = 1, col = "darkgreen",
       labels = paste0("mean = ", round(mean(x), 3), ", var = ", round(var(x), 3)), cex = 2)
}

#basic addition
stan_program <- '
parameters {
    real a;
    real b;
}
transformed parameters {
    real c = a + b;
}
model {
    // priors
    a ~ std_normal();
    c ~ normal(0,sqrt(2));
    target += log(fabs(1));
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, parallel_chains = 4, adapt_delta = 0.85, refresh = 10, init = 1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
pairs(samps[,c("a", "b", "c")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)


#basic addition, alternative bijection
stan_program <- '
parameters {
    real a;
    real b;
}
transformed parameters {
    real c = a + b;
}
model {
    // priors
    a ~ std_normal();
    c ~ normal(0,sqrt(2));
    target += log(fabs(a-b));
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, parallel_chains = 4, adapt_delta = 0.85, refresh = 10, init = 1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
pairs(samps[,c("a", "b", "c")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)

#squaring
stan_program <- '
parameters {
    real a;
}
transformed parameters {
    real b = a^2;
}
model {
    b ~ chi_square(1);
    target += log(fabs(2*a));
}
'
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E4, iter_warmup = 1E4, parallel_chains = 4, adapt_delta = 0.85, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
pairs(samps[,c("a", "b")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)
mean(rchisq(1E4, 1))
var(rchisq(1E4, 1))

#product of two parameters? see https://www.tamaspapp.eu/post/jacobian-chain/
stan_program <- '
parameters {
    real<lower=0> a;
    real<lower=0> b;
}
transformed parameters {
    real<lower=0> c = a*b;
}
model {
    a ~ lognormal(0, 1);
    c ~ lognormal(0, sqrt(2));
    target += log(fabs(b));
}
'
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, parallel_chains = 4, adapt_delta = 0.85, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
pairs(log(samps[,c("a", "b", "c")]), diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)


#product of two parameters, w/ normal priors?
stan_program <- '
parameters {
    real a;
    real b;
}
transformed parameters {
    real c = a*b;
}
model {
    a ~ std_normal();
    b ~ uniform(-1000000, 1000000);
    c ~ std_normal();
    target += log(fabs(a));
}
'
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, parallel_chains = 4, adapt_delta = 0.99, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
pairs(samps[,c("a", "b", "c")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)
pairs(log(abs(samps[,c("a", "b", "c")])), diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)



#product of three parameters?
stan_program <- '
parameters {
    real a;
    real b;
    real c;
}
transformed parameters {
    real d = a*b*c;
}
model {
    a ~ std_normal();
    b ~ std_normal();
    d ~ std_normal();
    target += log(fabs(a*b));
}
'
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, parallel_chains = 4, adapt_delta = 0.85, refresh = 10, init = 0.1, max_treedepth = 15)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
meanhist(samps$a)
meanhist(samps$b)
meanhist(samps$c)
meanhist(samps$d)
