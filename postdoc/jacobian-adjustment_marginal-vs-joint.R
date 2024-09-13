library(cmdstanr)
library(posterior)

#### functions ####

meanhist <- function(x, ...){hist(x, freq = F, main = paste0("mean = ", round(mean(x), 3), ", var = ", round(var(x), 3)), ...)}
panel.pts <- function(x, y, ...){
  points(x,y,pch = 19, col = adjustcolor(1, min(1, 1 / length(x) * 25)))
  rij <- round(cor(x, y), 3)
  lab <- paste0("cor = ", rij)
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
  text(x = par("usr")[1], y = par("usr")[4] - diff(par("usr")[3:4])/4, pos = 4, col = "darkgreen",
       labels = paste0("mean = ", round(mean(x), 3), "\nvar = ", round(var(x), 3)), cex = 1.25)
}


#### basic addition ####
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
    
    //jacobian adjustment
    target += log(abs(1)); // aka not necessary
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(chains = 4, 
                  iter_sampling = 1E4, 
                  iter_warmup = 1E4, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(fit$draws()))
pairs(samps[,c("a", "b", "c")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)


#### basic addition, MVN prior ####
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
    matrix[2,2] Sigma = [[1, 1], [1, 2]];
    [a, c] ~ multi_normal([0,0], Sigma);
    
    //jacobian adjustment
    target += log(abs(1)); // aka not necessary
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(chains = 4, 
                  iter_sampling = 5E3, 
                  iter_warmup = 5E3, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(fit$draws()))
pairs(samps[,c("a", "b", "c")], diag.panel = panel.hist, lower.panel = panel.pts, upper.panel = NULL, cex.labels = 5)


#### basic addition, fixing one marginal ####
n <- 1E2
a <- rnorm(n)
dat <- list(n=n, a=a)
stan_program <- '
data {
  int<lower=0> n;
  vector[n] a;
}
parameters {
  vector[n] b;
}
transformed parameters {
  vector[n] c = a + b;
}
model {
  // implicit marginal distribution if a ~ std_normal()
  c ~ normal(0,sqrt(2));
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(data = dat,
                  chains = 4, 
                  iter_sampling = 1E3, 
                  iter_warmup = 1E3, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, 
                  max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
badcols <- c(".chain", ".iteration", ".draw")
bsamps <- data.frame(as_draws_df(fit$draws("b")))
csamps <- data.frame(as_draws_df(fit$draws("c")))
bms <- apply(bsamps[,!(colnames(bsamps) %in% badcols)], 2, mean)
cms <- apply(csamps[,!(colnames(csamps) %in% badcols)], 2, mean)
out <- data.frame(a = a, b = bms, c = cms)
pairs(out, 
      diag.panel = panel.hist, 
      lower.panel = panel.pts, upper.panel = NULL, cex.labels = 3)



#### basic addition, fixing one marginal, prior on sample r ####
n <- 1E2
a <- rnorm(n)
dat <- list(n=n, a=a)
stan_program <- '
data {
  int<lower=0> n;
  vector[n] a;
}
parameters {
  vector[n] b;
}
transformed parameters {
  vector[n] c = a + b;
  real r = dot_product(a, b) / (n - 1) / sd(a) / sd(b);
  real t = r * sqrt((n-2.0)/1.0-r^2);
}
model {
  // implicit marginal distribution if a ~ std_normal()
  c ~ normal(0,sqrt(2));
  t ~ student_t(n-2, 0, sqrt((n-2.0)/(n-4.0)));
  //t ~ student_t(n-2, 0, sqrt((n-2.0)/(n-4.0))/100); //really force it to 0
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(data = dat,
                  chains = 4, 
                  iter_sampling = 1E3, 
                  iter_warmup = 1E3, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, 
                  max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
badcols <- c(".chain", ".iteration", ".draw")
samps <- data.frame(as_draws_df(fit$draws()))
bsamps <- data.frame(as_draws_df(fit$draws("b")))
csamps <- data.frame(as_draws_df(fit$draws("c")))
bms <- apply(bsamps[,!(colnames(bsamps) %in% badcols)], 2, mean)
cms <- apply(csamps[,!(colnames(csamps) %in% badcols)], 2, mean)
out <- data.frame(a = a, b = bms, c = cms)
pairs(out, 
      diag.panel = panel.hist, 
      lower.panel = panel.pts, upper.panel = NULL, cex.labels = 3)



#### basic addition, MVN prior, fixing one marginal ####
n <- 2E2
a <- rnorm(n)
dat <- list(n=n, a=a)
stan_program <- '
data {
  int<lower=0> n;
  vector[n] a;
}
parameters {
  vector[n] b;
}
transformed parameters {
  vector[n] c = a + b;
  array[n] row_vector[2] ac;
  for(i in 1:n){
    ac[i] = [a[i], c[i]];
  }
}
model {
  matrix[2,2] Sigma = [[1, 1], [1, 2]];
  ac ~ multi_normal([0,0], Sigma);
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(data = dat,
                  chains = 4, 
                  iter_sampling = 1E3, 
                  iter_warmup = 1E3, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, 
                  max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
badcols <- c(".chain", ".iteration", ".draw")
bsamps <- data.frame(as_draws_df(fit$draws("b")))
csamps <- data.frame(as_draws_df(fit$draws("c")))
bms <- apply(bsamps[,!(colnames(bsamps) %in% badcols)], 2, mean)
cms <- apply(csamps[,!(colnames(csamps) %in% badcols)], 2, mean)
out <- data.frame(a = a, b = bms, c = cms)
pairs(out, 
      diag.panel = panel.hist, 
      lower.panel = panel.pts, upper.panel = NULL, cex.labels = 3)

#### basic addition, noncentered parameterization, fixing one marginal ####
n <- 2E2
a <- rnorm(n)
dat <- list(n=n, a=a)
stan_program <- '
data {
  int<lower=0> n;
  vector[n] a;
}
transformed data {
  matrix[2,2] sigma = [[1, 1], [1, 2]];
  matrix[2,2] sigma_chol_inv = [[1, -1], [0, 1]];
}
parameters {
  vector[n] b;
}
transformed parameters {
  vector[n] c = a + b;
  matrix[n, 2] ac;
  ac[,1] = a;
  ac[,2] = c;
  matrix[n, 2] ac_norm = ac * sigma_chol_inv;
}
model {
  ac_norm[,1] ~ std_normal();
  ac_norm[,2] ~ std_normal();
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
fit <- mod$sample(data = dat,
                  chains = 4, 
                  iter_sampling = 1E3, 
                  iter_warmup = 1E3, 
                  parallel_chains = 4, 
                  adapt_delta = 0.85, 
                  refresh = 10, 
                  init = 1, 
                  max_treedepth = 15)
summ <- fit$summary()
summ[order(summ$ess_bulk),]

#check inference
badcols <- c(".chain", ".iteration", ".draw")
bsamps <- data.frame(as_draws_df(fit$draws("b")))
bsamps <- bsamps[,!(colnames(bsamps) %in% badcols)]
csamps <- data.frame(as_draws_df(fit$draws("c")))
csamps <- csamps[,!(colnames(csamps) %in% badcols)]

bms <- apply(bsamps, 2, mean)
cms <- apply(csamps, 2, mean)
out <- data.frame(a = a, b = bms, c = cms)


ab_corr <- sapply(1:nrow(bsamps), function(i) cor(a, unlist(bsamps[i,])))
ac_corr <- sapply(1:nrow(csamps), function(i) cor(a, unlist(csamps[i,])))
hist(apply(csamps, 1, var))
hist(apply(bsamps, 1, var))
pairs(out, 
      diag.panel = panel.hist, 
      lower.panel = panel.pts, upper.panel = NULL, cex.labels = 3)
