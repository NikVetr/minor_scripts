library(cmdstanr)
library(posterior)

#simple pairwise difference check
n <- 1E3
x <- rbinom(n, 1, 0.5)
y1 <- rnorm(n)
y2 <- rnorm(n)

f1 <- summary(lm(y1 ~ 1 + x))$coefficients
f2 <- summary(lm(y2 ~ 1 + x))$coefficients
f3 <- summary(lm(c(y1, y2) ~ 1 + c(x, x) * rep(c(0,1), each = n)))$coefficients

c(est = f2[2,1] - f1[2,1], SE = sqrt(f1[2,2]^2 + f2[2,2]^2))
f3[4,1:2]

#pairwise multi-category calibration check
n <- 5E2 #number of samples / individuals
p <- 1E2 #number of input variables, x
m <- 1E1 #number of categories eg cell types by which input variables affect output variables
x <- matrix(rbinom(n * p, 1, 0.5), n, p)
b_mu <- rnorm(p)
b_sd <- 0
b <- b_mu + matrix(rnorm(k*p) * b_sd, p, m)
bv <- c(b)
y_mu <- x %*% b
e <- matrix(rnorm(n*m), n, m) %*% diag(sqrt(apply(y_mu, 2, var)))
y <- y_mu + e

dat <- list(n = n,
            p = p,
            m = m,
            x = x,
            y = y)

#fit model with pooling within b_mu

stan_model_1 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> m;
  matrix[n,p] x;
  matrix[n,m] y;
}
parameters {
  vector[p] b_mu;
  real<lower=0> b_mu_sd;
  real<lower=0> b_sd;
  matrix[p, m] b;
  real<lower=0> sigma;
}
transformed parameters {
  
}
model {
  // priors //

  //scale parameters
  sigma ~ normal(0,3);
  b_mu_sd ~ std_normal();
  b_sd ~ std_normal();
  
  //coefficients
  b_mu ~ normal(0, b_mu_sd);
  for(k in 1:p){
    b[k,] ~ normal(b_mu[k], b_sd);
  }
  
  // model //
  for(i in 1:n){
    for(j in 1:m){
      real y_mu = sum(b[,j] .* to_vector(x[i,]));
      y[i,j] ~ normal(y_mu, sigma);  
    }
  }
}
"


mod_1 <- cmdstan_model(write_stan_file(stan_model_1))
fit_1 <- mod_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 1, init = 0.1, data = dat)

summ_1 <- fit_1$summary()
print(summ_1[order(summ_1$rhat, decreasing = T),])
samps_1 <- data.frame(as_draws_df(fit_1$draws("b")))

#fit model with pooling across b_mu

stan_model_2 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> m;
  matrix[n,p] x;
  matrix[n,m] y;
}
parameters {
  vector[m] b_mu;
  real<lower=0> b_mu_sd;
  real<lower=0> b_sd;
  matrix[p, m] b;
  real<lower=0> sigma;
}
transformed parameters {
  
}
model {
  // priors //

  //scale parameters
  sigma ~ normal(0,3);
  b_mu_sd ~ std_normal();
  b_sd ~ std_normal();
  
  //coefficients
  b_mu ~ normal(0, b_mu_sd);
  for(j in 1:m){
    b[,j] ~ normal(b_mu[j], b_sd);
  }
  
  // model //
  for(i in 1:n){
    for(j in 1:m){
      real y_mu = sum(to_row_vector(b[,j]) .* x[i,]);
      y[i,j] ~ normal(y_mu, sigma);  
    }
  }
}
"


mod_2 <- cmdstan_model(write_stan_file(stan_model_2))
fit_2 <- mod_2$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.9, parallel_chains = 4,
                      refresh = 100, max_treedepth = 10, 
                      thin = 1, init = 0.1, data = dat)

summ_2 <- fit_2$summary()
print(summ_2[order(summ_2$rhat, decreasing = T),])
samps_2 <- data.frame(as_draws_df(fit_2$draws("b")))

#evaluate calibration
bad_samps <- c(".chain", ".iteration", ".draw")
bs_1 <- samps_1[,!c(colnames(samps_1) %in% bad_samps)]
bs_2 <- samps_2[,!c(colnames(samps_2) %in% bad_samps)]
pairs(cbind(bv, bs1 = colMeans(bs_1), bs2 = colMeans(bs_2)))

plot(apply(bs_1, 2, sd), apply(bs_2, 2, sd)); abline(0,1)

breaks <- 0:20/20

#main effects
bcal_1 <- sapply(1:length(bv), function(i) mean(bv[i] > bs_1[,i]))
hist(bcal_1, breaks = breaks)
bcal_2 <- sapply(1:length(bv), function(i) mean(bv[i] > bs_2[,i]))
hist(bcal_2, breaks = breaks)

#pairwise comparisons
pcs <- t(combn(1:k, 2))
npcs <- nrow(pcs)
bcal_1 <- sapply(1:npcs, function(i){
  ms <- pcs[i,]
  sapply(1:p, function(j){
    true_diff <- (b[j,ms[1]] - b[j,ms[2]])  
    post_diff <-  (bs_1[,paste0("b.", j, ".", ms[1], ".")] - 
                     bs_1[,paste0("b.", j, ".", ms[2], ".")])
    mean(true_diff > post_diff)
  })
})
hist(bcal_1, breaks = breaks)

bcal_2 <- sapply(1:npcs, function(i){
  ms <- pcs[i,]
  sapply(1:p, function(j){
    true_diff <- (b[j,ms[1]] - b[j,ms[2]])  
    post_diff <-  (bs_2[,paste0("b.", j, ".", ms[1], ".")] - 
                     bs_2[,paste0("b.", j, ".", ms[2], ".")])
    mean(true_diff > post_diff)
  })
})
hist(bcal_2, breaks = breaks)
