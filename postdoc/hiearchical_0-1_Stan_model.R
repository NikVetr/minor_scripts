library(cmdstanr)
library(posterior)
library(caret)
library(MASS)
library(mvtnorm)


k <- 3
n_k <- 500
n <- n_k * k
n_j <- 4
mean_corrs <- seq(-0.7, 0.7, length.out = k)
mean_Betas <- (mean_corrs + 1) / 2
concentrations <- 1:k*5 + 10
corrs <- do.call(c, lapply(1:k, function(ki) rbeta(n = n_k, 
                                                       shape1 = mean_Betas[ki] * concentrations[ki], 
                                                       shape2 = (1-mean_Betas[ki]) * concentrations[ki]) * 2 - 1))
x_list <- lapply(1:n, function(indiv) rmvnorm(n_j, c(0,0), sigma = diag(2) * (1-corrs[indiv]) + corrs[indiv]))
x <- abind::abind(x_list, along = 3)

par(mfrow = c(3,3))
plot(corrs, sapply(1:dim(x)[3], function(xi) cor(x[,,xi])[1,2]),
     main = paste0("Pearson's r = ", round(cor(corrs, sapply(1:dim(x)[3], function(xi) cor(x[,,xi])[1,2])), 3)),
     xlab = "true correlations", ylab = "sample correlations")
abline(0,1,col=2,lwd=2)
plot(-1000:1000/1000, apply(sapply(1:k, function(ki) dbeta(0:2000/2000, shape1 = mean_Betas[ki] * concentrations[ki], 
                                                    shape2 = (1-mean_Betas[ki]) * concentrations[ki])), 1, mean), type = "l",
     xlab = "correlation", ylab = "true mixture density")
hist(sapply(1:dim(x)[3], function(xi) cor(x[,,xi])[1,2]), breaks = -10:10/10, xlab = "correlation", ylab = "sample density",
     main = "sample correlations", freq = F)

d <- list(n = n,
          k = k,
          n_j = n_j,
          x = x_list)

# STAN model
stan_program <- '
data {
    int<lower=1> n; //total number of individuals
    int<lower=3> n_j; //number of observations per individual
    int<lower=1> k; //number of components in the Beta mixture
    //matrix[n_j,2] x[n];
    array[n, n_j] vector[2] x;
}
transformed data {
    vector<lower=0>[k] mixing_alphas = rep_vector(1,k);
    vector[2] mu = rep_vector(0,2);
    vector<lower=0>[2] sd = rep_vector(1,2);
    real<lower=0} concentration_offset = 10;
}
parameters {
    vector<lower=0, upper=1>[n] squished_corrs;
    ordered[k] relative_mean_squished_corrs;
    vector<lower=0>[k] concentration_squished_corrs;
    simplex[k] theta;
}
transformed parameters {
    vector<lower=0,upper=1>[k] mean_squished_corrs = inv_logit(relative_mean_squished_corrs);
    vector<lower=-1, upper=1>[n] corrs = squished_corrs * 2 - 1;
    cholesky_factor_cov[2] L_Sigma[n];
    for(i in 1:n){
      L_Sigma[i] = diag_pre_multiply(sd, cholesky_decompose([[1, corrs[i]], [corrs[i], 1]]));
    }
    
    vector[k] log_theta = log(theta);
    vector[k] alpha = mean_squished_corrs .* (concentration_offset+concentration_squished_corrs);
    vector[k] beta = (1 - mean_squished_corrs) .* (concentration_offset+concentration_squished_corrs);
}

model {
    //priors
    relative_mean_squished_corrs ~ normal(0, 1.8); //1.8 set bc it implies an approx Uniform(-1,1) means
    theta ~ dirichlet(mixing_alphas);
    concentration_squished_corrs ~ exponential(0.1);
    
    //likelihood for obs
    for(i in 1:n){
      x[i] ~ multi_normal_cholesky(mu, L_Sigma[i]);
    }
    
    //Beta mixture model
    for (i in 1:n) {
      vector[k] lps = log_theta;
      for(l in 1:k){
        lps[l] += beta_lpdf(squished_corrs[i] | alpha[l], beta[l]);
      }
      target += log_sum_exp(lps);
    }
}
'


# stan_program <- stan_program_direct
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
# out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = lapply(1:4, function(x) list(sigma = rep(3,n_dim))))
out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = 0.1, max_treedepth = 20)
# out <- mod$variational(data = d)
summ <- out$summary()
summ[order(summ$ess_bulk),]
samps <- data.frame(as_draws_df(out$draws()))

plot(corrs, apply(samps[,paste0("corrs.", 1:(k*n_k), ".")], 2, mean),
     main = paste0("Pearson's r = ", round(cor(corrs, apply(samps[,paste0("corrs.", 1:(k*n_k), ".")], 2, mean)), 3)),
     xlab = "true correlations", ylab = "posterior mean correlations")
abline(0,1,col=2,lwd=2)

apply(samps[,paste0("alpha")], 2, mean)


hist(samps$mean_squished_corrs.1. * 2 - 1, xlim = c(-1,1), breaks = -20:20/20, col = adjustcolor(1, 0.5),
     xlab = "correlation", ylab = "density", freq = F, main = "posterior means of beta component means")
abline(v = mean_corrs[1], lwd = 2, col = 1)
for(i in 2:k){
  hist(samps[,paste0("mean_squished_corrs.",i,".")] * 2 - 1, add = T, breaks = -20:20/20, col = adjustcolor(i, 0.5), freq = F)
  abline(v = mean_corrs[i], lwd = 2, col = i)  
}

hist(samps$concentration_squished_corrs.1., 
     xlim = range(sapply(1:k, function(ki) samps[,paste0("concentration_squished_corrs.", ki, ".")] )), 
     breaks = seq(-1 + min(sapply(1:k, function(ki) samps[,paste0("concentration_squished_corrs.", ki, ".")] )),
                  1 + max(sapply(1:k, function(ki) samps[,paste0("concentration_squished_corrs.", ki, ".")] )),
                  length.out = 20), col = adjustcolor(1, 0.5), ylim = c(0,0.1),
     xlab = "concentration", ylab = "density", freq = F, main = "posterior means of beta component conc.")
abline(v = concentrations[1], lwd = 2, col = 1)  
for(i in 2:k){
  hist(samps[,paste0("concentration_squished_corrs.",i,".")], add = T, 
       breaks = seq(-1 + min(sapply(1:k, function(ki) samps[,paste0("concentration_squished_corrs.", ki, ".")] )),
                    1 + max(sapply(1:k, function(ki) samps[,paste0("concentration_squished_corrs.", ki, ".")] )), length.out = 20),
                    col = adjustcolor(i, 0.5), freq = F)
  abline(v = concentrations[i], lwd = 2, col = i)  
}


# plot(sapply(1:dim(x)[3], function(xi) cor(x[,,xi])[1,2]), apply(samps[,paste0("corrs.", 1:(k*n_k), ".")], 2, mean),
#      main = paste0("Pearson's r = ", round(cor(sapply(1:dim(x)[3], function(xi) cor(x[,,xi])[1,2]), apply(samps[,paste0("corrs.", 1:(k*n_k), ".")], 2, mean)), 3)),
#      xlab = "sample correlations", ylab = "posterior mean correlations")
# abline(0,1,col=2,lwd=2)
