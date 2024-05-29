library(cmdstanr)
library(posterior)

#set a seed for reproducibility
set.seed <- 123

#simulate coefficients
r <- 0.75 #true correlation between coefficients
R <- diag(2) * (1-r) + r #corresponding correlation matrix
n <- 1E3 #total number of samples
x <- matrix(rnorm(n*2), ncol = 2) %*% chol(R) #sample true coefficients
e_sd_sd <- c(0.1, 1) #define differential power across two dimensions
e_sd <- matrix(rexp(n*2), ncol = 2) %*% diag(e_sd_sd) #sample element-wise error
e <- matrix(rnorm(n*2), ncol = 2) * e_sd #generate error matrix
x_obs <- x + e #generate observed coefficients

#evaluate sample correlations
cor(x)[1,2] #should be approximately r
cor(x_obs)[1,2] #should be regressed from r towards 0

#pack data into a list
dat <- list(n=n, x_obs=x_obs, e_sd=e_sd)

#specify model with uncertainty
model_string <-  "
data{
    int n;
    matrix[n, 2] x_obs;
    matrix<lower=0>[n, 2] e_sd;
}
parameters{
    array[n] vector[2] x;
    cholesky_factor_corr[2] L;
}
model{
    x_obs[,1] ~ normal(x[,1], e_sd[,1]);
    x_obs[,2] ~ normal(x[,2], e_sd[,2]);
    L ~ lkj_corr_cholesky(1); //not strictly necessary (already implied)
    x ~ multi_normal_cholesky(rep_vector(0, 2), L);
}
generated quantities{
    real r = multiply_lower_tri_self_transpose(L)[1,2];
}
"

#fit model
mod <- cmdstan_model(write_stan_file(model_string))
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  data = dat, adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 1)

#check convergence
summ <- fit$summary()
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws("r")))
mean(samps$r)

#plot marginal posterior
hist(samps$r, breaks = -50:50/50, main = "posterior distribution\nof correlation", freq = F, xlab = "correlation")

#label lines
abline(v = r, col = "red", lwd = 3)
text(x = r, y = par("usr")[4], pos = 3, labels = "true\ncorr.", xpd = NA, col = "red")
abline(v = cor(x_obs)[1,2], col = "blue", lwd = 3)
text(x = cor(x_obs)[1,2], y = par("usr")[4], 
     pos = 3, labels = "sample\ncorr.", xpd = NA, col = "blue")
