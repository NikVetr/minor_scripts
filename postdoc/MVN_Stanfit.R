library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#specify a few functions
rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

#### start of with known MVN observations ####

# Simulate some data
n_dim <- 10
mu <- rnorm(n_dim, 0, sd = 3)
R <- rlkj(n_dim, 1)
sds <- rexp(n_dim, rate = 0.5)

# R <- diag(n_dim) + 0.9 - diag(n_dim) * 0.9
# sds <- rep(2, n_dim)
# mu <- c(rep(0, n_dim-1), rep(5, 1))

Sig <- diag(sds) %*% R %*% diag(sds)
N <- 50
dat <- mvrnorm(N, mu, Sig)

# Randomly remove some data
frac_missing <- 0
n_miss <- N*n_dim*frac_missing
miss <- sample.int(prod(dim(dat)), size = n_miss)
dat[miss] <- NA

# Extract the missing values into a VECTOR
dat_complete <- dat[!is.na(dat)]

# Extract the missing and present values as MATRICES
ind_pres <- which(!is.na(dat), arr.ind = TRUE)
ind_miss <- which(is.na(dat), arr.ind = TRUE)

# STAN model
stan_program <- '
data {

    int<lower=0> nrow;
    int<lower=0> ncol;
    int<lower=0> n_comp; // Number of non-missing values
    int<lower=0> n_miss; // Number of missing values
    vector[n_comp] dat_complete;   // Vector of non-missing values
    array[n_comp, 2] int ind_pres;     // Matrix (row, col) of non-missing value indices
    array[n_miss, 2] int ind_miss;     // Matrix (row, col) of missing value indices
    
}
parameters {

    // MVN distribution parameters
    vector[ncol] mu;
    cholesky_factor_corr[ncol] L;
    vector<lower=0>[ncol] sigma; 
    
    // missing values
    vector[n_miss] ymiss;
    
}
transformed parameters {

    matrix[nrow, ncol] y;
    // Fill y with non-missing values 
    for(n in 1:n_comp) {
        y[ind_pres[n,1],ind_pres[n,2]] = dat_complete[n];
    }
    
    // Fill the rest of y with missing value "parameters"
    for(n in 1:n_miss){
        y[ind_miss[n,1],ind_miss[n,2]] = ymiss[n];
    }
    
    //matrix[nrow, ncol] y_noncentered_meanSD;   // now try non-centering these data
    //matrix[nrow, ncol] y_noncentered;   // now try non-centering these data
    //for(i in 1:nrow){
    //    y_noncentered_meanSD[i,] = (y[i,] - mu) ./ sigma;
    //}
    //y_noncentered = mdivide_left_tri_low(L, y_noncentered_meanSD)

}
model {

    sigma ~ exponential(1); // prior on the standard deviations
    L ~ lkj_corr_cholesky(2);
    for(i in 1:nrow){
      y[i,] ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L));
    }

    
}
generated quantities {

  matrix[ncol, ncol] corr_mat; 
  corr_mat = L * L\';
  
}
'

# # STAN model
# stan_program_direct <- '
# data {
#     int<lower=0> nrow;
#     int<lower=0> ncol;
#     int<lower=0> n_comp; // Number of non-missing values
#     int<lower=0> n_miss; // Number of missing values
#     real dat_complete[n_comp];   // Vector of non-missing values
#     int ind_pres[n_comp, 2];     // Matrix (row, col) of non-missing value indices
#     int ind_miss[n_miss, 2];     // Matrix (row, col) of missing value indices
# }
# parameters {
#     // Multivariate normal distribution parameters
#     vector[ncol] mu;
#     corr_matrix[ncol] corr_mat; 
#     vector<lower=0>[ncol] sigma; 
#     
#     // missing data
#     real ymiss[n_miss];      
# }
# transformed parameters {
#     
#     // generate covariance matrix
#     cov_matrix[ncol] Sigma; 
#     Sigma = quad_form_diag(corr_mat, sigma);
#     
#     matrix[nrow, ncol] y;
#     // Fill y with non-missing values 
#     for(n in 1:n_comp) {
#         y[ind_pres[n,1],ind_pres[n,2]] = dat_complete[n];
#     }
#     
#     // Fill the rest of y with missing value "parameters"
#     for(n in 1:n_miss){
#         y[ind_miss[n,1],ind_miss[n,2]] = ymiss[n];
#     }
# }
# model {
#     sigma ~ cauchy(0, 5); // prior on the standard deviations
#     corr_mat ~ lkj_corr(1); // LKJ prior on the correlation matrix
#     for(i in 1:nrow){
#         y[i] ~ multi_normal(mu, Sigma);
#     }
# }
# '

d <- list(nrow = nrow(dat),
                 ncol = ncol(dat),
                 n_comp = length(dat_complete),
                 n_miss = sum(is.na(dat)),
                 dat_complete = dat_complete,
                 ind_pres = ind_pres,
                 ind_miss = ind_miss)

# stan_program <- stan_program_direct
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
# out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = lapply(1:4, function(x) list(sigma = rep(3,n_dim))))
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, adapt_delta = 0.95, refresh = 10, init = 0.1, max_treedepth = 20)
# out <- mod$variational(data = d)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
par(mfrow = c(2,2), mar = c(4,4,1,4))
samps <- data.frame(as_draws_df(out$draws()))
corr_mats <- lapply(1:nrow(samps), function(ri) matrix(c(as.matrix(samps[ri,grep("corr_mat", colnames(samps))])), ncol = n_dim))
corr_mats_array <- do.call(abind::abind, c(corr_mats, list(along = 3)))
mean_corr_mat <- sapply(1:nrow(corr_mats_array), function(ri) sapply(1:ncol(corr_mats_array), function(ci) mean(corr_mats_array[ri,ci,])))

plot(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean), ylab = "posterior means for params", pch = 19, col = adjustcolor(1,0.4)); 
cor(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean))
cor(mu, apply(dat, 2, mean, na.rm = T))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
plot(sds, apply(samps[,grep("sigma", colnames(samps))], 2, mean), 
     ylab = "posterior means for params", pch = 19, col = adjustcolor(1,0.4)); 
cor(sds, apply(samps[,grep("sigma", colnames(samps))], 2, mean))
cor(sds, apply(dat, 2, sd, na.rm = T))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
plot(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)], xlim = c(-1,1), 
     ylim = c(-1,1), ylab = "posterior means for params", pch = 19, col = adjustcolor(1,0.4)); 
cor(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)])
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
plot(R[upper.tri(R)], cor(dat, use = "pair")[upper.tri(R)], xlim = c(-1,1), 
     ylim = c(-1,1), ylab = "pairwise sample Pearson correlation", pch = 19, col = adjustcolor(1,0.4)); 
cor(R[upper.tri(R)], cor(dat, use = "pair")[upper.tri(R)])
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)


#### implementing a basic multivariate ordinal probit ####

n_dim <- 10
mu <- rnorm(n_dim, 0, sd = 3)
R <- rlkj(n_dim, 1)
sds <- rexp(n_dim, rate = 0.5)
Sig <- diag(sds) %*% R %*% diag(sds)
N <- 50
dat <- mvrnorm(N, mu, Sig)

# specify and apply cutpoints
n_cut <- 1
cutpoints <- cumsum(c(0,rexp(n_cut-1, 1)))
# cutpoints <- c(0,1:(n_cut-1)/n_cut*2)
gthan <- lapply(1:n_cut, function(ci) dat > cutpoints[ci])
dat_discr <- gthan[[1]]
if(n_cut >= 2){
  for(i in 2:n_cut){dat_discr <- dat_discr + gthan[[i]]}  
} else {
  dat_discr <- 0 + dat_discr
}

# unobserve some fraction of the data
frac_missing <- 0.0
n_miss <- N*n_dim*frac_missing
miss <- sample.int(prod(dim(dat)), size = n_miss)
dat_discr[miss] <- NA

# get all the intermediate data
dat_inter <- dat_discr[!is.na(dat_discr) & (dat_discr != n_cut) & (dat_discr != 0)]

# find indices of the different data types
ind_inter <- which(!is.na(dat_discr) & (dat_discr != n_cut) & (dat_discr != 0), arr.ind = TRUE)
ind_min <- which(dat_discr == 0, arr.ind = TRUE)
ind_max <- which(dat_discr == n_cut, arr.ind = TRUE)
ind_miss <- which(is.na(dat_discr), arr.ind = TRUE)

# STAN model
stan_program <- '
data {
    int<lower=0> nrow; // number of rows in discrete data matrix (ie the # of individuals / observations)
    int<lower=0> ncol; // number of columns in discrete data matrix (ie the # of dimensions of the MVN)
    int<lower=0> n_cut; // number of cutpoints
    int<lower=0> n_inter; // Number of intermediate or (upper-and-lower-bounded) interval values
    int<lower=0> n_min; // Number of minimum values (so only bounded above)
    int<lower=0> n_max; // Number of maximum values (so only bounded below)
    int<lower=0> n_miss; // Number of missing values (so unbounded)
    array[n_inter] int<lower=0,upper=n_cut> dat_inter; // the intermediate states
    array[n_inter, 2] int ind_inter; // the indices of obs. in the data matrix w/ intermediate states
    array[n_miss, 2] int ind_miss; // the indices of missing obs. in the data matrix
    array[n_min, 2] int ind_min; // the indices of obs. in the data matrix w/ minimum states
    array[n_max, 2] int ind_max; // the indices of obs. in the data matrix w/ maximum states
}

parameters {
    // MVN distribution parameters
    vector[ncol] mu;
    cholesky_factor_corr[ncol] L;
    vector<lower=0>[ncol] sigma; 
    
    // cutpoints, constrained for identifiability
    positive_ordered[n_cut-1] later_cutpoints;
    
    // liabilities with different bounds
    vector<lower=0,upper=1>[n_inter] y_inter; // corresponds to proportion within interval between cutpoints
    vector<lower=0>[n_max] y_max; // corresponds to distance above the very last cutpoint
    vector<upper=0>[n_min] y_min; // first cutpoint is later constrained to 0
    vector[n_miss] ymiss;
}

transformed parameters {
    // construct full cutpoints vector, constraining location of first cutpoint to 0
    vector<lower=0>[n_cut] cutpoints = append_row(0, later_cutpoints);

    matrix[nrow, ncol] liabilities;
    // fill in liability matrix
    for(n in 1:n_inter) {
        liabilities[ind_inter[n,1],ind_inter[n,2]] = cutpoints[dat_inter[n]] + (cutpoints[dat_inter[n]+1] - cutpoints[dat_inter[n]]) * y_inter[n];
    }
    
    for(n in 1:n_min) {
        liabilities[ind_min[n,1],ind_min[n,2]] = y_min[n];
    }
    
    for(n in 1:n_max) {
        liabilities[ind_max[n,1],ind_max[n,2]] = cutpoints[n_cut] + y_max[n];
    }
    
    for(n in 1:n_miss){
        liabilities[ind_miss[n,1],ind_miss[n,2]] = ymiss[n];
    }
}

model {
    //priors & augmented likelihood
    sigma ~ exponential(0.5);
    L ~ lkj_corr_cholesky(2);
    for(i in 1:nrow){
      liabilities[i,] ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L));
    }
    
    //jacobian adjustment for non-linear cutpoint transformation
    for(n in 1:n_inter) {
      target += log(abs(cutpoints[dat_inter[n]+1] - cutpoints[dat_inter[n]]));
    }
}

generated quantities {
  matrix[ncol, ncol] corr_mat; 
  corr_mat = L * L\';
}
'

d <- list(nrow = nrow(dat_discr),
          ncol = ncol(dat_discr),
          n_inter = nrow(ind_inter),
          n_miss = nrow(ind_miss),
          n_max = nrow(ind_max),
          n_min = nrow(ind_min),
          dat_inter = dat_inter,
          ind_inter = ind_inter,
          ind_min = ind_min,
          ind_max = ind_max,
          ind_miss = ind_miss,
          n_cut = n_cut)

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, data = d, 
                  parallel_chains = 4, adapt_delta = 0.98, max_treedepth = 15, refresh = 10, init = 0.1)
# out <- mod$pathfinder(d, init = 0.1)
# out <- mod$variational(data = d, iter = 1E5)
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
# pairs(cbind(samps[,grep("mu", colnames(samps))], samps[,setdiff(grep("cutpoints", colnames(samps)), grep("later_cutpoints", colnames(samps)))]),
#       col = adjustcolor(1, 0.1), pch = 19)

#do some plotting plot
par(mfrow = c(2,2))
plot(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean), ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean)), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
legend(lty = 2, lwd = 2, col = adjustcolor(2,0.5), legend = "1-to-1 line", x = "topleft")
plot(sds, apply(samps[,grep("sigma", colnames(samps))], 2, mean), ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(sds, apply(samps[,grep("sigma", colnames(samps))], 2, mean)), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
plot(cutpoints, apply(samps[,setdiff(grep("cutpoints", colnames(samps)), grep("later_cutpoints", colnames(samps)))], 2, mean), ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(cutpoints, apply(samps[,setdiff(grep("cutpoints", colnames(samps)), grep("later_cutpoints", colnames(samps)))], 2, mean)), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)

corr_mats <- lapply(1:nrow(samps), function(ri) matrix(c(as.matrix(samps[ri,grep("corr_mat", colnames(samps))])), ncol = n_dim))
corr_mats_array <- do.call(abind::abind, c(corr_mats, list(along = 3)))
mean_corr_mat <- sapply(1:nrow(corr_mats_array), function(ri) sapply(1:ncol(corr_mats_array), function(ci) mean(corr_mats_array[ri,ci,])))
plot(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)], ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)]), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)

