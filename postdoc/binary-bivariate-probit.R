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

# simulate data
D <- 10
mu <- rnorm(D, 0, sd = 1)
R <- rlkj(D, 1)
sds <- rep(1,D)
Sig <- diag(sds) %*% R %*% diag(sds)
N <- 500
liab <- mvrnorm(N, mu, Sig)

# specify and apply cutpoints
y <- (liab > 0) + 0.0

# construct data object
dat <- list(D = D,
            N = N,
            y = y)

# STAN model
stan_loc <- "~/scripts/minor_scripts/postdoc/binary-bivariate-probit.stan"
stan_program <- paste0(readLines(stan_loc), collapse = "\n")
mod <- cmdstan_model(stan_loc)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, data = dat, 
                  parallel_chains = 4, adapt_delta = 0.95, max_treedepth = 10, refresh = 100, init = 0.1)
out_pf <- mod$pathfinder(history_size=100, max_lbfgs_iters=100,
                      data = dat, num_threads = 2, init = 0.2)

#mcmc diagnostics
summ <- out$summary()
summ[order(summ$ess_bulk),]

#check inference
samps <- data.frame(as_draws_df(out$draws()))
# pairs(cbind(samps[,grep("mu", colnames(samps))], samps[,setdiff(grep("cutpoints", colnames(samps)), grep("later_cutpoints", colnames(samps)))]),
#       col = adjustcolor(1, 0.1), pch = 19)

#do some plotting plot
par(mfrow = c(1,2))
plot(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean), ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(mu, apply(samps[,grep("mu", colnames(samps))], 2, mean)), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
legend(lty = 2, lwd = 2, col = adjustcolor(2,0.5), legend = "1-to-1 line", x = "topleft")

corr_mats <- lapply(1:nrow(samps), function(ri) matrix(c(as.matrix(samps[ri,grep("R", colnames(samps))])), ncol = n_dim))
corr_mats_array <- do.call(abind::abind, c(corr_mats, list(along = 3)))
mean_corr_mat <- sapply(1:nrow(corr_mats_array), function(ri) sapply(1:ncol(corr_mats_array), function(ci) mean(corr_mats_array[ri,ci,])))
plot(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)], ylab = "posterior means for params",
     main = paste0("Pearson Correlation = ", round(cor(R[upper.tri(R)], mean_corr_mat[upper.tri(mean_corr_mat)]), 2)))
abline(0,1, col = adjustcolor(2,0.5), lty = 2, lwd = 2)
