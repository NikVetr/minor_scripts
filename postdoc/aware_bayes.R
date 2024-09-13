#exploring extensions to naive bayes under feature interdependence
library(MASS)
library(data.table)
library(mvtnorm)
source("~/repos/Stan2R/R/functions.R")

#### functions ####
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

#### simulate data ####

test_methods <- function(){

#mvn corrmat
p <- 15
n <- 50
# Rs <- replicate(100, rlkj(p), simplify = F)
# R <- Rs[[which.min(sapply(Rs, det))]]
R <- diag(p) #meet naive bayes assumption with identity matrix
r <- 0.5
R <- diag(p) + r - diag(p) * r #worst case scenario

#get lower cholesky factor for MVN sampling
L <- t(chol(R))

#mvn means
diff_x_mu <- L %*% rnorm(p)
x1_mu <- diff_x_mu/2
x2_mu <- -diff_x_mu/2

#individual liabilies (train)
x1 <- t(replicate(n, c(L %*% rnorm(p) + x1_mu)))
x2 <- t(replicate(n, c(L %*% rnorm(p) + x2_mu)))

#individual liabilies (test)
ntest <- 500
x1t <- t(replicate(ntest, c(L %*% rnorm(p) + x1_mu)))
x2t <- t(replicate(ntest, c(L %*% rnorm(p) + x2_mu)))

#individual binaries (train)
x1b <- (x1 > 0) * 1
x2b <- (x2 > 0) * 1

#individual binaries (test)
x1tb <- (x1t > 0) * 1
x2tb <- (x2t > 0) * 1

#population sample frequencies
sf1 <- apply(rbind(x1b, 0, 1), 2, mean)
sf2 <- apply(rbind(x2b, 0, 1), 2, mean)

#### naive bayes ####

#fit naive bayes and apply to test
prob_x1 <- 0.5
prob_x2 <- 1 - prob_x1

unnorm_prob_x1t_1 <- apply(t(t(x1tb) * sf1) + t(t(1-x1tb) * (1-sf1)), 1, prod) * prob_x1
unnorm_prob_x1t_2 <- apply(t(t(x1tb) * sf2) + t(t(1-x1tb) * (1-sf2)), 1, prod) * prob_x2
prob_x1t_1 <- unnorm_prob_x1t_1 / (unnorm_prob_x1t_1 + unnorm_prob_x1t_2)

unnorm_prob_x2t_1 <- apply(t(t(x2tb) * sf1) + t(t(1-x2tb) * (1-sf1)), 1, prod) * prob_x1
unnorm_prob_x2t_2 <- apply(t(t(x2tb) * sf2) + t(t(1-x2tb) * (1-sf2)), 1, prod) * prob_x2
prob_x2t_1 <- unnorm_prob_x2t_1 / (unnorm_prob_x2t_1 + unnorm_prob_x2t_2)

#binarize with flat loss
test_results_NB <- matrix(NA, 2, 2, dimnames = list(c("A", "B"), c("est_A", "est_B")))
test_results_NB[1,1] <- mean(prob_x1t_1 > 0.5)
test_results_NB[1,2] <- mean(prob_x1t_1 < 0.5)
test_results_NB[2,1] <- mean(prob_x2t_1 > 0.5)
test_results_NB[2,2] <- mean(prob_x2t_1 < 0.5)
test_results_NB

#### try PCA-NB-Normal ####

# lda_out <- MASS::lda(x = rbind(x1b, x2b), grouping = rep(1:2, each = n))
prop_var <- 0.9

#get training set PCs
within_vcov <- (cov(x1b) + cov(x2b)) / 2
total_vcov <- cov(rbind(x1b, x2b))
cov2use <- total_vcov
eigendecomp <- eigen(cov2use)
U <- eigendecomp$values
nPCs <- max(min(sum(cumsum(U) / sum(U) < prop_var) + 1, p), 1)
V <- eigendecomp$vectors
x1_PCs <- t(t((x1b %*% V %*% diag(1/sqrt(U)))[,1:nPCs]))
x2_PCs <- t(t((x2b %*% V %*% diag(1/sqrt(U)))[,1:nPCs]))

#estimate PC distributions (basic version)
x1_PCs_mu <- apply(x1_PCs, 2, mean)
x1_PCs_sigma <- apply(x1_PCs, 2, sd)
x2_PCs_mu <- apply(x2_PCs, 2, mean)
x2_PCs_sigma <- apply(x2_PCs, 2, sd)

#fit naive bayes on test set
x1t_PCs <- t(t((x1tb %*% V %*% diag(1/sqrt(U)))[,1:nPCs]))
x2t_PCs <- t(t((x2tb %*% V %*% diag(1/sqrt(U)))[,1:nPCs]))

unnorm_logprob_x1t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(logprob_1 = sum(dnorm(x = x1t_PCs[i,], mean = x1_PCs_mu, sd = x1_PCs_sigma, log = T)) + log(prob_x1),
           logprob_2 = sum(dnorm(x = x1t_PCs[i,], mean = x2_PCs_mu, sd = x2_PCs_sigma, log = T)) + log(prob_x2))
  out <- out - max(out)
  return(out)
})))
prob_x1t_1 <- exp(unnorm_logprob_x1t$logprob_1) / apply(exp(unnorm_logprob_x1t), 1, sum)

unnorm_logprob_x2t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(logprob_1 = sum(dnorm(x = x2t_PCs[i,], mean = x1_PCs_mu, sd = x1_PCs_sigma, log = T)) + log(prob_x1),
           logprob_2 = sum(dnorm(x = x2t_PCs[i,], mean = x2_PCs_mu, sd = x2_PCs_sigma, log = T)) + log(prob_x2))
  out <- out - max(out)
  return(out)
})))
prob_x2t_1 <- exp(unnorm_logprob_x2t$logprob_1) / apply(exp(unnorm_logprob_x2t), 1, sum)

#binarize with flat loss
test_results_NB_PCA <- matrix(NA, 2, 2, dimnames = list(c("A", "B"), c("est_A", "est_B")))
test_results_NB_PCA[1,1] <- mean(prob_x1t_1 > 0.5)
test_results_NB_PCA[1,2] <- mean(prob_x1t_1 < 0.5)
test_results_NB_PCA[2,1] <- mean(prob_x2t_1 > 0.5)
test_results_NB_PCA[2,2] <- mean(prob_x2t_1 < 0.5)

#### try PCA-NB-t ####

#estimate PC distributions (student's t version)
ll_t <- function(params, data) {
  m <- params[1] # location (mean)
  s <- params[2] # scale (standard deviation)
  df <- params[3] # degrees of freedom
  df_penalty <- dnorm(df - 2, sd = 20, log = T)
  -sum(dt((data - m) / s, df = df, log = TRUE) - log(s)) - df_penalty 
}

x1_PCs_t.params <- data.frame(t(apply(x1_PCs, 2, function(d){
  init <- c(mean(d), sd(d), 5) # mean, sd, and initial df
  fit <- optim(init, ll_t, data = d, method = "L-BFGS-B", lower = c(-Inf, 0.0001, 2))$par  
  return(c(m = fit[1], s = fit[2], df = fit[3]))
})))

x2_PCs_t.params <- data.frame(t(apply(x2_PCs, 2, function(d){
  init <- c(mean(d), sd(d), 5) # mean, sd, and initial df
  fit <- optim(init, ll_t, data = d, method = "L-BFGS-B", lower = c(-Inf, 0.0001, 2))$par  
  return(c(m = fit[1], s = fit[2], df = fit[3]))
})))

#fit naive bayes on test set
unnorm_logprob_x1t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(logprob_1 = sum(dt(x = (x1t_PCs[i,] - x1_PCs_t.params$m) / x1_PCs_t.params$s, 
                              df = x1_PCs_t.params$df, log = T)) + log(prob_x1),
           logprob_2 = sum(dt(x = (x1t_PCs[i,] - x2_PCs_t.params$m) / x2_PCs_t.params$s, 
                              df = x2_PCs_t.params$df, log = T)) + log(prob_x2))
  
  #alternatively
  out <- c(logprob_1 = sum(dt(x = (x1t_PCs[i,] - x1_PCs_mu) / x1_PCs_sigma, 
                              df = n-2, log = T)) + log(prob_x1),
           logprob_2 = sum(dt(x = (x1t_PCs[i,] - x2_PCs_mu) / x2_PCs_sigma, 
                              df = n-2, log = T)) + log(prob_x2))
  
  out <- out - max(out)
  return(out)
})))
prob_x1t_1 <- exp(unnorm_logprob_x1t$logprob_1) / apply(exp(unnorm_logprob_x1t), 1, sum)

unnorm_logprob_x2t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(logprob_1 = sum(dt(x = (x2t_PCs[i,] - x1_PCs_t.params$m) / x1_PCs_t.params$s, 
                              df = x1_PCs_t.params$df, log = T)) + log(prob_x1),
           logprob_2 = sum(dt(x = (x2t_PCs[i,] - x2_PCs_t.params$m) / x2_PCs_t.params$s, 
                              df = x2_PCs_t.params$df, log = T)) + log(prob_x2))
  
  #alternatively
  out <- c(logprob_1 = sum(dt(x = (x2t_PCs[i,] - x1_PCs_mu) / x1_PCs_sigma, 
                              df = n-2, log = T)) + log(prob_x1),
           logprob_2 = sum(dt(x = (x2t_PCs[i,] - x2_PCs_mu) / x2_PCs_sigma, 
                              df = n-2, log = T)) + log(prob_x2))
  
  out <- out - max(out)
  return(out)
})))
prob_x2t_1 <- exp(unnorm_logprob_x2t$logprob_1) / apply(exp(unnorm_logprob_x2t), 1, sum)

#binarize with flat loss
test_results_NB_PCA_t <- matrix(NA, 2, 2, dimnames = list(c("A", "B"), c("est_A", "est_B")))
test_results_NB_PCA_t[1,1] <- mean(prob_x1t_1 > 0.5)
test_results_NB_PCA_t[1,2] <- mean(prob_x1t_1 < 0.5)
test_results_NB_PCA_t[2,1] <- mean(prob_x2t_1 > 0.5)
test_results_NB_PCA_t[2,2] <- mean(prob_x2t_1 < 0.5)

#### try PCA-NB-KDE ####

#estimate PC distributions (KDE version)
range_test.x1t <- apply(x1t_PCs, 2, range, simplify = F)
range_test.x2t <- apply(x2t_PCs, 2, range, simplify = F)
x1_PCs_KDE <- lapply(1:ncol(x1_PCs), function(i){
  d <- x1_PCs[,i]
  kde <- density(d, from = range_test.x1t[[i]][1], to = range_test.x1t[[i]][2])
  kde_norm <- median(kde$y) / 10 #normalization constant
  return(list(x = kde$x, y = kde$y + kde_norm))
})

x2_PCs_KDE <- lapply(1:ncol(x2_PCs), function(i){
  d <- x2_PCs[,i]
  kde <- density(d, from = range_test.x2t[[i]][1], to = range_test.x2t[[i]][2])
  kde_norm <- median(kde$y) / 10 #normalization constant
  return(list(x = kde$x, y = kde$y + kde_norm))
})

#fit naive bayes on test set
unnorm_logprob_x1t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  d <- x1t_PCs[i,]
  out <- c(
    logprob_1 = sum(sapply(1:length(d), function(j){
      log(x1_PCs_KDE[[j]]$y[which.min(abs(x1_PCs_KDE[[j]]$x - d[j]))])
    })) + log(prob_x1),
    logprob_2 = sum(sapply(1:length(d), function(j){
      log(x2_PCs_KDE[[j]]$y[which.min(abs(x2_PCs_KDE[[j]]$x - d[j]))])
    })) + log(prob_x2)
  )
  out <- out - max(out)
  return(out)
})))
prob_x1t_1 <- exp(unnorm_logprob_x1t$logprob_1) / apply(exp(unnorm_logprob_x1t), 1, sum)

unnorm_logprob_x2t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  d <- x2t_PCs[i,]
  out <- c(
    logprob_1 = sum(sapply(1:length(d), function(j){
      log(x1_PCs_KDE[[j]]$y[which.min(abs(x1_PCs_KDE[[j]]$x - d[j]))])
    })) + log(prob_x1),
    logprob_2 = sum(sapply(1:length(d), function(j){
      log(x2_PCs_KDE[[j]]$y[which.min(abs(x2_PCs_KDE[[j]]$x - d[j]))])
    })) + log(prob_x2)
  )
  out <- out - max(out)
  return(out)
})))
prob_x2t_1 <- exp(unnorm_logprob_x2t$logprob_1) / apply(exp(unnorm_logprob_x2t), 1, sum)

#binarize with flat loss
test_results_NB_PCA_KDE <- matrix(NA, 2, 2, dimnames = list(c("A", "B"), c("est_A", "est_B")))
test_results_NB_PCA_KDE[1,1] <- mean(prob_x1t_1 > 0.5)
test_results_NB_PCA_KDE[1,2] <- mean(prob_x1t_1 < 0.5)
test_results_NB_PCA_KDE[2,1] <- mean(prob_x2t_1 > 0.5)
test_results_NB_PCA_KDE[2,2] <- mean(prob_x2t_1 < 0.5)


#### try multi-probit ####
integrate_over_joint_posterior <- F

# construct data objects
dat_1 <- list(D = p,
              N = n,
              y = x1b)

dat_2 <- list(D = p,
              N = n,
              y = x2b)

# read in the model
stan_loc <- "~/scripts/minor_scripts/postdoc/binary-bivariate-probit.stan"
stan_program <- paste0(readLines(stan_loc), collapse = "\n")
mod <- suppressMessages(cmdstan_model(stan_loc))

#fit model
sink(tempfile())
out_1 <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat_1, 
                    parallel_chains = 4, adapt_delta = 0.95, max_treedepth = 10, 
                    refresh = 100, init = 0.1, show_exceptions = F, show_messages = F)
out_2 <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat_2, 
                    parallel_chains = 4, adapt_delta = 0.95, max_treedepth = 10, 
                    refresh = 100, init = 0.1, show_exceptions = F, show_messages = F)
sink()

#can inspect mcmc diagnostics later
# summ_1 <- out_1$summary()
# summ_1[order(summ_1$ess_bulk),]
# summ_2 <- out_2$summary()
# summ_2[order(summ_2$ess_bulk),]

#extract posterior draws
samps_1 <- data.table(as_draws_df(out_1$draws()))
samps_2 <- data.table(as_draws_df(out_2$draws()))

mu_1 <- data.frame(subset_samps("mu", samps_1))
mu_2 <- data.frame(subset_samps("mu", samps_2))

R_1 <- abind::abind(munge_samps("R", subset_samps("R", samps_1)), along = 3)
R_2 <- abind::abind(munge_samps("R", subset_samps("R", samps_2)), along = 3)

#easy first pass, extract posterior means
mu_mu_1 <- apply(mu_1, 2, mean)
mu_mu_2 <- apply(mu_2, 2, mean)
mu_R_1 <- apply(R_1, c(1,2), mean)
mu_R_2 <- apply(R_2, c(1,2), mean)

#compute orthant probabilities
x1tb_uniq <- x1tb[!duplicated(x1tb),]
dupe_inds_x1tb <- match(apply(x1tb, 1, paste0, collapse = "-"), apply(x1tb_uniq, 1, paste0, collapse = "-"))
x2tb_uniq <- x2tb[!duplicated(x2tb),]
dupe_inds_x2tb <- match(apply(x2tb, 1, paste0, collapse = "-"), apply(x2tb_uniq, 1, paste0, collapse = "-"))

traits_to_bounds <- function(xb, thresh = 0){
  xb_ub <- xb
  xb_ub[xb == 1] <- Inf
  xb_ub[xb == 0] <- thresh
  xb_lb <- xb
  xb_lb[xb == 1] <- thresh
  xb_lb[xb == 0] <- -Inf
  return(list(
    upper = xb_ub,
    lower = xb_lb
  ))
}

x1tb_bounds <- traits_to_bounds(x1tb_uniq)
x2tb_bounds <- traits_to_bounds(x2tb_uniq)

op_1t <- do.call(rbind, lapply(1:nrow(x1tb_uniq), function(i){
  if(integrate_over_joint_posterior){
    in_1 <- sapply(1:nrow(samps_1), function(j){
      pmvnorm(lower = x1tb_bounds$lower[i,], 
              upper = x1tb_bounds$upper[i,], 
              mean = mu_1[j,], corr = R_1[,,j], algorithm = GenzBretz())[[1]]
    })
    in_2 <- sapply(1:nrow(samps_1), function(j){
      pmvnorm(lower = x1tb_bounds$lower[i,], 
              upper = x1tb_bounds$upper[i,], 
              mean = mu_2[j,], corr = R_2[,,j], algorithm = GenzBretz())[[1]]
    })
  } else {
    in_1 <- pmvnorm(lower = x1tb_bounds$lower[i,], 
                    upper = x1tb_bounds$upper[i,], 
                    mean = mu_mu_1, corr = mu_R_1, algorithm = GenzBretz())[[1]]
    in_2 <- pmvnorm(lower = x1tb_bounds$lower[i,], 
                    upper = x1tb_bounds$upper[i,], 
                    mean = mu_mu_2, corr = mu_R_2, algorithm = GenzBretz())[[1]]
  }
  return(c(in_1 = mean(in_1), in_2 = mean(in_2)))
}))[dupe_inds_x1tb,]

op_2t <- do.call(rbind, lapply(1:nrow(x2tb_uniq), function(i){
  if(integrate_over_joint_posterior){
    in_1 <- sapply(1:nrow(samps_1), function(j){
      pmvnorm(lower = x2tb_bounds$lower[i,], 
              upper = x2tb_bounds$upper[i,], 
              mean = mu_1[j,], corr = R_1[,,j], algorithm = GenzBretz())[[1]]
    })
    in_2 <- sapply(1:nrow(samps_1), function(j){
      pmvnorm(lower = x2tb_bounds$lower[i,], 
              upper = x2tb_bounds$upper[i,], 
              mean = mu_2[j,], corr = R_2[,,j], algorithm = GenzBretz())[[1]]
    })
  } else {
    in_1 <- pmvnorm(lower = x2tb_bounds$lower[i,], 
                    upper = x2tb_bounds$upper[i,], 
                    mean = mu_mu_1, corr = mu_R_1, algorithm = GenzBretz())[[1]]
    in_2 <- pmvnorm(lower = x2tb_bounds$lower[i,], 
                    upper = x2tb_bounds$upper[i,], 
                    mean = mu_mu_2, corr = mu_R_2, algorithm = GenzBretz())[[1]]
  }
  return(c(in_1 = mean(in_1), in_2 = mean(in_2)))
}))[dupe_inds_x2tb,]

#fit bayes on test set
unnorm_logprob_x1t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(
    logprob_1 = log(op_1t[i,1]) + log(prob_x1),
    logprob_2 = log(op_1t[i,2]) + log(prob_x2)
  )
  out <- out - max(out)
  return(out)
})))
prob_x1t_1 <- exp(unnorm_logprob_x1t$logprob_1) / apply(exp(unnorm_logprob_x1t), 1, sum)

unnorm_logprob_x2t <- data.frame(do.call(rbind, lapply(1:ntest, function(i){
  out <- c(
    logprob_1 = log(op_2t[i,1]) + log(prob_x1),
    logprob_2 = log(op_2t[i,2]) + log(prob_x2)
  )
  out <- out - max(out)
  return(out)
})))
prob_x2t_1 <- exp(unnorm_logprob_x2t$logprob_1) / apply(exp(unnorm_logprob_x2t), 1, sum)

#binarize with flat loss
test_results_mvProbit <- matrix(NA, 2, 2, dimnames = list(c("A", "B"), c("est_A", "est_B")))
test_results_mvProbit[1,1] <- mean(prob_x1t_1 > 0.5)
test_results_mvProbit[1,2] <- mean(prob_x1t_1 < 0.5)
test_results_mvProbit[2,1] <- mean(prob_x2t_1 > 0.5)
test_results_mvProbit[2,2] <- mean(prob_x2t_1 < 0.5)

#### results ####
list(test_results_NB = test_results_NB,
     test_results_NB_PCA = test_results_NB_PCA,
     test_results_NB_PCA_t = test_results_NB_PCA_t,
     test_results_NB_PCA_KDE = test_results_NB_PCA_KDE,
     test_results_mvProbit = test_results_mvProbit)

}

n_sim_reps <- 100
n_results <- 5
results <- lapply(1:n_sim_reps, function(i){print(i); test_methods()})
results <- lapply(1:n_results, function(i) lapply(results, function(res) res[[i]]))
results_v <- lapply(1:n_results, function(i) unlist(lapply(results[[i]], function(res) diag(res))))
names(results_v) <- c("Naive Bayes", 
                      "Aware Bayes (PCA, Gaussian PCs)", 
                      "Aware Bayes (PCA, t PCs)", 
                      "Aware Bayes (PCA, KDE PCs)",
                      "Aware Bayes (multivariate Probit)")
means <- round(sapply(results_v, mean), 3)

par(mfrow = c(n_results,1))
for(i in 1:n_results){
  hist(results_v[[i]], breaks = 0:20/20, main = names(results_v)[i], 
       xlab = "Proportion Correct in Both Groups")
  abline(v = means[i], col = 2, lwd = 3)
  text(x = means[i], y = par("usr")[4], pos = 3, col = 2, 
       labels = paste0("mean = ", means[i]), xpd = NA)
}
