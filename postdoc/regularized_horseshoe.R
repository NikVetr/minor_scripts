# set.seed(1)

library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
# source("~/repos/Stan2R/R/functions.R")

#### simulate data  ####

#specify high level params
bad_tau0 <- F
lasso_col <- "orange"
refit_stan_model <- T
use_unilasso <- T
unilasso.lower.limits <- 0
use_centered <- F
use_s2z <- F
ps_step <- F
rop <- 0.8
p <- 200
prop_n0_b <- 0.05
posterior_inclusion_threshold <- 0.8
B_pow <- 0
Bn_mult <- 0
invert_effects <- F
error_prop_var <- 1
rs <- c(1, rop^(1:p))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)

#sample inputs
n <- 100
x <- matrix(rnorm(n*p), n, p)

#discretize?
invlogit <- function(x) exp(x) / (1 + exp(x))
x <- t(apply(x, 1, function(ri) rbinom(p, 2, prob = invlogit(ri))))
if(n>p){
  x <- x %*% solve(chol(cov(x))) %*% chol(R)   
} else {
  x <- x %*% chol(R) #whitening will just put 0s on the diag 
}

#invert sign of x
if(invert_effects){
  invertinds <- sample(1:p, floor(p/2), replace = F)
  x[,invertinds] <- -x[,invertinds]
}

#simulate outcomes
nB <- ceiling(prop_n0_b * p)
nzi <- sample(1:p, nB, replace = F)
nnzi <- setdiff(1:p, nzi)
B <- qnorm(p = rbeta(nB, 0.2, 0.2))^B_pow
Bn <- qnorm(p = rbeta(p-nB, 50, 50))*Bn_mult
y <- apply(x[,nzi] %*% diag(B), 1, sum) + 
  apply(x[,nnzi] %*% diag(Bn), 1, sum)
e_sigma <- sd(y) * error_prop_var
e <- rnorm(n) * e_sigma
y <- y + e


#### marginal fits ####

#fit all params marginally
estimates <- do.call(rbind, apply(x, 2, function(xv) 
  summary(lm(y ~ xv))$coef[2,c(1,4)], simplify = F))

#### lasso fits ####


if (use_unilasso) {
  library(uniLasso)
  cv_fit   <- cv.uniLasso(x, y, nfolds = 10, lower.limits = unilasso.lower.limits)
  est_coef <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]   # γ̂’s
} else {
  lambda   <- glmnet::cv.glmnet(x, y)$lambda.min
  fit      <- glmnet::glmnet(x, y, lambda = lambda, alpha = 1)
  est_coef <- as.numeric(coef(fit))[-1]
}

nzi_est <- which(abs(est_coef) > 1E-5)      # selected variables

if (ps_step && length(nzi_est) > 0) {
  alt_psf_fit <- summary(lm(y ~ x[, nzi_est]))$coef[-1, c(1, 4)]
  coef_plot   <- alt_psf_fit[, 1]                 # OLS β̂
  p_plot      <- alt_psf_fit[, 2]                 # p-values
} else {
  coef_plot   <- est_coef[nzi_est]                # raw γ̂ (UniLasso) or β̂ (Lasso)
  p_plot      <- rep(NA_real_, length(coef_plot)) # no p-values
}

#recover point sizes
if (ps_step) {
  # size ∝ –log10(p)
  ps_cex <- -log10(p_plot)
  ps_cex[is.infinite(ps_cex)] <- max(ps_cex[is.finite(ps_cex)])
  ps_cex <- 1.5 + (ps_cex - min(ps_cex)) / (max(ps_cex) - min(ps_cex) + 1e-6) * 3
} else {
  # size ∝ |coefficient|
  mag     <- abs(coef_plot)
  mag     <- abs(coef_plot) * 0 + 1
  ps_cex  <- 1.5 + (mag - min(mag)) / (max(mag) - min(mag) + 1e-6) * 3
}

# #fit lasso regression to see if we can recover nzi
# lambda <- glmnet::cv.glmnet(x, y)$lambda.min
# lasso_fit <- glmnet::glmnet(x, y, alpha = 1, 
#                             lambda = lambda)
# est_coefs <- coef(lasso_fit)[-1]
# nzi_est <- which(abs(est_coefs) > 1E-6)
# # ps_fit <- suppressWarnings(selectiveInference::fixedLassoInf(x, y, 
# #                         beta = est_coefs, 
# #                         lambda = lambda))
# # ps_fit <- cbind(ps_fit$coef0, ps_fit$pv)
# alt_psf_fit <- summary(lm(y ~ x[,nzi_est]))$coef[-1,c(1,4)]

#### horseshoe fits  ####

#pack data into a list

# regularized horseshoe prior params
p0 <- nB #approx number of nonzero params
pseudosigma <- sqrt(abs(1/mean(y)/(1-mean(y))))
tau0 <- p0/(p-p0)*pseudosigma/sqrt(n)
if(bad_tau0){
  tau0 <- tau0 * 3
}

#generate Q
ZerosumBasis <- function(n) {
  v <- rep(1, n)
  M <- diag(n) - (1/n) * (v %*% t(v))
  Q_full <- qr.Q(qr(M), complete = FALSE)
  Q_part <- Q_full[, 1:(n - 1)]
  Q_part
}
Q <- ZerosumBasis(p)

## data
dat <- list(
  n = n,
  p = p,
  y = y,
  x = x,
  scale_icept = 5,
  scale_global = tau0,
  nu_global = 1,
  nu_local = 1,
  slab_scale = 2,
  slab_df = 100,
  Q = Q
)

#fit model

if(refit_stan_model | !exists("lmean")){

#compile all models
mod_ns2z <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe.stan")
mod_nc_ns2z <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe_nc.stan")
mod_s2z <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe_s2z.stan")
mod_nc_s2z <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe_nc_s2z.stan")

#accommodate sum to zero toggle
if(use_s2z){
  mod <- mod_s2z
  mod_nc <- mod_nc_s2z
} else {
  mod <- mod_ns2z
  mod_nc <- mod_nc_ns2z
}

#fit centered or non-centered parameterization
if(use_centered){
  hs_fit <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2,
                        data = dat, adapt_delta = 0.85, parallel_chains = 4,
                        refresh = 100, max_treedepth = 15, 
                        thin = 1)
} else {
  hs_fit <- mod_nc$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2,
                        data = dat, adapt_delta = 0.85, parallel_chains = 4,
                        refresh = 100, max_treedepth = 15, 
                        thin = 1)
}

#check convergence
check_mcmc_diag <- T
if(check_mcmc_diag){
  summ <- hs_fit$summary()
  print(summ[order(summ$ess_bulk),])
  print(summ[order(summ$rhat, decreasing = T),])
}


#### inspect posterior  ####

#extract samples and inspect
samps <- as_draws_df(hs_fit$draws())
bsamps <- munge_samps("beta", subset_samps("beta", data.table(samps)))
b <- do.call(rbind, bsamps)
lsamps <- munge_samps("lambda", subset_samps("lambda", data.table(samps)))
l <- do.call(rbind, lsamps)
ltsamps <- munge_samps("lambda_tilde", subset_samps("lambda_tilde", data.table(samps)))
lt <- do.call(rbind, ltsamps)

#process these
not0 <- function(x) max(mean(x > 0), 1 - mean(x > 0))
bprob <- apply(b, 2, not0)
bmean <- apply(b, 2, mean)
lprob <- apply(l, 2, not0)
lmean <- apply(l, 2, mean)
ltprob <- apply(lt, 2, not0)
ltmean <- apply(lt, 2, mean)

}

hist(bmean, breaks = 100)
points(bmean[nzi], rep(0, nB), pch = 18, col = 3, cex = 2)
hist(bprob, breaks = 100)
points(bprob[nzi], rep(0, nB), pch = 18, col = 3, cex = 2)
plot(bmean, ltmean)
points(bmean[nzi], ltmean[nzi], pch = 19, col = 2, cex = 2)
text(bmean[nzi], ltmean[nzi], labels = 1:nB, col = "white", cex = 0.75)

#compute posterior predictive
mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}
npreds <- 10
sample_indices <- sample(1:nrow(samps), npreds, replace = F)
model_string <-  paste0(mod_nc$code(), collapse = "\n")
ppfits <- parallel::mclapply(seq_along(sample_indices), function(i){
  
  mcprint(paste0(i, " "))
  sample_index <- sample_indices[i]
  r_code <- parse_Stan(stan_code = strsplit(model_string, "\n")[[1]], 
                       dat = dat, 
                       samps = as.data.table(samps), 
                       output_file = NA, 
                       sample_index = sample_index, 
                       post_pred_sim = T, 
                       sim = F)
  # open_script(r_code)
  
  stan_env <- new.env()
  eval(parse(text = r_code), envir = stan_env)
  posterior_predictive_mass <- stan_env$out
  
  return(posterior_predictive_mass)
}, mc.cores = 8)


#### plotting ####

#first add to plot
cexes <- -log10(estimates[,2])
cexes[cexes == Inf] <- max(cexes[cexes != Inf])
cexes <- 0.5 + ((cexes - min(cexes)) / max(cexes - min(cexes))) * 3
par(mar = c(5,5,6,14))
plot.new()
plot.window(ylim = c(-max(abs(estimates[,1])), max(abs(estimates[,1]))) * 1.1, xlim = c(0,p))
axis(1); mtext("parameter index", 1, line = 3, cex = 1.5)
axis(2); mtext("estimate", 2, line = 3, cex = 1.5)
abline(v = nzi, col = 2, lty = 2, lwd = 2)
points(estimates[,1], cex = cexes, xlab = "index", ylab = "estimate",
     ylim = c(-max(abs(estimates[,1])), max(abs(estimates[,1]))) * 1.1,
     pch = 19, col = adjustcolor(1,0.5))
# text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
#      labels = paste0("coef #", nzi, " (", c(-1,1)[nzi %in% invertinds + 1], ")"), 
#      xpd = NA, srt = 60, pos = 4, col = 2)
text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
     labels = paste0("coef #", nzi, " (", round(B, 2),")"), 
     xpd = NA, srt = 60, pos = 4, col = 2)
# text(x = par("usr")[2], y = 1, 
#      labels = "true effect", 
#      xpd = NA, pos = 4, col = 2)


abline(h = 0, col = 1)
# legend(xpd = NA, x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/2, 
#        legend = c("marginal est effect", 
#        paste0(ifelse(use_unilasso, "unilasso", "lasso"), 
#               " estimate ", 
#               ifelse(ps_step, "(w/ post-selection)", "")), 
#               "Bayes P((X>0)⊕(X<0))>0.95", "true +", "false +", "circle diam. ∝ -log(p-val)"),
#        pch = c(19, 19, 18, 1, 1, NA), col = c(adjustcolor(c(1,lasso_col,4), 0.5), 1, 2, NA), pt.cex = 1.5, bty = "n")

ifelse2 <- function(test, yes, no) if(test) return(yes) else return(no)
legend(xpd = NA, x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/2, 
       legend = c("marginal est effect", ifelse2(use_unilasso, "uni-lasso coef", "lasso coef"), 
                  paste0("Bayes P((X>0)⊕(X<0))>", posterior_inclusion_threshold), 
                  "true +", "false +", "circle diam. ∝ -log(p-val)", ifelse2(use_unilasso, "(no p-vals from uni-lasso)", NULL)),
       pch = c(19, 19, 18, 1, 1, NA, ifelse2(use_unilasso, NA, NULL)), 
       col = c(adjustcolor(c(1,lasso_col,4), 0.5), 1, 2, NA, ifelse2(use_unilasso, NA, NULL)), pt.cex = 1.5, bty = "n")

#numbers are > 1 (true beta vals) bc correlated predictors all influence the outcome

#plot the lowest pvals?
# lowest_pvals <- which(rank(alt_psf_fit[,2]) < (nB + 1E-6))
# lowest_pvals <- 1:nrow(alt_psf_fit)
# nzi_est_sub <- nzi_est[lowest_pvals]
# ps_cexes <- -log10(alt_psf_fit[lowest_pvals,2])
# ps_cexes[ps_cexes == Inf] <- max(ps_cexes[ps_cexes != Inf])
# ps_cexes <- 1.5 + ((ps_cexes - min(ps_cexes)) / max(ps_cexes - min(ps_cexes))) * 3
# 
# points(nzi_est_sub, alt_psf_fit[lowest_pvals,1], 
#        pch = 19, col = adjustcolor(lasso_col, 0.5), cex = ps_cexes)
# for(cex_mult in c(1, 0.95, 0.9, 0.85)){
#   points(nzi_est_sub[nzi_est_sub %in% nzi], alt_psf_fit[lowest_pvals,1][nzi_est_sub %in% nzi], 
#          pch = 1, col = 1, cex = ps_cexes[nzi_est_sub %in% nzi] * cex_mult)
#   points(nzi_est_sub[!(nzi_est_sub %in% nzi)], alt_psf_fit[lowest_pvals,1][!(nzi_est_sub %in% nzi)], 
#          pch = 1, col = 2, cex = ps_cexes[!(nzi_est_sub %in% nzi)] * cex_mult)
# }

#plot the (uni)lasso estimates
plot_idx <- nzi_est
## add points
points(plot_idx, coef_plot,
       pch = 19, col = adjustcolor(lasso_col, 0.5), cex = ps_cex)
for (s in c(1, 0.95, 0.9, 0.85)) {
  tp <- plot_idx %in% nzi           # nzi = true non-zero indices
  points(plot_idx[ tp], coef_plot[ tp], pch = 1, col = 1, cex = ps_cex[ tp] * s)
  points(plot_idx[!tp], coef_plot[!tp], pch = 1, col = 2, cex = ps_cex[!tp] * s)
}


#plot top hits
mean_mix_params <- bprob
mean_coefs <- bmean
mean_coefs_present <- bprob
# nzi_est_sub_bayes <- which(rank(mean_mix_params) > (d-nB+1E-6))
# nzi_est_sub_bayes <- which(mean_mix_params > 0.9)
# nzi_est_sub_bayes <- order(mean_mix_params, decreasing = T)[1:nB]
nzi_est_sub_bayes <- which(mean_mix_params > posterior_inclusion_threshold)

points(nzi_est_sub_bayes, mean_coefs_present[nzi_est_sub_bayes], 
       pch = 18, col = adjustcolor(4, 0.5), cex = 2)
for(cex_mult in c(1, 0.95, 0.9, 0.85, 0.8, 0.75)){
  points(nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi], mean_coefs_present[nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi]], 
         pch = 5, col = 1, cex = 2 * cex_mult)
  points(nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)], mean_coefs_present[nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)]], 
         pch = 5, col = 2, cex = 2 * cex_mult)
}

if(abs(B_pow) < 1E-4){
  text(x = par("usr")[2], y = 1, 
       labels = "true effect", 
       xpd = NA, pos = 4, col = 2)
  abline(h = 1, col = 2)
}

if(abs(Bn_mult) < 1E-4){
  text(x = par("usr")[2], y = par("usr")[4] + diff(par("usr")[3:4])/60, 
       labels = "all unlabeled true effects = 0", 
       xpd = NA, pos = 4, col = 1)  
}
text(x = par("usr")[2], y = par("usr")[4] - diff(par("usr")[3:4])/30, 
     labels = latex2exp::TeX(paste0("|corr(x$_i$, x$_j$)| = ", rop, "^|i-j|")), 
     xpd = NA, pos = 4, col = 1)

#then finds odds ratio

#first for bayesian analysis
ctab2x2_bayes <- matrix(c(length(intersect(nzi, nzi_est_sub_bayes)),
                    length(intersect(nzi, setdiff(1:p, nzi_est_sub_bayes))),
                    length(intersect(setdiff(1:p, nzi), nzi_est_sub_bayes)),
                    length(intersect(setdiff(1:p, nzi), setdiff(1:p, nzi_est_sub_bayes)))), 2, 2)
rownames(ctab2x2_bayes) <- c("estimated positive", "estimated negative")
colnames(ctab2x2_bayes) <- c("positive", "negative")

exact_test_out_bayes <- fisher.test(ctab2x2_bayes)
ctab2x2_bayes
exact_test_out_bayes$estimate


#then for lasso
ctab2x2_lasso <- matrix(c(length(intersect(nzi, nzi_est)),
                    length(intersect(nzi, setdiff(1:p, nzi_est))),
                    length(intersect(setdiff(1:p, nzi), nzi_est)),
                    length(intersect(setdiff(1:p, nzi), setdiff(1:p, nzi_est)))), 2, 2)
rownames(ctab2x2_lasso) <- c("estimated positive", "estimated negative")
colnames(ctab2x2_lasso) <- c("positive", "negative")

exact_test_out_lasso <- fisher.test(ctab2x2_lasso)
ctab2x2_lasso
exact_test_out_lasso$estimate

#plot Bayesian OR as a function of probability threshold?

# Set up the range for posterior_inclusion_threshold values
posterior_inclusion_thresholds <- seq(0.5, 1, by = 0.01)
odds_ratios <- numeric(length(posterior_inclusion_thresholds))
sample_sizes <- numeric(length(posterior_inclusion_thresholds))

# Calculate odds ratios and sample sizes for each threshold
for (i in seq_along(posterior_inclusion_thresholds)) {
  threshold <- posterior_inclusion_thresholds[i]
  nzi_est_sub_bayes <- which(mean_mix_params > threshold)
  ctab2x2_bayes <- matrix(c(
    length(intersect(nzi, nzi_est_sub_bayes)),  # True Positives
    length(intersect(nzi, setdiff(1:p, nzi_est_sub_bayes))),  # False Negatives
    length(intersect(setdiff(1:p, nzi), nzi_est_sub_bayes)),  # False Positives
    length(intersect(setdiff(1:p, nzi), setdiff(1:p, nzi_est_sub_bayes)))  # True Negatives
  ), 2, 2)
  rownames(ctab2x2_bayes) <- c("estimated positive", "estimated negative")
  colnames(ctab2x2_bayes) <- c("positive", "negative")
  
  # Perform Fisher's exact test to get the odds ratio
  exact_test_out_bayes <- fisher.test(ctab2x2_bayes)
  odds_ratios[i] <- exact_test_out_bayes$estimate
  
  # Store the sample size (number of estimated positives)
  sample_sizes[i] <- length(nzi_est_sub_bayes)
}

# Plot the odds ratios against the threshold values
par(mar = c(5,5,3,5))
plot(posterior_inclusion_thresholds, odds_ratios, type = "l", col = "darkblue", lwd = 2,
     xlab = "Posterior Inclusion Threshold", ylab = "Odds Ratio (TP / FN) / (FP / TN)",
     main = "", log = "y")
abline(h = exact_test_out_lasso$estimate, col = "darkred", lwd = 2)
top_or <- ceiling(max(odds_ratios[odds_ratios!=Inf]))
axis(2, at = 1:top_or, tick = T, labels = rep("", top_or), , 
     lwd.ticks = 0.5, tcl=-0.4)

# Add the sample size on the right margin
par(new = TRUE)
plot(posterior_inclusion_thresholds, sample_sizes, type = "l", col = "skyblue", lwd = 2,
     axes = FALSE, xlab = "", ylab = "", log = "y")
axis(4, col = 1, col.axis = 1)
mtext("Number of Positive IDs", side = 4, line = 3, col = 1)
abline(h = exact_test_out_lasso$estimate, col = "pink", lwd = 2)
ss_ord_mag <- ceiling(log10(ceiling(max(sample_sizes))))
ss_ticks <- c(outer(1:9, 10^(0:ss_ord_mag)))
axis(4, at = ss_ticks, tick = T, labels = rep("", length(ss_ticks)), 
     lwd.ticks = 0.5, tcl=-0.4)


legend("topright", legend = c("Horseshoe OR", "uniLasso OR", "Horseshoe # Hits", "uniLasso # Hits"),
       col = c("darkblue", "darkred", "skyblue", "pink"), lwd = 2, cex = 0.8)


#compare MSE

#simulate new data
beta_true <- numeric(p)
beta_true[nzi]  <- B          # non-zero effects
beta_true[nnzi] <- Bn
n_test <- 1000                                # pick any size
x_test <- matrix(rnorm(n_test * p), n_test, p)
x_test <- t(apply(x_test, 1, function(ri)
  rbinom(p, 2, prob = invlogit(ri))))
if (n > p) {
  x_test <- x_test %*% solve(chol(cov(x_test))) %*% chol(R)
} else {
  x_test <- x_test %*% chol(R)
}
if (invert_effects) x_test[, invertinds] <- -x_test[, invertinds]

# generate y_test with the same SNR recipe
signal_test <- x_test %*% beta_true
e_sigma_test <- sd(signal_test) * error_prop_var
y_test <- signal_test + rnorm(n_test, 0, e_sigma_test)

#make predictions

#uni(lasso)
if (use_unilasso) {
  intercept_hat <- as.numeric(coef(cv_fit,  s = "lambda.min"))[1]
} else {
  intercept_hat <- as.numeric(coef(fit))[1]  # glmnet fit
}
y_hat_lasso <- drop(intercept_hat + x_test %*% est_coef)
mse_lasso <- mean( (y_test - y_hat_lasso   )^2 )
  
#horseshoe
hs_draws <- as_draws_df(hs_fit$draws(c("beta0", "beta", "sigma")))
beta0_draws <- hs_draws$beta0
beta_draws  <- as.matrix(hs_draws[, grep("^beta\\[", names(hs_draws))])
sigma_draws <- hs_draws$sigma
y_hat_hs <- rowMeans(matrix(beta0_draws, nrow = n_test, ncol = nrow(beta_draws), byrow = TRUE) + 
                       x_test %*% t(beta_draws))
mse_hs    <- mean( (y_test - y_hat_hs   )^2 )

#prediction metrics
r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum( (y - mean(y))^2 )

cat(paste0("\n===  Prediction performance on fresh test data  ===\n",
sprintf("(Uni)-Lasso  :  MSE = %.3f   |  R^2 = %.3f |  r = %.3f\n",
            mse_lasso, r2(y_test, y_hat_lasso), cor(y_test, y_hat_lasso)),
sprintf("Horseshoe+   :  MSE = %.3f   |  R^2 = %.3f |  r = %.3f\n",
            mse_hs, r2(y_test, y_hat_hs), cor(y_test, y_hat_hs))
))
