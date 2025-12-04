set.seed(1)

library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
library(uniLasso)
source("~/repos/Stan2R/R/functions.R")

#generate high level params
lasso_col <- "orange"
refit_stan_model <- T
use_unilasso <- T
unilasso.lower.limits <- 0
ps_step <- F
bayes_model <- c("horseshoe", "spike_and_slab", "NEG")[1]
rop <- 0.8
d <- 200
invert_effects <- T
error_prop_var <- 1
rs <- c(1, rop^(1:d))
R <- outer(1:d, 1:d, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)

#sample inputs
n <- 300
x <- matrix(rnorm(n*d), n, d)
if(n>d){
  x <- x %*% solve(chol(cov(x))) %*% chol(R)   
} else {
  x <- x %*% chol(R) #whitening will just put 0s on the diag 
}

#invert sign of x
if(invert_effects){
  invertinds <- sample(1:d, floor(d/2), replace = F)
  x[,invertinds] <- -x[,invertinds]
}

#simulate outcomes w/ B = 1
nB <- 20
nzi <- sample(1:d, nB, replace = F)
y <- apply(x[,nzi], 1, sum)
e_sigma <- sd(y) * error_prop_var
e <- rnorm(n) * e_sigma
y <- y + e

#fit all params marginally
estimates <- do.call(rbind, apply(x, 2, function(xv) 
  summary(lm(y ~ xv))$coef[2,c(1,4)], simplify = F))
cexes <- -log10(estimates[,2])
cexes[cexes == Inf] <- max(cexes[cexes != Inf])
cexes <- 0.5 + ((cexes - min(cexes)) / max(cexes - min(cexes))) * 3



#fit spike and slab model in Stan

#pack data into a list
dat <- list(n=n, d=d, x=x, y=y)

#specify model with uncertainty
model_string_spike_and_slab_hack <-  "
data{
  int n;
  int d;
  matrix[n, d] x;
  vector[n] y;
}
parameters{
  real a;
  vector[d] b_raw;
  vector[d] logit_p;
  real<lower=0> sd_b;
  real<lower=0> sd_e;
  //real<lower=0> sd_p;
  real mu_b;
  real mu_p;
}
model{
  // define parameters without sparsity
  a ~ std_normal();
  b_raw ~ std_normal();
  logit_p ~ std_normal();
  sd_b ~ std_normal();
  sd_e ~ std_normal();
  //sd_p ~ std_normal();
  mu_b ~ std_normal();
  mu_p ~ std_normal();

  // define parameters with sparsity
  //vector[d] p = inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5);
  //vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5) - 0.5)));
  vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p) - 0.5)));
  vector[d] b = mu_b + b_raw .* p * sd_b;

  //likelihood
  y ~ normal(a + x * b, sd_e);
}
"

model_string_horseshoe <-  "
data{
  int n;
  int d;
  matrix[n, d] x;
  vector[n] y;
}
parameters{
  real a;
  vector[d] b_raw;
  vector[d] logit_p;
  real<lower=0> sd_b;
  real<lower=0> sd_e;
  //real<lower=0> sd_p;
  real mu_b;
  real mu_p;
}
model{
  // define parameters without sparsity
  a ~ std_normal();
  b_raw ~ std_normal();
  logit_p ~ std_normal();
  sd_b ~ std_normal();
  sd_e ~ std_normal();
  //sd_p ~ std_normal();
  mu_b ~ std_normal();
  mu_p ~ std_normal();

  // define parameters with sparsity
  //vector[d] p = inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5);
  //vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5) - 0.5)));
  vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p) - 0.5)));
  vector[d] b = mu_b + b_raw .* p * sd_b;

  //likelihood
  y ~ normal(a + x * b, sd_e);
}
"

model_string_NEG <-  "
data{
  int n;
  int d;
  matrix[n, d] x;
  vector[n] y;
}
parameters{
  real a;
  vector[d] b_raw;
  vector[d] logit_p;
  real<lower=0> sd_b;
  real<lower=0> sd_e;
  //real<lower=0> sd_p;
  real mu_b;
  real mu_p;
}
model{
  // define parameters without sparsity
  a ~ std_normal();
  b_raw ~ std_normal();
  logit_p ~ std_normal();
  sd_b ~ std_normal();
  sd_e ~ std_normal();
  //sd_p ~ std_normal();
  mu_b ~ std_normal();
  mu_p ~ std_normal();

  // define parameters with sparsity
  //vector[d] p = inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5);
  //vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p * sd_p * 5) - 0.5)));
  vector[d] p = 1 / (1 + exp(-100 * (inv_logit(mu_p * 2 - 2 + logit_p) - 0.5)));
  vector[d] b = mu_b + b_raw .* p * sd_b;

  //likelihood
  y ~ normal(a + x * b, sd_e);
}
"

model_string <- switch(EXPR = bayes_model,
                         horseshoe = model_string_horseshoe,
                         spike_and_slab = model_string_spike_and_slab_hack,
                         NEG = model_string_NEG
                       )

#separately estimate sd for center spike? and multiply together with b in likelihood, effectively non-centering
#can get that sd by taking a weighted average of the slab sd and 0, like in a horseshoe

#fit model
mod <- cmdstan_model(write_stan_file(model_string))

if(refit_stan_model | !exists("pp")){

fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  data = dat, adapt_delta = 0.85, parallel_chains = 4,
                  refresh = 100, max_treedepth = 15, 
                  thin = 1)

#check convergence
# summ <- fit$summary()
# print(summ[order(summ$ess_bulk),])
# print(summ[order(summ$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws()))

#compute posterior predictive
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}
npreds <- 500
sample_indices <- sample(1:nrow(samps), npreds, replace = F)
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

pp <- do.call(rbind, lapply(ppfits, function(x) x$p))
pb <- do.call(rbind, lapply(ppfits, function(x) x$b))

}


#fit lasso or unilasso regression
# if(use_unilasso){
#   cv_uni   <- cv.uniLasso(x, y, nfolds = 10)     # cross-validated path
#   lambda   <- cv_uni$lambda.min                 # same idea as before
#   uni_fit  <- cv_uni$glmnet.fit                 # the internal glmnet object
#   est_coefs <- as.numeric(coef(cv_uni, s = "lambda.min"))[-1]
#   nzi_est <- which(est_coefs != 0)
#   alt_psf_fit <- summary(lm(y ~ x[, nzi_est]))$coef[-1, c(1, 4)]
# } else {
#   #fit lasso regression to see if we can recover nzi
#   lambda <- glmnet::cv.glmnet(x, y)$lambda.min
#   lasso_fit <- glmnet::glmnet(x, y, alpha = 1, 
#                               lambda = lambda)
#   est_coefs <- coef(lasso_fit)[-1]
#   nzi_est <- which(abs(est_coefs) > 1E-6)
#   # ps_fit <- suppressWarnings(selectiveInference::fixedLassoInf(x, y, 
#   #                         beta = est_coefs, 
#   #                         lambda = lambda))
#   # ps_fit <- cbind(ps_fit$coef0, ps_fit$pv)
#   alt_psf_fit <- summary(lm(y ~ x[,nzi_est]))$coef[-1,c(1,4)]
# }

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


#then finds odds ratio
# ctab2x2 <- matrix(c(length(intersect(nzi, nzi_est_sub)),
#                     length(intersect(nzi, setdiff(1:d, nzi_est_sub))),
#                     length(intersect(setdiff(1:d, nzi), nzi_est_sub)),
#                     length(intersect(setdiff(1:d, nzi), setdiff(1:d, nzi_est_sub)))), 2, 2)
# rownames(ctab2x2) <- c("estimated positive", "estimated negative")
# colnames(ctab2x2) <- c("positive", "negative")
# 
# exact_test_out <- fisher.test(ctab2x2)
# ctab2x2
# exact_test_out$estimate


#do some plotting
cairo_pdf("~/sparsity_comparison.pdf", width = 800/72, height = 800/72)
par(mar = c(5,5,6,12), mfrow = c(2,1))
plot(estimates[,1], cex = cexes, xlab = "index", ylab = "estimate",
     ylim = c(-max(abs(estimates[,1])), max(abs(estimates[,1]))) * 1.1,
     pch = 19, col = adjustcolor(1,0.15))
abline(v = nzi, col = 2, lty = 2, lwd = 2)
# text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
#      labels = paste0("coef #", nzi, " (", c(-1,1)[nzi %in% invertinds + 1], ")"), 
#      xpd = NA, srt = 60, pos = 4, col = 2)
text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
     labels = paste0("coef #", nzi, " (1)"), 
     xpd = NA, srt = 60, pos = 4, col = 2)
text(x = par("usr")[2], y = 1, 
     labels = "true effect", 
     xpd = NA, pos = 4, col = 2)
text(x = par("usr")[2], y = par("usr")[4] + diff(par("usr")[3:4])/60, 
     labels = "all unlabeled true effects = 0", 
     xpd = NA, pos = 4, col = 1)
text(x = par("usr")[2], y = par("usr")[4] - diff(par("usr")[3:4])/30, 
     labels = latex2exp::TeX("|corr(x$_i$, x$_j$)| = 0.9^|i-j|"), 
     xpd = NA, pos = 4, col = 1)

abline(h = 0, col = 1)
ifelse2 <- function(test, yes, no) if(test) return(yes) else return(no)
legend(xpd = NA, x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/2, 
       legend = c("marginal est effect", ifelse2(use_unilasso, "uni-lasso coef", "lasso coef"), "bayes Pr(p > 0.5)", 
                  "true +", "false +", "circle diam. ∝ -log(p-val)", ifelse2(use_unilasso, "(no p-vals from uni-lasso)", NULL)),
       pch = c(19, 19, 18, 1, 1, NA, ifelse2(use_unilasso, NA, NULL)), 
       col = c(adjustcolor(c(1,lasso_col,4), 0.5), 1, 2, NA, ifelse2(use_unilasso, NA, NULL)), pt.cex = 1.5, bty = "n")

#numbers are > 1 (true beta vals) bc correlated predictors all influence the outcome


## indices of the variables we are about to plot
plot_idx <- nzi_est                     # same as before
## add points
points(plot_idx, coef_plot,
       pch = 19, col = adjustcolor(lasso_col, 0.5), cex = ps_cex)

## outline true vs false positives
for (s in c(1, 0.95, 0.9, 0.85)) {
  tp <- plot_idx %in% nzi           # nzi = true non-zero indices
  points(plot_idx[ tp], coef_plot[ tp], pch = 1, col = 1, cex = ps_cex[ tp] * s)
  points(plot_idx[!tp], coef_plot[!tp], pch = 1, col = 2, cex = ps_cex[!tp] * s)
}


#plot top hits from bayesian analysis
mean_mix_params <- apply(pp, 2, mean)
mean_coefs <- apply(pb, 2, mean)
mean_coefs_present <- apply(pb * pp, 2, function(x) mean(x[x>1E-3]))
# nzi_est_sub_bayes <- which(rank(mean_mix_params) > (d-nB+1E-6))
nzi_est_sub_bayes <- which(mean_mix_params > 0.5)
points(nzi_est_sub_bayes, mean_coefs_present[nzi_est_sub_bayes], 
       pch = 18, col = adjustcolor(4, 0.5), cex = 2)
for(cex_mult in c(1, 0.95, 0.9, 0.85, 0.8, 0.75)){
  points(nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi], mean_coefs_present[nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi]], 
         pch = 5, col = 1, cex = 2 * cex_mult)
  points(nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)], mean_coefs_present[nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)]], 
         pch = 5, col = 2, cex = 2 * cex_mult)
}

dev.off()
