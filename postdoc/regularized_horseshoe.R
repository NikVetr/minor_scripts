
library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#### simulate data  ####

#generate high level params
rop <- 0.95
p <- 200
invert_effects <- T
error_prop_var <- 1
rs <- c(1, rop^(1:p))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)

#sample inputs
n <- 50
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
nB <- 10
nzi <- sample(1:p, nB, replace = F)
nnzi <- setdiff(1:p, nzi)
B <- qnorm(p = rbeta(nB, 0.2, 0.2))
Bn <- qnorm(p = rbeta(p-nB, 50, 50))
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

#fit lasso regression to see if we can recover nzi
lambda <- glmnet::cv.glmnet(x, y)$lambda.min
lasso_fit <- glmnet::glmnet(x, y, alpha = 1, 
                            lambda = lambda)
est_coefs <- coef(lasso_fit)[-1]
nzi_est <- which(abs(est_coefs) > 1E-6)
# ps_fit <- suppressWarnings(selectiveInference::fixedLassoInf(x, y, 
#                         beta = est_coefs, 
#                         lambda = lambda))
# ps_fit <- cbind(ps_fit$coef0, ps_fit$pv)
alt_psf_fit <- summary(lm(y ~ x[,nzi_est]))$coef[-1,c(1,4)]

#### horseshoe fits  ####

#pack data into a list

# regularized horseshoe prior params
p0 <- nB #approx number of nonzero params
pseudosigma <- sqrt(abs(1/mean(y)/(1-mean(y))))
tau0 <- p0/(p-p0)*pseudosigma/sqrt(n)

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
  slab_df = 100
)

#fit model
mod <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe.stan")
mod_alt <- cmdstan_model("~/scripts/minor_scripts/postdoc/regularized_horseshoe_alt.stan")
fit <- mod_alt$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2,
                  data = dat, adapt_delta = 0.85, parallel_chains = 4,
                  refresh = 100, max_treedepth = 15, 
                  thin = 1)

#check convergence
# summ <- fit$summary()
# print(summ[order(summ$ess_bulk),])
# print(summ[order(summ$rhat, decreasing = T),])


#### inspect posterior  ####

#extract samples and inspect
samps <- as_draws_df(fit$draws())
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
model_string <-  paste0(mod_alt$code(), collapse = "\n")
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
plot(estimates[,1], cex = cexes, xlab = "index", ylab = "estimate",
     ylim = c(-max(abs(estimates[,1])), max(abs(estimates[,1]))) * 1.1,
     pch = 19, col = adjustcolor(1,0.5))
abline(v = nzi, col = 2, lty = 2, lwd = 2)
# text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
#      labels = paste0("coef #", nzi, " (", c(-1,1)[nzi %in% invertinds + 1], ")"), 
#      xpd = NA, srt = 60, pos = 4, col = 2)
text(x = nzi - 1, y = par("usr")[4] + diff(par("usr")[3:4])/100, 
     labels = paste0("coef #", nzi, " (", round(B, 2),")"), 
     xpd = NA, srt = 60, pos = 4, col = 2)
# text(x = par("usr")[2], y = 1, 
#      labels = "true effect", 
#      xpd = NA, pos = 4, col = 2)
text(x = par("usr")[2], y = par("usr")[4] + diff(par("usr")[3:4])/60, 
     labels = "all unlabeled true effects = 0", 
     xpd = NA, pos = 2, col = 1)
text(x = par("usr")[2], y = par("usr")[4] - diff(par("usr")[3:4])/30, 
     labels = latex2exp::TeX(paste0("|corr(x$_i$, x$_j$)| = ", rop, "^|i-j|")), 
     xpd = NA, pos = 2, col = 1)


abline(h = 0, col = 1)
legend(xpd = NA, x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/2, 
       legend = c("marginal est effect", "laplace MAP w/ PS (lasso)", "Bayes P((X>0)⊕(X<0))>0.95", "true +", "false +", "circle diam. ∝ -log(p-val)"),
       pch = c(19, 19, 18, 1, 1, NA), col = c(adjustcolor(c(1,2,4), 0.5), 1, 2, NA), pt.cex = 1.5, bty = "n")

#numbers are > 1 (true beta vals) bc correlated predictors all influence the outcome


lowest_pvals <- which(rank(alt_psf_fit[,2]) < (nB + 1E-6))
lowest_pvals <- 1:nrow(alt_psf_fit)
nzi_est_sub <- nzi_est[lowest_pvals]
ps_cexes <- -log10(alt_psf_fit[lowest_pvals,2])
ps_cexes[ps_cexes == Inf] <- max(ps_cexes[ps_cexes != Inf])
ps_cexes <- 1.5 + ((ps_cexes - min(ps_cexes)) / max(ps_cexes - min(ps_cexes))) * 3
points(nzi_est_sub, alt_psf_fit[lowest_pvals,1], 
       pch = 19, col = adjustcolor(3, 0.5), cex = ps_cexes)

for(cex_mult in c(1, 0.95, 0.9, 0.85)){
  points(nzi_est_sub[nzi_est_sub %in% nzi], alt_psf_fit[lowest_pvals,1][nzi_est_sub %in% nzi], 
         pch = 1, col = 1, cex = ps_cexes[nzi_est_sub %in% nzi] * cex_mult)
  points(nzi_est_sub[!(nzi_est_sub %in% nzi)], alt_psf_fit[lowest_pvals,1][!(nzi_est_sub %in% nzi)], 
         pch = 1, col = 2, cex = ps_cexes[!(nzi_est_sub %in% nzi)] * cex_mult)
}

#plot top hits
mean_mix_params <- bprob
mean_coefs <- bmean
mean_coefs_present <- bprob
# nzi_est_sub_bayes <- which(rank(mean_mix_params) > (d-nB+1E-6))
# nzi_est_sub_bayes <- which(mean_mix_params > 0.9)
nzi_est_sub_bayes <- order(mean_mix_params, decreasing = T)[1:nB]
points(nzi_est_sub_bayes, mean_coefs_present[nzi_est_sub_bayes], 
       pch = 18, col = adjustcolor(4, 0.5), cex = 2)
for(cex_mult in c(1, 0.95, 0.9, 0.85, 0.8, 0.75)){
  points(nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi], mean_coefs_present[nzi_est_sub_bayes[nzi_est_sub_bayes %in% nzi]], 
         pch = 5, col = 1, cex = 2 * cex_mult)
  points(nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)], mean_coefs_present[nzi_est_sub_bayes[!(nzi_est_sub_bayes %in% nzi)]], 
         pch = 5, col = 2, cex = 2 * cex_mult)
}

#then finds odds ratio
ctab2x2 <- matrix(c(length(intersect(nzi, nzi_est_sub)),
                    length(intersect(nzi, setdiff(1:p, nzi_est_sub))),
                    length(intersect(setdiff(1:p, nzi), nzi_est_sub)),
                    length(intersect(setdiff(1:p, nzi), setdiff(1:p, nzi_est_sub)))), 2, 2)
rownames(ctab2x2) <- c("estimated positive", "estimated negative")
colnames(ctab2x2) <- c("positive", "negative")

exact_test_out <- fisher.test(ctab2x2)
ctab2x2
exact_test_out$estimate
