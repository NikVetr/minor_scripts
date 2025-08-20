library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#simulate coefficients
p <- 1E2 #total number of coefficients
x <- rnorm(p) #sample true coefficients

#incorporate sparsity?
prop_null <- 0.8
if(prop_null != 0){
  p_null <- ceiling(p*prop_null)
  null_inds <- sample(1:p, p_null)
  x[null_inds] <- x[null_inds] / 10
}

#simulate data and fit OLS model
n <- 10 #define sample size across dimensions
e_sd_sd <- 1 #define differential power across two dimensions
e_sd <- rexp(p) * e_sd_sd #sample element-wise error
sim_and_fit_lm <- function(b, err_sd, nobs){
  asim <- rnorm(1)
  xsim <- rnorm(nobs)
  esim <- rnorm(n = nobs, mean = 0, sd = err_sd)
  ysim <- xsim * b + esim
  fit <- lm(ysim ~ xsim)
  summary(fit)$coefficients[2,1:2]
}

fits <- do.call(rbind, lapply(1:p, function(i){
  sim_and_fit_lm(b = x[i], err_sd = e_sd[i], nobs = n)
}))

x_err <- fits[,"Estimate"]
sd_x_err <- fits[,"Std. Error"]

model_loc <- "~/scripts/minor_scripts/postdoc/regularized_horseshoe_effect-recovery.stan"
# model_loc <- "~/scripts/minor_scripts/postdoc/regularized_horseshoe_effect-recovery_origparam.stan"
model_string <-  paste0(readLines(model_loc), collapse = "\n")

#pack data into a list
dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err, n=n, p0 = p - p_null)

#### fit model ####
mod <- cmdstan_model(model_loc)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  data = dat, adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 1, init = 0.1)

#### mcmc diagnostics ####

#check convergence
summ <- fit$summary()
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws()))
samps_x <- subset_samps("x", data.table(as_draws_df(fit$draws("x"))))
x_pm <- apply(samps_x, 2, mean)


par(mfrow = c(1,2), mar = c(5,5,5,1))
plot(x, x_err, xlab = "true coefficient", 
     ylab = latex2exp::TeX("\\textbf{OLS estimate} of coefficient"), 
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                  round(cor(cbind(c(x), c(x_err)))[1,2], 3))),
     col = adjustcolor(1, 0.5), pch = 19)
abline(0,1, col = 1, lty = 2, lwd = 2)

plot(x, x_pm, xlab = "true coefficient", 
     ylab = latex2exp::TeX("\\textbf{posterior mean} of coefficient"),
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                  round(cor(cbind(c(x), c(x_pm)))[1,2], 3))),
     col = adjustcolor(1, 0.5), pch = 19)
abline(0,1, col = 1, lty = 2, lwd = 2)
