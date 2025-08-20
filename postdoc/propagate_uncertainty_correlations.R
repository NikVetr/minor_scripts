library(limma)
library(wCorr)
library(ashr)
library(cmdstanr)
library(posterior)
library(statmod)

par_base <- par(no.readonly = TRUE)

#set a seed for reproducibility
# set.seed(1)

#### functions ####

#getting mixture weights
make_mix <- function(nu) {
  s2 <- c(0.29, 1.0, 12.68)
  if (nu <= 4) {
    # two-component match (variance only)
    w3 <- 0
    w1 <- (nu/(nu-2) - s2[2]) / (s2[1] - s2[2])
    w2 <- 1 - w1
  } else {
    w3 <- (nu/(nu-2) * s2[1]^2 - 
             3*nu^2/((nu-2)*(nu-4))*s2[1]) /
      ((s2[3]-s2[1]) * (s2[3]-s2[2]))
    w1 <- (nu/(nu-2) * (s2[3]^2 - s2[2]^2) - 
             3*nu^2/((nu-2)*(nu-4)) * (s2[3] - s2[2])) /
      ((s2[1]-s2[2]) * (s2[1]-s2[3]))
    w2 <- 1 - w1 - w3
  }
  list(w = c(w1, w2, w3),
       s2_extra = s2)   # add these to sd_x_err^2 inside Stan
}

# example for ν = 8
make_mix(8)

# heteroskedastic error-in-variables correlation estimator
hetero_cor <- function(x, y,
                       var_x,           # vector of measurement-error variances σ²_{x,i}
                       var_y,           # vector of measurement-error variances σ²_{y,i}
                       cov_xy = NULL,   # optional vector of error covariances τ_i (default 0)
                       ci = FALSE,      # logical: return 95 % bootstrap CI?
                       n_boot = 999,    # bootstrap replicates
                       seed = NULL) {
  
  # basic checks
  if (!(length(x) == length(y) &&
        length(x) == length(var_x) &&
        length(x) == length(var_y))) {
    stop("x, y, var_x, and var_y must have the same length")
  }
  n <- length(x)
  if (n < 3) stop("need at least 3 observations")
  
  # allow for optional correlated errors
  if (is.null(cov_xy)) cov_xy <- rep(0, n)
  if (length(cov_xy) != n) stop("cov_xy must have length n or be NULL")
  
  # sample moments (unbiased, denominator n-1)
  Sxx <- stats::var(x)
  Syy <- stats::var(y)
  Sxy <- stats::cov(x, y)
  
  # average measurement-error terms
  bar_sig2_x <- mean(var_x)
  bar_sig2_y <- mean(var_y)
  bar_tau     <- mean(cov_xy)
  
  # debias
  tilde_Sxx <- Sxx - bar_sig2_x
  tilde_Syy <- Syy - bar_sig2_y
  tilde_Sxy <- Sxy - bar_tau
  
  # guard against negative debiased variances
  if (tilde_Sxx <= 0 || tilde_Syy <= 0){
    # stop("debiased variance ≤ 0; measurement noise overwhelms the signal")
    rho_hat <- cor(x, y)
    boot_vals <- replicate(n_boot, {
      idx <- sample.int(n, replace = TRUE)
      cor(x[idx], y[idx])
    })
    ci_bounds <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
    out <- list(rho = rho_hat,
                ci95 = as.numeric(ci_bounds),
                boot_values = boot_vals)
    
    return(out)
  }
  
  rho_hat <- tilde_Sxy / sqrt(tilde_Sxx * tilde_Syy)
  
  # optional bootstrap CI
  if (!ci) return(rho_hat)
  
  if (!is.null(seed)) set.seed(seed)
  boot_vals <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    Sxx_b <- stats::var(x[idx])
    Syy_b <- stats::var(y[idx])
    Sxy_b <- stats::cov(x[idx], y[idx])
    rho_b <- (Sxy_b - mean(cov_xy[idx])) /
      sqrt((Sxx_b - mean(var_x[idx])) *
             (Syy_b - mean(var_y[idx])))
    rho_b
  })
  boot_vals <- boot_vals[!is.na(boot_vals)]
  
  ci_bounds <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
  
  list(rho = rho_hat,
       ci95 = as.numeric(ci_bounds),
       boot_values = boot_vals)
}

ccc <- function(x, y) {
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  
  #components
  mean_x <- mean(x)
  mean_y <- mean(y)
  var_x <- var(x)
  var_y <- var(y)
  cov_xy <- cov(x, y)
  
  #ccc formula
  numerator <- 2 * cov_xy
  denominator <- var_x + var_y + (mean_x - mean_y)^2
  
  rho_c <- numerator / denominator
  return(rho_c)
}

#### simulation parameters ####

use_ols_fit <- T
check_diag <- T

#for ols fit
n <- c(10, 5) #define sample size across dimensions
e_sd_sd <- c(4, 8) #define differential power across two dimensions? or just get from sample size
use_shrinkage_estimator <- F

#for static error fit
err_var <- 1
# err_var <- 16 * matrix(exp(rnorm(p*2)), ncol = 2, nrow = p)
error_r <- 0
error_R <- diag(2) * (1-error_r) + error_r
x_var <- 1
plot_scatter <- F

#for all fits
p <- p_base <- 5E2 #total number of samples
bivariate_laplace = T #simulate true coefficients from bivariate laplace?

#try pseudoreplication experiment?
pseudoreplicate <- T
pseudorep_times <- 10
pseudorep_sd <- 0.2
if(pseudoreplicate){
  p <- p_base * pseudorep_times
} else {
  p <- p_base
}

#specify model with uncertainty
model_index <- 15
model_dir <- "~/scripts/minor_scripts/postdoc/"
model_name <- list("correlation_uncertainty.stan", #1
                   "correlation_uncertainty_horseshoe-unpooled.stan", #2
                   "correlation_uncertainty_horseshoe-partially-pooled.stan", #3
                   "correlation_uncertainty_horseshoe-pooled.stan", #4
                   "correlation_uncertainty_logit-rescale.stan", #5
                   "correlation_uncertainty_bivariate-laplace.stan", #6
                   "bivariate_correlation_uncertainty_fast.stan", #7
                   "bivariate_correlation_uncertainty_fast_centered.stan", #8
                   "bivariate_correlation_uncertainty_fast_marginalize-out-latents.stan", #9
                   "bivariate_correlation_uncertainty_fast_scaled.stan", #10
                   "bivariate-t_correlation_uncertainty_fast.stan", #11
                   "bivariate_correlation_uncertainty_fast_chi2error.stan", #12
                   "bivariate_correlation_uncertainty_fast_marginalize-out-latents_vectorized.stan", #13
                   "bivariate-t_correlation_uncertainty_fast_marginalize-out-latents.stan", #14
                   "bivariate_var-infl-normal_correlation_uncertainty_fast_marginalize-out-latents_vectorized.stan", #15
                   "mixture-approx_biv-t_vectorized.stan" #16
)[[model_index]]
model_path <- paste0(model_dir, model_name)

#can maybe add a further term to the model to represent distribution of estimated errors?

#run multiple times?
nruns <- 100
true_corrs <- numeric(nruns)
sample_corrs <- numeric(nruns)
sample_corrs_noerror <- numeric(nruns)
posterior_corrs <- as.list(true_corrs)
summs <- as.list(true_corrs)
hetero_corrs <- numeric(nruns)
hetero_corrs_bs <- as.list(hetero_corrs)
w_corrs <- numeric(nruns)
ashr_corrs <- numeric(nruns)
w_corrs_bs <- as.list(hetero_corrs)
ashr_corrs_bs <- as.list(hetero_corrs)

for(run_i in 1:nruns){
  
  #print progress  
  cat(paste0("(", run_i, ") "))
  
  #### simulate input ####
  
  #simulate coefficients
  r <- -0.7 #true correlation between coefficients
  r <- runif(1,-1,1)
  R <- diag(2) * (1-r) + r #corresponding correlation matrix
  inflate_reported_error_by <- matrix(1/2, ncol = 2, nrow = p)
  inflate_reported_error_by <- matrix(exp(rnorm(p*2)/10), ncol = 2, nrow = p)
  inflate_reported_error_by <- matrix(1 + runif(p*2, -1, 1) * 0, ncol = 2, nrow = p)
  x <- matrix(rnorm(p*2) * sqrt(x_var), ncol = 2) %*% chol(R) #sample true coefficients
  
  #modify data according to model
  if(grepl("horseshoe", model_path)){
    prop_null <- 0
    if(prop_null != 0){
      # null_inds <- sample(1:(2*n), ceiling(n*prop_null*2))
      # x[null_inds] <- x[null_inds] / 10
      
      n_null <- ceiling(p*prop_null)
      null_inds <- sample(1:p, n_null)
      which_null <- 1 + rbinom(n_null, 2, prob = 0.5)
      x[null_inds[which_null == 2],] <- x[null_inds[which_null == 2],] / 10
      x[null_inds[which_null == 1], 1] <- x[null_inds[which_null == 1], 1] / 10
      x[null_inds[which_null == 3], 2] <- x[null_inds[which_null == 3], 2] / 10
    }
  }
  
  if(bivariate_laplace){
    exp_rv <- rexp(p, 1)
    x[,1] <- x[,1] * exp_rv
    x[,2] <- x[,2] * exp_rv
  }
  
  if(use_ols_fit){
    #simulate data and fit OLS model
    e_sd <- matrix(rexp(p*2), ncol = 2) %*% diag(e_sd_sd) #sample element-wise error
    sim_and_fit_lm <- function(b, err_sd, nobs){
      asim <- rnorm(1)
      xsim <- rnorm(nobs)
      esim <- rnorm(n = nobs, mean = 0, sd = err_sd)
      ysim <- asim + xsim * b + esim
      fit <- lm(ysim ~ xsim)
      summary(fit)$coefficients[2,1:2]
    }
    
    sim_data <- function(b, err_sd, nobs){
      asim <- rnorm(1)
      xsim <- rnorm(nobs)
      esim <- rnorm(n = nobs, mean = 0, sd = err_sd)
      ysim <- asim + xsim * b + esim
      return(list(
        xsim = xsim,
        ysim = ysim
      ))
    }
    
    fit_lm <- function(simdat){
      fit <- lm(simdat$ysim ~ simdat$xsim)
      summary(fit)$coefficients[2,1:2]
    }
      
    fits_1 <- do.call(rbind, lapply(1:p_base, function(i){
      if(!pseudoreplicate){
        sim_and_fit_lm(b = x[i,1], err_sd = e_sd[i,1], nobs = n[1])  
      } else {
        base_simdat <- sim_data(b = x[i,1], err_sd = e_sd[i,1], nobs = n[1])
        do.call(rbind, lapply(1:pseudorep_times, function(j){
          pr_simdat <- base_simdat
          pr_simdat$x + rnorm(n[1]) * pseudorep_sd
          pr_simdat$y + rnorm(n[1]) * pseudorep_sd
          fit_lm(pr_simdat)
        }))
      }
    }))
    
    fits_2 <- do.call(rbind, lapply(1:p_base, function(i){
      if(!pseudoreplicate){
        sim_and_fit_lm(b = x[i,2], err_sd = e_sd[i,2], nobs = n[2])
      } else {
        base_simdat <- sim_data(b = x[i,2], err_sd = e_sd[i,2], nobs = n[2])
        do.call(rbind, lapply(1:pseudorep_times, function(j){
          pr_simdat <- base_simdat
          pr_simdat$x + rnorm(n[2]) * pseudorep_sd
          pr_simdat$y + rnorm(n[2]) * pseudorep_sd
          fit_lm(pr_simdat)
        }))
      }
    }))
    
    x_err <- cbind(fits_1[,"Estimate"], fits_2[,"Estimate"])
    sd_x_err <- cbind(fits_1[,"Std. Error"], fits_2[,"Std. Error"])
    df <- cbind(rep(n[1] - 2, p), rep(n[2] - 2, p))
    
    if(use_shrinkage_estimator){
      v      <- as.numeric(sd_x_err^2)
      df_obs <- as.numeric(df)
      
      sq      <- squeezeVar(v, df = df_obs)
      sd_post <- matrix(sqrt(sq$var.post), nrow = nrow(sd_x_err))
      df_post <- matrix(sq$df.prior + df_obs, ncol = 2)
      
      beta_vec <- as.numeric(x_err)            # p×2 → long vector
      sd_vec   <- as.numeric(sd_post)
      
      beta_post_vec <- unlist(
        lapply(seq_along(beta_vec), function(i) 
          sign(beta_vec[i]) * max(abs(beta_vec[i]) - sd_vec[i]^2, 0))
      )
      beta_post <- matrix(beta_post_vec, nrow = nrow(x_err))
      
      #assign to unshrunk data objects
      x_err <- beta_post
      df <- df_post
      sd_x_err <- sd_post
    }
    
  } else {
    e <- matrix(rnorm(p * 2, sd = sqrt(err_var)), ncol = 2) %*% chol(error_R)
    x_err <- x + e
    sd_x_err <- matrix(err_var, nrow = p, ncol = 2)
    df <- cbind(rep(1E6, p), rep(1E6, p))
  }
  
  #plot these
  if(plot_scatter){
    
    par(mfrow = c(1,2), mar = c(5,5,5,1))
    plot(x, xlab = "true coefficient 1", 
         ylab = "true coefficient 2", 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(x)[1,2], 3))),
         col = adjustcolor(1, 0.3), pch = 19, asp = 1)
    abline(0,1, col = 2, lty = 2, lwd = 2)
    
    plot(x_err, xlab = "OLS estimated coefficient 1", 
         ylab = "OLS estimated coefficient 2", 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(x_err)[1,2], 3))),
         col = adjustcolor(1, 0.3), pch = 19, asp = 1)
    abline(0,1, col = 2, lty = 2, lwd = 2)
    
    
    #plot teaching example
    par(mfrow = c(1,3), mar = c(5,5,5,1))
    plot(x * 0.5 + 2.5, xlab = "True Length of Left Leg (ft.)", 
         ylab = "True Length of Right Leg (ft.)", 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(x)[1,2], 3))),
         col = adjustcolor(1, 0.3), pch = 19, asp = 1)
    abline(0,1, col = 2, lty = 2, lwd = 2)
    legend("topleft", lty = c(NA,2), col = c(adjustcolor(1, 0.5), 2), pch = c(19, NA),
           legend = c("indiv. obs", "1-to-1 line"), lwd = c(NA, 2))
    
    plot(x_err * 0.5 + 2.5, xlab = "Obs. Length of Left Leg (ft.)", 
         ylab = "Obs. Length of Right Leg (ft.)", 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(x_err)[1,2], 3))),
         col = adjustcolor(1, 0.3), pch = 19, asp = 1)
    abline(0,1, col = 2, lty = 2, lwd = 2)
    legend("topleft", lty = c(NA,2), col = c(adjustcolor(1, 0.5), 2), pch = c(19, NA),
           legend = c("indiv. obs", "1-to-1 line"), lwd = c(NA, 2))
    
    plot(x_err * 0.5 + 2.5, xlab = "Obs. Length of Left Leg (ft.)", 
         ylab = "Obs. Length of Right Leg (ft.)", 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(x_err)[1,2], 3))),
         col = adjustcolor(1, 0.3), pch = 19, asp = 1)
    for(i in 1:p){
      segments(x0 = x_err[i,1] * 0.5 + 2.5 - sd_x_err[i,1] / 2,
               x1 = x_err[i,1] * 0.5 + 2.5 + sd_x_err[i,1] / 2,
               y0 = x_err[i,2] * 0.5 + 2.5,
               y1 = x_err[i,2] * 0.5 + 2.5, lwd = 0.5, col = adjustcolor(1, 0.75))
      segments(x0 = x_err[i,1] * 0.5 + 2.5,
               x1 = x_err[i,1] * 0.5 + 2.5,
               y0 = x_err[i,2] * 0.5 + 2.5 - sd_x_err[i,2] / 2,
               y1 = x_err[i,2] * 0.5 + 2.5 + sd_x_err[i,2] / 2, 
               lwd = 0.5, col = adjustcolor(1, 0.75))
    }
    abline(0,1, col = 2, lty = 2, lwd = 2)
    legend("topleft", lty = c(NA, 2, 1), col = c(adjustcolor(1, 0.5), 2, adjustcolor(1, 0.75)), pch = c(19, NA, NA),
           legend = c("indiv. obs", "1-to-1 line", "±1SE"), lwd = c(NA, 2, 0.74))
    
  }
  
  #evaluate sample correlations
  cov(x_err)
  cor(x)[1,2] #should be approximately r
  cor(x_err)[1,2] #should be regressed from r towards 0
  
  #### preprocess model ####
  
  #get mixture approximation weights
  unique_df <- sort(unique(c(df)))           # all distinct dfs
  K <- 3
  w_tab  <- array(NA, dim = c(length(unique_df), K))
  s2_tab  <- array(NA, dim = c(length(unique_df), K))
  for (j in seq_along(unique_df)) {
    m <- make_mix(unique_df[j])
    w_tab[j, ] <- m$w
    s2_tab[j, ] <- m$s2_extra
  }
  
  # Map each row/column df to its table index
  df_index <- match(df, unique_df)
  
  
  #pack data into a list
  dat <- list(p=p, 
              x_err=x_err, 
              sd_x_err=sd_x_err * inflate_reported_error_by,
              df=df)
  
  #or for the horseshoe version
  if(grepl("horseshoe", model_path)){
    dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err, n = mean(n),
                p0 = c(rep(p, ifelse(grepl("unpooled", model_path), 2, 1))) * (1-prop_null))
  }
  
  #or for the logit rescale version
  if(grepl("logit-rescale", model_path)){
    dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err, 
                p0 = c(rep(p, ifelse(grepl("unpooled", model_path), 2, 1))) * (1-prop_null))
  }
  
  #### fit model ####
  mod <- cmdstan_model(model_path)
  fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 2E3,
                    data = dat, adapt_delta = 0.9, parallel_chains = 4,
                    refresh = 100, max_treedepth = 12, 
                    thin = 1, init = 0.1)
  
  #### mcmc diagnostics ####
  
  #check convergence
  if(check_diag){
    summ <- fit$summary()
    print(summ[order(summ$ess_bulk),])
    print(summ[order(summ$rhat, decreasing = T),])
  }
  
  #extract samples and inspect
  samps_r <- data.frame(as_draws_df(fit$draws("r")))
  
  #### plotting ####
  par(mfrow = c(1,1), oma = c(1,1,3,1))
  cols <- adjustcolor(c(2,4), 0.5)
  cols_vec <- rep(cols, each = p)
  shuffi <- sample(1:(p*2))
  
  #plot marginal posterior
  range_breaks <- range(c(samps_r$r, r, cor(x_err)[1,2]))
  breaks <- floor(range_breaks[1]*20):ceiling(range_breaks[2]*20)/20
  hist(samps_r$r, breaks = breaks, 
       main = "", freq = F, 
       xlab = "correlation", xlim = c(-1,1))
  mtext(paste0("true corr. = ", round(r, 3), ", sample size = ", p, 
               ", range errors = [", round(min(sd_x_err), 2), ", ", round(max(sd_x_err), 2), "]"), 
        side = 3, line = 4, cex = 1.25)
  
  #label lines
  abline(v = r, col = "red", lwd = 3)
  text(x = r, y = par("usr")[4], pos = 3, labels = "true\ncorr.", xpd = NA, col = "red")
  abline(v = cor(x_err)[1,2], col = "blue", lwd = 3)
  text(x = cor(x_err)[1,2], y = par("usr")[4], 
       pos = 3, labels = "sample\ncorr.", xpd = NA, col = "blue")
  segments(x0 = mean(samps_r$r), x1 = mean(samps_r$r), 
           y0 = par("usr")[3], y1 = par("usr")[3] -  diff(par("usr")[3:4])/5,
           col = "purple", lwd = 3, xpd = NA)
  text(x = mean(samps_r$r), par("usr")[3] -  diff(par("usr")[3:4])/5, 
       pos = 1, labels = "pmean\ncorr.", xpd = NA, col = "purple")
  
  plot_extra_plots <- F
  if(plot_extra_plots){
    
    samps_x <- data.frame(as_draws_df(fit$draws("x")))
    
    
    #look at estimates of x
    est_x_pmeans <- apply(samps_x, 2, mean)
    x_pmeans <- cbind(est_x_pmeans[1:p], est_x_pmeans[(1+p):(2*p)])
    est_x_psds <- apply(samps_x, 2, sd)
    x_psds <- cbind(est_x_psds[1:p], est_x_psds[(1+p):(2*p)])
    
    #plot means
    par(mfrow = c(1,2), mar = c(5,5,5,1))
    plot(x[shuffi], x_err[shuffi], xlab = "true coefficient", 
         ylab = latex2exp::TeX("\\textbf{OLS estimate} of coefficient"), 
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(cbind(c(x), c(x_err)))[1,2], 3))),
         col = cols_vec[shuffi], pch = 19)
    abline(0,1, col = 1, lty = 2, lwd = 2)
    
    plot(x[shuffi], x_pmeans[shuffi], xlab = "true coefficient", 
         ylab = latex2exp::TeX("\\textbf{posterior mean} of coefficient"),
         main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                      round(cor(cbind(c(x), c(x_pmeans)))[1,2], 3))),
         col = cols_vec[shuffi], pch = 19)
    abline(0,1, col = 1, lty = 2, lwd = 2)
    
    #plot sds
    par(par_base)
    
    #scatterplot
    plot(log10(sd_x_err[shuffi]), log10(x_psds[shuffi]), pch = 19, 
         col = cols_vec[shuffi],
         xlab = latex2exp::TeX("log$_{10}$(OLS Standard Error)"), 
         ylab = latex2exp::TeX("log$_{10}$Posterior SD)"))
    legend("topleft", legend=c("Var 1", "Var 2"), 
           col = cols, pch = 19)
    abline(0,1,lwd=2,lty=2)
    
    #get histogram data
    sd_diffs <- log10(sd_x_err) - log10(x_psds) #positive means posterior sd is smaller
    breaks <- seq(min(sd_diffs), max(sd_diffs), length.out = 20)
    histall <- hist(sd_diffs, breaks = breaks, plot=FALSE)
    hist1 <- hist(sd_diffs[,1], breaks = breaks, plot=FALSE)
    hist2 <- hist(sd_diffs[,2], breaks = breaks, plot=FALSE)
    
    #munge to pseudolog for frequency
    hist1$counts <- log10(hist1$counts + 1)
    hist2$counts <- log10(hist2$counts + 1)
    histall$counts <- log10(histall$counts + 1)
    
    # Plot the stacked histogram
    par(fig=c(0.5, 0.95, 0.25, 0.5), mar = c(0,0,0,0), new=TRUE)
    plot(hist1, col = cols[1], xlab="", ylab="", main = "", xpd = NA, cex.axis=0.7)
    mtext(1, text = "Difference", line = 2, cex = 0.7)
    mtext(2, text = latex2exp::TeX("log$_{10}$(Frequency)"), line = 2, cex = 0.7)
    plot(hist2, col=cols[2], add=TRUE)
    
  }
  
  #reset par
  par(par_base)
  
  #### try out other methods ####
  
  # helper 
  iv_weights <- 1 / rowMeans(sd_x_err^2)          # inverse variance weights
  
  # inverse-variance-weighted Pearson 
  cor_weighted <- wCorr::weightedCorr(x_err[,1], x_err[,2],
                               weights = iv_weights,
                               method  = "Pearson")
  nboot <- 4E2
  cor_weighted_boot <- replicate(n = nboot, expr = {
    idx <- sample(1:nrow(x_err), replace = T)
    cor_weighted <- wCorr::weightedCorr(x_err[idx,1], x_err[idx,2],
                                        weights = iv_weights[idx],
                                        method  = "Pearson")
    cor_weighted
  }, simplify = T)
  
  # 4. empirical-Bayes shrinkage with ashr 
  fit1  <- ash(beta = x_err[,1], se = sd_x_err[,1])
  fit2  <- ash(beta = x_err[,2], se = sd_x_err[,2])
  beta1 <- get_pm(fit1);
  beta2 <- get_pm(fit2)
  cor_ashr <- cor(beta1, beta2)
  nboot <- 1E1
  cor_ashr_boot <- replicate(n = nboot, expr = {
    idx <- sample(1:nrow(x_err), replace = T)
    fit1  <- ash(beta = x_err[idx,1], se = sd_x_err[idx,1])
    fit2  <- ash(beta = x_err[idx,2], se = sd_x_err[idx,2])
    beta1 <- get_pm(fit1);
    beta2 <- get_pm(fit2)
    cor_ashr <- cor(beta1, beta2)
    cor_ashr
  }, simplify = T)
  
  #try the hetero_corr function too?
  hc_out <- hetero_cor(x = x_err[,1], y = x_err[,2], 
                       var_x = sd_x_err[,1]^2, var_y = sd_x_err[,2]^2, ci = T)
  
  #record results
  w_corrs[run_i] <- cor_weighted
  w_corrs_bs[[run_i]] <- cor_weighted_boot
  ashr_corrs[run_i] <- cor_ashr
  ashr_corrs_bs[[run_i]] <- cor_ashr_boot
  hetero_corrs[run_i] <- hc_out$rho
  hetero_corrs_bs[[run_i]] <- hc_out$boot_values
  true_corrs[run_i] <- r
  sample_corrs[run_i] <- cor(x_err)[1,2]
  sample_corrs_noerror[run_i] <- cor(x)[1,2]
  posterior_corrs[[run_i]] <- samps_r$r
  if(check_diag){
    summs[[run_i]] <- summ
  }
  
}

#### examine aggregated results ####
pmean_corrs <- sapply(posterior_corrs, mean)
pmode_corrs <- sapply(posterior_corrs, function(ps){
  ps_atanh <- atanh(ps)
  dens <- density(ps_atanh)
  tanh(dens$x[which.max(dens$y)])
})
calib_q <- sapply(1:nruns, function(i) mean(posterior_corrs[[i]] > true_corrs[i]))
hc_calib_q <- sapply(1:nruns, function(i) mean(hetero_corrs_bs[[i]] > true_corrs[i]))
wc_calib_q <- sapply(1:nruns, function(i) mean(w_corrs_bs[[i]] > true_corrs[i]))
ac_calib_q <- sapply(1:nruns, function(i) mean(ashr_corrs_bs[[i]] > true_corrs[i]))
ess <- sapply(posterior_corrs, ess_basic)
rhats <- sapply(posterior_corrs, rhat)
ess_thresh <- 100
rhat_thresh <- 1.02
good_runs <- rep(T, nruns)
good_runs <- abs(true_corrs) < 0.65
good_runs <- ess > ess_thresh & rhats < rhat_thresh
if(check_diag){
  min_ess <- sapply(summs, function(si) min(c(si$ess_bulk, si$ess_tail)))
  max_rhat <- sapply(summs, function(si) max(si$rhat))
  good_runs <- min_ess > ess_thresh & max_rhat < rhat_thresh
}

#subset to only those with passing MCMC diags
tc <- true_corrs[good_runs]
sc <- sample_corrs[good_runs]
scne <- sample_corrs_noerror[good_runs]
pc <- pmean_corrs[good_runs]
pmc <- pmode_corrs[good_runs]
cq <- calib_q[good_runs]
hc <- hetero_corrs[good_runs]
hcq <- hc_calib_q[good_runs]
cq <- calib_q[good_runs]

wc <- w_corrs[good_runs]
ac <- ashr_corrs[good_runs]
wcq <- wc_calib_q[good_runs]
acq <- ac_calib_q[good_runs]

#### plot all results ####

par(mfcol = c(3,3), mar = c(4,5,4,3), oma = c(1,1,1,1))

plot(tc, sc, asp = 1, main = paste0("a) sample corrs (error), ccc = ", round(ccc(tc, sc), 2)),
     xlab = "true correlation", ylab = "sample correlation (error)")
abline(0,1,lwd=2,lty=2,col=2)

plot(tc, scne, asp = 1, main = paste0("b) sample corrs (no error), ccc = ", round(ccc(tc, scne), 2)),
     xlab = "true correlation", ylab = "sample correlation (no error)", xpd = NA)
abline(0,1,lwd=2,lty=2,col=2)

plot(tc, pc, asp = 1, main = paste0("c) ccc = ", round(ccc(tc, pc), 2)),
     xlab = "true correlation", ylab = "MY posterior mean correlation")
abline(0,1,lwd=2,lty=2,col=2)

hist(cq, breaks = 0:10/10, freq = F, xlab = "posterior quantile of true correlation", 
     main = "d) calibration of MY\ndisattenuated correlation")
abline(h=1,lwd=2,lty=2,col=2)

#get average correlation in there
interval_inds <- cut(cq, breaks = 0:10/10, include.lowest = T, labels = F)
cqmeans_in_interval <- sapply(split(cq, interval_inds), mean)
tcorrmeans_in_interval <- sapply(split(abs(tc), interval_inds), mean)
rax_ticks <- pretty(tcorrmeans_in_interval)
proj_to_axis <- function(rawx){
  range_rawx <- range(rawx)
  yax <- par("usr")[3:4]
  rawx01 <- (rawx - min(rawx)) / diff(range_rawx)
  trawx <- rawx01 * diff(yax) + min(yax)
  return(trawx)
}
proj_yval <- proj_to_axis(c(rax_ticks, tcorrmeans_in_interval))
lines(cqmeans_in_interval, proj_yval[-seq_along(rax_ticks)], type = "l", col = "blue")
axis(4, at = proj_yval[seq_along(rax_ticks)], labels = rax_ticks, col = "blue")
# mtext(side = 4, text = "E(|true corr|)", col = 1, cex = 0.85, line = 3)


# plot(tc, hc, asp = 1, main = paste0("e) algebraic correction, ccc = ", round(ccc(tc, hc), 2)),
#      xlab = "true correlation", ylab = "algebraically corrected correlation")
# abline(0,1,lwd=2,lty=2,col=2)

# hist(hcq, breaks = 0:10/10, freq = F, xlab = "posterior quantile of true correlation", 
#      main = "f) calibration of algebra-\nically corrected correlation")
# abline(h=1,lwd=2,lty=2,col=2)

plot(tc, wc, asp = 1, main = paste0("e) inverse variance\nweighted correction, ccc = ", round(ccc(tc, wc), 2)),
     xlab = "true correlation", ylab = "weighted correlation")
abline(0,1,lwd=2,lty=2,col=2)

hist(wcq, breaks = 0:10/10, freq = F, xlab = "posterior quantile of true correlation",
     main = "f) calibration of inverse-\nvariance weighted correlation")
abline(h=1,lwd=2,lty=2,col=2)

#add in bin-wise correlation (to check if "conservative" vs "anti-conservative")
interval_inds <- cut(wcq, breaks = 0:10/10, include.lowest = T, labels = F)
wcqmeans_in_interval <- sapply(split(wcq, interval_inds), mean)
tcorrmeans_in_interval <- sapply(split(abs(tc), interval_inds), mean)
rax_ticks <- pretty(tcorrmeans_in_interval)
proj_to_axis <- function(rawx){
  range_rawx <- range(rawx)
  yax <- par("usr")[3:4]
  rawx01 <- (rawx - min(rawx)) / diff(range_rawx)
  trawx <- rawx01 * diff(yax) + min(yax)
  return(trawx)
}
proj_yval <- proj_to_axis(c(rax_ticks, tcorrmeans_in_interval))
lines(wcqmeans_in_interval, proj_yval[-seq_along(rax_ticks)], type = "l", col = "blue")
axis(4, at = proj_yval[seq_along(rax_ticks)], labels = rax_ticks, col = "blue")
mtext(side = 4, text = "E(|true corr|)", col = 1, line = 2, cex = 0.6)

#is calibration concentrated in particular regions of correlation space?
# plot.new()
# plot.new()
# plot.new()
# plot(abs(tc), cq, xlab = "absolute true correlation |r|", 
#      ylab = "calibration of\ndisattenuated correlation", xpd = NA,
#      main = "g) does miscalibration systematically\nreflect particular true correlations?")
# breaks <- 0:10/10
# bmidp <- breaks[-1] - diff(breaks) / 2
# quantiles_of_cq <- sapply(seq_along(breaks)[-1], function(i){
#   bounds <- c(breaks[i-1], breaks[i])
#   quantile(cq[tc > bounds[1] & tc < bounds[2]], c(0.25, 0.5, 0.75))
# })
# lines(bmidp, quantiles_of_cq[2,], col = 2, lwd = 2)
# lines(bmidp, quantiles_of_cq[1,], col = 2, lwd = 2)
# lines(bmidp, quantiles_of_cq[3,], col = 2, lwd = 2)

# par(oma = c(0,1,0,0))Ç=

plot(abs(tc), abs(cq-0.5), xlab = "absolute true correlation |r|", 
     ylab = "calibration of\ndisattenuated correlation", xpd = NA,
     main = "g) does miscalibration systematically\nreflect particular true correlations?")
breaks <- 0:5/5
bmidp <- breaks[-1] - diff(breaks) / 2
quantiles_of_cq <- sapply(seq_along(breaks)[-1], function(i){
  bounds <- c(breaks[i-1], breaks[i])
  quantile(abs(cq-0.5)[tc > bounds[1] & tc < bounds[2]], c(0.25, 0.5, 0.75))
})
lines(bmidp, quantiles_of_cq[1,], col = 2, lwd = 2)
lines(bmidp, quantiles_of_cq[2,], col = 2, lwd = 5)
lines(bmidp, quantiles_of_cq[3,], col = 2, lwd = 2)
polygon(x = c(bmidp, rev(bmidp)), y = c(quantiles_of_cq[1,], rev(quantiles_of_cq[3,])), border = NA, col = adjustcolor(2, 0.25))


plot(abs(pc), abs(cq-0.5), xlab = "absolute posterior mean correlation", 
     ylab = "calibration of\ndisattenuated correlation", xpd = NA,
     main = "h) does miscalibration systematically\nreflect particular posterior geometries?")
breaks <- 0:5/5
bmidp <- breaks[-1] - diff(breaks) / 2
quantiles_of_cq <- sapply(seq_along(breaks)[-1], function(i){
  bounds <- c(breaks[i-1], breaks[i])
  quantile(abs(cq-0.5)[pc > bounds[1] & pc < bounds[2]], c(0.25, 0.5, 0.75))
})
lines(bmidp, quantiles_of_cq[1,], col = 2, lwd = 2)
lines(bmidp, quantiles_of_cq[2,], col = 2, lwd = 5)
lines(bmidp, quantiles_of_cq[3,], col = 2, lwd = 2)
polygon(x = c(bmidp, rev(bmidp)), 
        y = c(quantiles_of_cq[1,], rev(quantiles_of_cq[3,])), 
        border = NA, col = adjustcolor(2, 0.25))

#legend
plot.new()
legend("topleft", legend = c("perfect calibration", 
                             "average true correlation in window",
                             "quartiles of calibration"),
       bty = "n", lwd = 2, lty = c(2,1,NA), pch = c(NA, NA, 15), pt.cex =  c(NA, NA, 2),
       col = c(2,"blue",2))
