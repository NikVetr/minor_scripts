library(cmdstanr)
library(posterior)
#specify functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1 + exp(x))
mod <- cmdstan_model("~/scripts/minor_scripts/postdoc/causal-recovery.stan")
simulate_causal_inference <- function(cholR, p = 50, n = 300, n_pred = 500, p_liab, 
                                      prop_yvar = 0.4, prop_zvar = 0.4, sd_b_xz = 1, b_yz = 0,
                                      standardize_vars = F, same_x_for_yz = F, incl_intercept = T){
  
  #simulate population of ternary ordinal correlated random variables for x -> y relationship
  x_xy <- Reduce("+", lapply(1:2, function(i){
    liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
    liabs <- t(t(liabs) + p_liab)
    probs <- pnorm(liabs)
    apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
  }))
  
  #simulate second population for x -> z relationship (or reuse the first)
  if(same_x_for_yz){
    
    x_xz <- x_xy
    
  } else {
    
    x_xz <- Reduce("+", lapply(1:2, function(i){
      liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
      liabs <- t(t(liabs) + p_liab)
      probs <- pnorm(liabs)
      apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
    }))  
    
  }
  
  #simulate independent population for prediction
  x_pred <- Reduce("+", lapply(1:2, function(i){
    liabs <- matrix(rnorm(p * n_pred), nrow = n_pred) %*% cholR
    liabs <- t(t(liabs) + p_liab)
    probs <- pnorm(liabs)
    apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
  }))
  
  #check these graphically for desired properties
  # hist(cor(x_xy)[upper.tri(cor(x_xy))])
  # plot(p_freqs, apply(x_xy, 2, mean) / 2); abline(0,1)
  
  #specify mediator, outcome, and exposure relationships
  b_xy <- rnorm(p)
  b_xz <- rnorm(p) * sd_b_xz #when sd_b_xz = 0, no horizontal pleiotropy 
  
  #simulate data for
  # z <- x -> y -> z
  mu_y <- c(x_xy %*% t(t(b_xy)))
  e_y <- rnorm(n) * c(sqrt(var(mu_y) / prop_yvar * (1-prop_yvar)))
  y <- mu_y + e_y
  
  #simulate z through different set of y from different population of x
  mu_y_xz <- c(x_xz %*% t(t(b_xy)))
  e_y_xz <- rnorm(n) * c(sqrt(var(mu_y_xz) / prop_yvar * (1-prop_yvar)))
  y_xz <- mu_y_xz + e_y_xz
  
  mu_z <- y_xz * b_yz + c(x_xz %*% t(t(b_xz)))
  e_z <- rnorm(n) * c(sqrt(var(mu_z) / prop_zvar * (1-prop_zvar)))
  z <- mu_z + e_z
  
  #standardize variables?
  if(standardize_vars){
    x_xy <- as.matrix(scale(x_xy))
    x_xz <- as.matrix(scale(x_xz))
    x_pred <- as.matrix(scale(x_pred))
    y <- c(scale(y))
    z <- c(scale(z))
  }
  
  #name the coefs explicitly
  snps <- paste0("SNP", 1:p)
  colnames(x_xy)  <- snps
  colnames(x_xz)  <- snps
  colnames(x_pred)  <- snps
  df_xy <- data.frame(y = y,  x_xy)          # 50 SNP columns + y
  df_xz <- data.frame(z = z,  x_xz)          # 50 SNP columns + z
  
  if(incl_intercept){
    fit_xy <- lm(y ~ ., data = df_xy)      # “- 1” drops a duplicate intercept
    fit_xz <- lm(z ~ ., data = df_xz)
  } else {
    fit_xy <- lm(y ~ . - 1, data = df_xy)      # “- 1” drops a duplicate intercept
    fit_xz <- lm(z ~ . - 1, data = df_xz)  
  }
  
  # #fit multiple regression models 
  # #(equivalent to sweeping univariate regression models)
  # fit_xy <- lm(y ~ 1 + x_xy)
  # fit_xz <- lm(z ~ 1 + x_xz)
  
  #retrieve coefficients and their standard errors
  inds <- setdiff(1:(p+incl_intercept), incl_intercept)
  coefs_xy <- fit_xy$coefficients[inds]
  SE_xy <- summary(fit_xy)$coefficients[inds,2]
  coefs_xz <- fit_xz$coefficients[inds]
  SE_xz <- summary(fit_xz)$coefficients[inds,2]
  
  # summary(lm(coefs_xz ~ 1 + coefs_xy))$coefficients[2,] # recovers b_yz
  # summary(lm(coefs_xy ~ 1 + coefs_xz))$coefficients[2,] # "recovers" reverse causal relationship, y <- x -> z -> y
  coefreg_fit <- summary(lm(coefs_xz ~ 1 + coefs_xy))$coefficients
  #these also get attenuated by estimation error!
  #try Bayesian method
  dat <- list(
    p = p,
    b_xy_est = coefs_xy, 
    b_xy_se = SE_xy,
    b_xz_est = coefs_xz, 
    b_xz_se = SE_xz
  )
  fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                    adapt_delta = 0.9, parallel_chains = 4,
                    refresh = 10, max_treedepth = 10, 
                    thin = 1, init = 0.1, data = dat)
  # summ <- fit$summary()
  # print(summ[order(summ$rhat, decreasing = T),c("variable", "rhat", "ess_bulk", "ess_tail")])
  samps <- data.frame(as_draws_df(fit$draws()))
  # hist(samps$tau, freq = F)
  # curve(dt(x=x, df = 4), from = 0, to = 4, add = T)
  # hist(samps$beta_yz)
  # plot(dat$b_xz_est, dat$b_xy_est)
  bayesian_b_yz <- samps$beta_yz
  bayesian_est <- mean(bayesian_b_yz)
  bayesian_pval <- min(mean(bayesian_b_yz > 0), mean(bayesian_b_yz < 0)) * 2
  bayesian_tau <- mean(samps$tau)
  
  #try the "imputed values" approach? predict y and z from x
  #following procedure here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552594/figure/F2/
  
  # y_pred <- predict(fit_xy, newdata = data.frame(x_pred))
  # y_pred_zpop <- predict(fit_xy, newdata = data.frame(x_xz))
  # z_pred <- predict(fit_xz, newdata = data.frame(x_pred))
  # indep_pop_prediction_fit <- summary(lm(z_pred ~ 1 + y_pred))$coefficients
  # predixcan_fit <- summary(lm(z ~ 1 + y_pred_zpop))$coefficients
  # 
  # #better predixcan estimate?
  # y_pred_zpop <- predict(fit_xy, newdata = data.frame(x_xz))
  # rho         <- cor(y_xz, y_pred_zpop)
  # beta_raw    <- coef(lm(z ~ y_pred_zpop))[2]
  # beta_cal    <- beta_raw / rho   # ≈ true b_yz
  # predixcan_fit[2,1] <- beta_cal
  
  #predixcan requires a new population to test this out in
  y_pred        <- predict(fit_xy, newdata = data.frame(x_pred))
  y_pred_zpop   <- predict(fit_xy, newdata = df_xz)      # <- **names match**
  z_pred        <- predict(fit_xz, newdata = data.frame(x_pred))
  indep_pop_prediction_fit <- summary(lm(z_pred ~ 1 + y_pred))$coefficients
  rho_predixcan       <- cor(y_xz, y_pred_zpop)        # correlation in the *test* sample
  beta_raw_predixcan  <- coef(lm(z ~ y_pred_zpop))[2]
  beta_cal_predixcan  <- beta_raw_predixcan / rho_predixcan                # consistent for b_yz
  se_raw_predixcan <- summary(lm(z ~ y_pred_zpop))$coefficients[2,2]
  se_cal_predixcan <- se_raw_predixcan / rho_predixcan                     # delta-method
  p_cal_predixcan  <- 2 * pnorm(-abs(beta_cal_predixcan / se_cal_predixcan))
  
  #try the s-predixcan formula (eq 11 here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940825/#Equ11)
  spredixcan_z <- sum(coefs_xy * (apply(x_pred, 2, var) / var(y_pred)) * (coefs_xz / SE_xz))
  spredixcan_p <- 2 * pnorm(-abs(spredixcan_z))
  # SNP SDs in the GWAS sample
  sigma_x  <- sqrt(apply(x_xz, 2, var))

  # predicted expression in the GWAS sample
  y_pred_zpop <- predict(fit_xy, newdata = data.frame(x_xz))
  sigma_g <- sd(y_pred_zpop)

  # GWAS Z-scores for each SNP
  Z_snps <- coefs_xz / SE_xz

  spredixcan_z <- sum(coefs_xy *  (sigma_x / sigma_g) * Z_snps)
  spredixcan_p <- 2 * pnorm(-abs(spredixcan_z))
  spredixcan_b  <- sum(coefs_xy * sigma_x^2 * coefs_xz) / (sigma_g^2)
  Z_snps <- vapply(1:p, function(j){
    fit <- lm(z ~ x_xz[, j])
    coef(summary(fit))[2, 3]        # t statistic = Z
  }, numeric(1))
  beta_snps <- vapply(1:p, function(j){
    coef(summary(lm(z ~ x_xz[, j])))[2, 1]
  }, numeric(1))
  spredixcan_z <- sum(coefs_xy * (sigma_x / sigma_g) * Z_snps)
  spredixcan_p <- 2 * pnorm(-abs(spredixcan_z))
  spredixcan_b <- sum(coefs_xy * sigma_x^2 * beta_snps) / sigma_g^2   # Eq. 5
  
  #try the TWMR formula (eq (1) here https://www.nature.com/articles/s41467-019-10936-0#Sec2 
  #or here: https://www.nature.com/articles/s41467-021-25805-y#Sec12)
  # twmr_a <- c(solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xz))) *
  #   solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xy))))
  # 
  # # Partial derivatives (quick numerical approximation)
  # eps <- 1e-6
  # twmr_a_xz_eps <- sapply(1:p, function(i){
  #   eps_vec <- c(rep(0,i-1), eps, rep(0,p-i))
  #   c(solve(t(coefs_xz + eps_vec) %*% solve(cor(x_pred)) %*% t(t(coefs_xz + eps_vec))) *
  #   solve(t(coefs_xz + eps_vec) %*% solve(cor(x_pred)) %*% t(t(coefs_xy))))
  #   })
  # twmr_a_xy_eps  <- sapply(1:p, function(i){
  #   eps_vec <- c(rep(0,i-1), eps, rep(0,p-i))
  #   c(solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xz))) *
  #   solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xy + eps_vec))))
  # })
  # da_db <- (twmr_a_xz_eps - twmr_a) / eps
  # da_dg <- (twmr_a_xy_eps - twmr_a) / eps
  # 
  # # Variance of alpha using the Delta method
  # var_twmr_a <- t(da_db) %*% vcov(fit_xz)[inds,inds] %*% t(t(da_db)) + t(da_dg) %*% vcov(fit_xy)[inds,inds] %*% t(t(da_dg))
  # # SE of alpha
  # SE_twmr_a <- sqrt(var_twmr_a)
  # twmr_a_std <- twmr_a / SE_twmr_a
  # twmr_p <- 2 * pnorm(-abs(twmr_a_std)) #under normal approx, can also use pt()
  
  #only meant to apply for independent SNPs, so "C" can be swapped with identity
  
  # --- TWMR (one-gene version) ----------------------------------------
  LD  <- cor(x_xz)
  E   <- coefs_xy            # SNP → expression
  G   <- coefs_xz            # SNP → trait
  S   <- solve(LD)
  
  alpha_hat <- as.numeric( t(E) %*% S %*% G / (t(E) %*% S %*% E) )
  
  # Delta-method SE
  denom      <- as.numeric(t(E) %*% S %*% E)
  dG         <- as.vector(S %*% E / denom)
  dE         <- -alpha_hat * dG + (S %*% G - alpha_hat * S %*% E) / denom
  var_alpha  <- sum(dE^2 * SE_xy^2) + sum(dG^2 * SE_xz^2)
  se_alpha   <- sqrt(var_alpha)
  p_alpha    <- 2 * pnorm(-abs(alpha_hat / se_alpha))
  
  return(
    list(
      pvals = list(
        true_coef = NA,
        coef_reg = coefreg_fit[2,4],
        pred_reg = indep_pop_prediction_fit[2,4],
        # predix_reg = predixcan_fit[2,4],
        predix_reg = p_cal_predixcan,
        s_predix = spredixcan_p,
        twmr = p_alpha,
        bayes = bayesian_pval,
        bayes_tau = NA
    ), 
      coefs = list(
        true_coef = b_yz,
        coef_reg = coefreg_fit[2,1],
        pred_reg = indep_pop_prediction_fit[2,1],
        # predix_reg = predixcan_fit[2,1],
        predix_reg = beta_cal_predixcan,
        s_predix = spredixcan_b,
        twmr = alpha_hat,
        bayes = bayesian_est,
        bayes_tau = bayesian_tau)
    )
  )
  
}

#### do simulations + analysis ####

#specify instrument / exposure (x) properties
p <- 50
n <- 300
p_freqs <- rbeta(p, 2, 2)
p_liab <- qnorm(p_freqs)

#simulate correlated counts in {0,1,2}
rop <- 0.5
rs <- c(1, rop^(1:p / 3))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)
cholR <- chol(R)

nrep <- 2E2
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}
true_b_yz <- 4
out <- parallel::mclapply(1:nrep, function(i){
  if(i %% 50 == 0){mcprint(paste0(i, " "))}
  do.call(rbind, simulate_causal_inference(cholR = chol(diag(p)), 
                                           p_liab = p_liab, 
                                           b_yz = true_b_yz, prop_zvar = 0.75, prop_yvar = 0.75, 
                                           standardize_vars = F, same_x_for_yz = F))}, mc.cores = 12)
pvals <- data.frame(do.call(rbind, lapply(1:nrep, function(i) unlist(out[[i]]["pvals",]))))
coefs <- data.frame(do.call(rbind, lapply(1:nrep, function(i) unlist(out[[i]]["coefs",]))))

par(mfrow = c(5,1))
breaks <- 0:20/20
hist(pvals$coef_reg, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "my method (freq)")
# hist(pvals$pred_reg, breaks = breaks, 
#      xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "PrediXcan (2-pop)")
hist(pvals$bayes, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "my method (bayes)")
hist(pvals$predix_reg, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "PrediXcan (1-pop)")
hist(pvals$s_predix, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "S-PrediXcan (3-pop)")
hist(pvals$twmr, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "TWMR")

par(mfrow = c(5,1))
breaks_range <- range(do.call(cbind, coefs)[,-1])
breaks <- seq(breaks_range[1], breaks_range[2], length.out = 20)

hist(coefs$coef_reg, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"), 
     main = paste0("my method (freq),\n mse = ", round(mean((coefs$coef_reg - true_b_yz)^2), 4)))
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

# hist(coefs$pred_reg, breaks = breaks, 
#      xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"), main = paste0("PrediXcan (2-pop)")
# abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
# text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

hist(coefs$bayes, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"),
     main = paste0("my method (bayes),\n mse = ", round(mean((coefs$bayes - true_b_yz)^2), 4)))
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

hist(coefs$predix_reg, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"), 
     main = paste0("PrediXcan (1-pop),\n mse = ", round(mean((coefs$predix_reg - true_b_yz)^2), 4)))
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

hist(coefs$s_predix, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"), 
     main = paste0("S-PrediXcan (2-pop),\n mse = ", round(mean((coefs$s_predix - true_b_yz)^2), 4)))
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

hist(coefs$twmr, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ =/= 0)"), 
     main = paste0("TWMR,\n mse = ", round(mean((coefs$twmr - true_b_yz)^2), 4)))
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

