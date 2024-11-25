#specify functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1 + exp(x))
simulate_causal_inference <- function(cholR, p = 50, n = 300, n_pred = 500, p_liab, 
                                      prop_yvar = 0.4, prop_zvar = 0.4, sd_b_xz = 1, b_yz = 0,
                                      standardize_vars = F, same_x_for_yz = F){
  
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
  e_y <- rnorm(p) * c(sqrt(var(mu_y) / prop_yvar * (1-prop_yvar)))
  y <- mu_y + e_y
  
  #simulate z through different set of y from different population of x
  mu_y_xz <- c(x_xz %*% t(t(b_xy)))
  e_y_xz <- rnorm(p) * c(sqrt(var(mu_y_xz) / prop_yvar * (1-prop_yvar)))
  y_xz <- mu_y_xz + e_y_xz
  
  mu_z <- y_xz * b_yz + c(x_xz %*% t(t(b_xz)))
  e_z <- rnorm(p) * c(sqrt(var(mu_z) / prop_zvar * (1-prop_zvar)))
  z <- mu_z + e_z
  
  #standardize variables?
  if(standardize_vars){
    x_xy <- as.matrix(scale(x_xy))
    x_xz <- as.matrix(scale(x_xz))
    x_pred <- as.matrix(scale(x_pred))
    y <- c(scale(y))
    z <- c(scale(z))
  }
  
  #fit multiple regression models 
  #(equivalent to sweeping univariate regression models)
  fit_xy <- lm(y ~ 1 + x_xy)
  fit_xz <- lm(z ~ 1 + x_xz)
  
  #retrieve coefficients and their standard errors
  coefs_xy <- fit_xy$coefficients[-1]
  SE_xy <- summary(fit_xy)$coefficients[-1,2]
  coefs_xz <- fit_xz$coefficients[-1]
  SE_xz <- summary(fit_xz)$coefficients[-1,2]
  
  # summary(lm(coefs_xz ~ 1 + coefs_xy))$coefficients[2,] # recovers b_yz
  # summary(lm(coefs_xy ~ 1 + coefs_xz))$coefficients[2,] # "recovers" reverse causal relationship, y <- x -> z -> y
  coefreg_fit <- summary(lm(coefs_xz ~ 1 + coefs_xy))$coefficients
  #these also get attenuated by estimation error!
  
  #try the "imputed values" approach? predict y and z from x
  #following procedure here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552594/figure/F2/
  y_pred <- predict(fit_xy, newdata = data.frame(x_pred))
  y_pred_zpop <- predict(fit_xy, newdata = data.frame(x_xz))
  z_pred <- predict(fit_xz, newdata = data.frame(x_pred))
  indep_pop_prediction_fit <- summary(lm(z_pred ~ 1 + y_pred))$coefficients
  predixcan_fit <- summary(lm(z ~ 1 + y_pred_zpop))$coefficients
  
  #try the s-predixcan formula (eq 11 here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940825/#Equ11)
  spredixcan_z <- sum(coefs_xy * (apply(x_pred, 2, var) / var(y_pred)) * (coefs_xz / SE_xz))
  spredixcan_p <- 2 * pnorm(-abs(spredixcan_z))
  
  #try the TWMR formula (eq (1) here https://www.nature.com/articles/s41467-019-10936-0#Sec2 
  #or here: https://www.nature.com/articles/s41467-021-25805-y#Sec12)
  twmr_a <- c(solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xz))) *
    solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xy))))
  
  # Partial derivatives (quick numerical approximation)
  eps <- 1e-6
  twmr_a_xz_eps <- sapply(1:p, function(i){
    eps_vec <- c(rep(0,i-1), eps, rep(0,p-i))
    c(solve(t(coefs_xz + eps_vec) %*% solve(cor(x_pred)) %*% t(t(coefs_xz + eps_vec))) *
    solve(t(coefs_xz + eps_vec) %*% solve(cor(x_pred)) %*% t(t(coefs_xy))))
    })
  twmr_a_xy_eps  <- sapply(1:p, function(i){
    eps_vec <- c(rep(0,i-1), eps, rep(0,p-i))
    c(solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xz))) *
    solve(t(coefs_xz) %*% solve(cor(x_pred)) %*% t(t(coefs_xy + eps_vec))))
  })
  da_db <- (twmr_a_xz_eps - twmr_a) / eps
  da_dg <- (twmr_a_xy_eps - twmr_a) / eps
  
  # Variance of alpha using the Delta method
  var_twmr_a <- t(da_db) %*% vcov(fit_xz)[-1,-1] %*% t(t(da_db)) + t(da_dg) %*% vcov(fit_xy)[-1,-1] %*% t(t(da_dg))
  # SE of alpha
  SE_twmr_a <- sqrt(var_twmr_a)
  twmr_a_std <- twmr_a / SE_twmr_a
  twmr_p <- 2 * pnorm(-abs(twmr_a_std)) #under normal approx, can also use pt()
  
  #only meant to apply for independent SNPs, so "C" can be swapped with identity
  
  return(
    list(
      pvals = list(
        coef_reg = coefreg_fit[2,4],
        pred_reg = indep_pop_prediction_fit[2,4],
        predix_reg = predixcan_fit[2,4],
        s_predix = spredixcan_p,
        twmr = twmr_p
    ), 
      coefs = list(
        true_coef = b_yz,
        coef_reg = coefreg_fit[2,1],
        pred_reg = indep_pop_prediction_fit[2,1],
        predix_reg = predixcan_fit[2,1],
        s_predix = spredixcan_z,
        twmr = twmr_a)
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

nrep <- 1E3
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}
true_b_yz <- 0
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
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "my method")
hist(pvals$pred_reg, breaks = breaks, 
     xlab = latex2exp::TeX("p-values under null (B$_{y,z}$ = 0)"), main = "PrediXcan (2-pop)")
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
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ = 5)"), main = "my method")
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)
hist(coefs$pred_reg, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ = 5)"), main = "PrediXcan (2-pop)")
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)
hist(coefs$predix_reg, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ = 5)"), main = "PrediXcan (1-pop)")
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)
hist(coefs$s_predix, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ = 5)"), main = "S-PrediXcan (3-pop)")
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)
hist(coefs$twmr, breaks = breaks, 
     xlab = latex2exp::TeX("estimated coefficients under alternate (B$_{y,z}$ = 5)"), main = "TWMR")
abline(v = true_b_yz, lwd = 2, lty = 2, col = 2)
text(x = true_b_yz, y = par("usr")[4], labels = latex2exp::TeX("True B$_{y,z}$"), col = 2, xpd = NA, pos = 3, cex = 1.5)

