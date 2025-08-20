# -------- Poisson–exposure simulation + model-fitting ------------------------

simulate_and_fit <- function(n1      = 100,          # size of indicator==0 group
                             n2      = 100,          # size of indicator==1 group
                             exp_mu  = c(2, 3),      # meanlog for log-normal exposures
                             exp_sd  = c(1, 1),  # sdlog  for log-normal exposures
                             B       = 0,            # true coefficient on indicator
                             beta0   = 0,            # baseline log-rate
                             nreps   = 100,          # Monte-Carlo repetitions
                             tau = 100,
                             maxit = 1E4,
                             useNB = T) {         # concentration parameter for gampois
  
  # pre-allocate results data frame
  out <- data.frame(rep = integer(), model = character(),
                    beta_hat = numeric(), p_value = numeric())
  row_index <- 1
  
  for (r in seq_len(nreps)) {
    if(r %% 10 == 0){cat(paste0(r, " "))}
    
    # ---------- simulate one data set ----------------------------------------
    g0_exp <- rlnorm(n1, meanlog = exp_mu[1], sdlog = exp_sd[1])
    g1_exp <- rlnorm(n2, meanlog = exp_mu[2], sdlog = exp_sd[2])
    
    exposure   <- c(g0_exp, g1_exp)
    indicator  <- c(rep(0, n1), rep(1, n2))
    
    lambda     <- exp(beta0 + B * indicator) * exposure
    overdispersed_lambda <- sapply(lambda, function(li){
      gp <- gamma_params(mean = li, concentration = tau)
      return(rgamma(n = 1, shape = gp$shape, rate = gp$rate))
    })
    counts     <- rpois(length(lambda), overdispersed_lambda)
    
    df <- data.frame(y = counts,
                     ind = indicator,
                     exposure = exposure,
                     log_exposure = log(exposure))
    
    # ---------- model 1: exposure as offset ----------------------------------
    if(useNB){
      fit1 <- MASS::glm.nb(y ~ 1 + ind + offset(log_exposure),
                  link = log, data = df, 
                  control = glm.control(maxit = maxit),
                  init.theta = 1.0,   
                  limit = maxit)      
    } else {
      fit1 <- glm(y ~ 1 + ind + offset(log_exposure),
                  family = poisson(link = "log"), data = df)
    }
    
    summ1 <- summary(fit1)$coefficients["ind", c("Estimate", "Pr(>|z|)")]
    
    # ---------- model 2: exposure (natural scale) as covariate --------------
    if(useNB){
      fit2 <- MASS::glm.nb(y ~ 1 + ind + exposure, link = log, data = df, 
                           control = glm.control(maxit = maxit),
                           init.theta = 1.0,   
                           limit = maxit)  
    } else {
      fit2 <- glm(y ~ 1 + ind + exposure,
                  family = poisson(link = "log"), data = df)  
    }
    
    summ2 <- summary(fit2)$coefficients["ind", c("Estimate", "Pr(>|z|)")]
    
    # ---------- model 3: log(exposure) as covariate --------------------------
    if(useNB){
      fit3 <- MASS::glm.nb(y ~ 1 + ind + log_exposure, link = log, data = df, 
                           control = glm.control(maxit = maxit),
                           init.theta = 1.0,   
                           limit = maxit)  
    } else {
      fit3 <- glm(y ~ 1 + ind + log_exposure,
                  family = poisson(link = "log"), data = df)  
    }
    
    summ3 <- summary(fit3)$coefficients["ind", c("Estimate", "Pr(>|z|)")]
    
    # ---------- store results ------------------------------------------------
    out <- rbind(out,
                 data.frame(rep = r, model = "offset",           t(summ1)),
                 data.frame(rep = r, model = "linear_exposure",  t(summ2)),
                 data.frame(rep = r, model = "log_exposure",     t(summ3)))
  }
  colnames(out) <- c("res", "model", "beta_hat", "p_value")
  out
}

# ------------------ example call ---------------------------------------------
lognormal_params <- function(mean, sd) {
  # Calculate the lognormal parameters
  sigma <- sqrt(log(1 + (sd / mean)^2))
  mu <- log(mean) - (sigma^2) / 2
  
  # Return as a named list
  list(mu = mu, sigma = sigma)
}

gamma_params <- function(mean, concentration, parameterization = "rate") {
  # Check for valid inputs
  if (mean <= 0) stop("Mean must be positive.")
  if (concentration <= 0) stop("Concentration must be positive.")
  
  # Calculate shape and rate (or scale)
  shape <- concentration
  scale <- mean / concentration
  
  if (parameterization == "rate") {
    rate <- 1 / scale
    return(list(shape = shape, rate = rate))
  } else if (parameterization == "scale") {
    return(list(shape = shape, scale = scale))
  } else {
    stop("Parameterization must be either 'rate' or 'scale'.")
  }
}

# Example usage
params_rate <- gamma_params(1100, 2, parameterization = "rate")
params_scale <- gamma_params(1100, 2, parameterization = "scale")

print(params_rate)
print(params_scale)


res <- simulate_and_fit(n1 = 508, 
                        n2 = 315, 
                        exp_mu = c(lognormal_params(1100, 4500)$mu, 
                                   lognormal_params(4200, 11200)$mu),
                        exp_sd = c(lognormal_params(1100, 4500)$sigma, 
                                   lognormal_params(4200, 11200)$sigma), 
                        B = 0, beta0 = 1,
                        nreps = 1000)

# quick checks ---------------------------------------------------------------
# Average estimate of B and proportion of p < 0.05 for each model
aggregate(beta_hat ~ model, res, mean)
aggregate(p_value < 0.05 ~ model, res, mean)

pvals <- split(res$p_value, res$model)
coefs <- split(res$beta_hat, res$model)

par(mfrow = c(3,1))
hist(pvals$offset, breaks = 0:100/100, xlab = "p-value", freq = F, 
     main = "correct approach (using log(exposure) as an offset)")
abline(h = 1, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "expected density under null")

hist(pvals$linear_exposure, breaks = 0:100/100, xlab = "p-value", freq = F, 
     main = "Goodman et al. 2014 approach (including exposure as a feature)")
abline(h = 1, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "expected density under null")

hist(pvals$log_exposure, breaks = 0:100/100, xlab = "p-value", freq = F, 
     main = "Hope-of-recovery approach (including log(exposure) as a a feature)")
abline(h = 1, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "expected density under null")

xlims <- range(unlist(coefs))
breaks <- seq(xlims[1], xlims[2], length.out = 100)
hist(coefs$offset, breaks = breaks, 
     xlab = "estimate of focal param (coef on how much group 2 has a higher rate than group 1)", 
     freq = F, 
     main = "correct approach (using log(exposure) as an offset)")
abline(v = 0, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "true parameter value")

hist(coefs$linear_exposure, breaks = breaks, 
     xlab = "estimate of focal param (coef on how much group 2 has a higher rate than group 1)", 
     freq = F, 
     main = "Goodman et al. 2014 approach (including exposure as a feature)")
abline(v = 0, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "true parameter value")

hist(coefs$log_exposure, breaks = breaks, 
     xlab = "estimate of focal param (coef on how much group 2 has a higher rate than group 1)", 
     freq = F, 
     main = "Hope-of-recovery approach (including log(exposure) as a a feature)")
abline(v = 0, col = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, col = 2, legend = "true parameter value")
