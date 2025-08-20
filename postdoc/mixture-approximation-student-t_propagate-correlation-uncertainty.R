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

#simulate data

#### simulate input ####
use_ols_fit <- T
check_diag <- T

#for ols fit
n <- c(10, 5) #define sample size across dimensions
e_sd_sd <- c(4, 8) #define differential power across two dimensions? or just get from sample size

#for static error fit
err_var <- 1
# err_var <- 16 * matrix(exp(rnorm(p*2)), ncol = 2, nrow = p)
error_r <- 0
error_R <- diag(2) * (1-error_r) + error_r
x_var <- 1
plot_scatter <- F

#simulate coefficients
r <- -0.7 #true correlation between coefficients
r <- runif(1,-1,1)
R <- diag(2) * (1-r) + r #corresponding correlation matrix
p <- 1E2 #total number of samples
x <- matrix(rnorm(p*2) * sqrt(x_var), ncol = 2) %*% chol(R) #sample true coefficients


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
  
  fits_1 <- do.call(rbind, lapply(1:p, function(i){
    sim_and_fit_lm(b = x[i,1], err_sd = e_sd[i,1], nobs = n[1])
  }))
  fits_2 <- do.call(rbind, lapply(1:p, function(i){
    sim_and_fit_lm(b = x[i,2], err_sd = e_sd[i,2], nobs = n[2])
  }))
  
  x_err <- cbind(fits_1[,"Estimate"], fits_2[,"Estimate"])
  sd_x_err <- cbind(fits_1[,"Std. Error"], fits_2[,"Std. Error"])
  df <- cbind(rep(n[1] - 2, p), rep(n[2] - 2, p))
  
} else {
  e <- matrix(rnorm(p * 2, sd = sqrt(err_var)), ncol = 2) %*% chol(error_R)
  x_err <- x + e
  sd_x_err <- matrix(err_var, nrow = p, ncol = 2)
  df <- cbind(rep(1E6, p), rep(1E6, p))
}

#### preprocess model ####

make_lag_mix <- function(df, K = 3) {
  a  <- df / 2                        # shape  = ν/2
  lg <- gauss.quad(K, kind = "laguerre", alpha = a - 1)
  
  ## nodes and raw weights for ∫ e^{-z} z^{a-1} f(z) dz
  z      <- lg$nodes                  # length-K
  w_star <- lg$weights                # positive
  
  ## convert to Gamma(shape=a, rate=b=a)  --------------------------
  b         <- a                      # rate
  lambda    <- z / b                  # λ_k  (scales for the Normals)
  w         <- w_star / gamma(a)      # normalise so Σ w_k = 1
  
  list(lambda = lambda, w = w)
}

unique_df <- sort(unique(c(df)))   # all different dfs in your matrix
mix       <- lapply(unique_df, make_lag_mix, K = 3)

idx1 <- match(df[,1], unique_df)
idx2 <- match(df[,2], unique_df)

fill <- function(index, slot) {
  out <- matrix(NA_real_, nrow(df), 3)
  for (j in seq_along(unique_df)){
    nmatch <- sum(index == j)
    if(nmatch > 0){
      out[index == j, ] <- matrix(mix[[j]][[slot]], nrow = nmatch, ncol = 3, byrow = T)
    }  
  }
  out
}

w1   <- fill(idx1, "w")
w2   <- fill(idx2, "w")
s2_1 <- 1 / fill(idx1, "lambda")
s2_2 <- 1 / fill(idx2, "lambda")

#pack data into a list
dat <- list(p=p, 
            x_err=x_err, 
            sd_x_err=sd_x_err,
            df=df,
            K=3,
            w1 = w1, w2 = w2, s2_1 = s2_1, s2_2 = s2_2)


#### fit model ####
mod <- cmdstan_model(model_path)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 2E3,
                  data = dat, adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 100, max_treedepth = 12, 
                  thin = 1, init = 0.1)

#### mcmc diagnostics ####

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
