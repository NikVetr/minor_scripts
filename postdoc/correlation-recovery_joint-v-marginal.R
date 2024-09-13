library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#specify a few functions
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

# simulate data
D <- 15
R <- rlkj(D, 1)
# r <- 0.5
# R <- diag(D) + r - diag(D) * r
N <- 20
y <- mvrnorm(N, rep(0,D), R)
cor(y)

# construct data object
dat <- list(D = D,
            N = N,
            y = y)

# construct data object from fisher z transform

#first evaluate pairwise correlations
combos <- t(apply(expand.grid(1:D, 1:D), 1, sort))
combos <- combos[!duplicated(combos),]
combos <- combos[!apply(combos, 1, function(x) x[1] == x[2]),]
pairwise_corrs <- data.frame(do.call(rbind, lapply(1:choose(D, 2), function(i){
  r <- cor(y[,combos[i,1]], y[,combos[i,2]])
  n <- nrow(y)
  z <- 0.5 * log((1 + r) / (1 - r))
  sigma_z <- 1 / sqrt(n - 3)
  return(c(z = z, sigma_z = sigma_z))
})))

dat_fisher_z <- list(p = D,
                     np = choose(D,2),
                     ip = combos,
                     z = pairwise_corrs$z,
                     sigma_z = pairwise_corrs$sigma_z)



# STAN model
stan_loc_joint <- "~/scripts/minor_scripts/postdoc/multi-normal_joint.stan"
stan_program_joint <- paste0(readLines(stan_loc_joint), collapse = "\n")
mod_joint <- cmdstan_model(stan_loc_joint)

stan_loc_marg <- "~/scripts/minor_scripts/postdoc/multi-normal_marginal.stan"
stan_program_marg <- paste0(readLines(stan_loc_marg), collapse = "\n")
mod_marg <- cmdstan_model(stan_loc_marg)

stan_loc_marg_fb <- "~/scripts/minor_scripts/postdoc/multi-normal_marginal_flat-beta.stan"
stan_program_marg_fb <- paste0(readLines(stan_loc_marg_fb), collapse = "\n")
mod_marg_fb <- cmdstan_model(stan_loc_marg_fb)

stan_loc_fisher_z <- "~/scripts/minor_scripts/postdoc/correlation-matrix_from-bivariate-Fisher-z.stan"
stan_program_fisher_z <- paste0(readLines(stan_loc_fisher_z), collapse = "\n")
mod_fisher_z <- cmdstan_model(stan_loc_fisher_z)

#fit model
fit_joint <- mod_joint$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, 
                            parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, refresh = 100, init = 0.1)

fit_marg <- mod_marg$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, 
                            parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, refresh = 100, init = 0.1)

fit_marg_fb <- mod_marg_fb$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, 
                                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, refresh = 100, init = 0.1)

fit_fisher_z <- mod_fisher_z$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat_fisher_z, 
                                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, refresh = 100, init = 0.1)

#mcmc diagnostics
summ_joint <- fit_joint$summary("R")
summ_joint[order(summ$ess_bulk),]

summ_marg <- fit_marg$summary("R")
summ_marg[order(summ_marg$ess_bulk),]

summ_marg_fb <- fit_marg_fb$summary("R")
summ_marg_fb[order(summ_marg_fb$ess_bulk),]

summ_fisher_z <- fit_fisher_z$summary("R")
summ_fisher_z[order(summ_fisher_z$ess_bulk),]

#retrieve samples
samps_joint <- data.frame(as_draws_df(fit_joint$draws("R")))
samps_marg <- data.frame(as_draws_df(fit_marg$draws("R")))
samps_marg_fb <- data.frame(as_draws_df(fit_marg_fb$draws("R")))
samps_fisher_z <- data.frame(as_draws_df(fit_fisher_z$draws("R")))

samples <- list(
  joint = samps_joint,
  marginal = samps_marg,
  marginal_fb = samps_marg_fb,
  fisher_z = samps_fisher_z
)


# process samples
post <- lapply(setNames(names(samples), names(samples)), function(samps_name){
  
  #retrieve object
  print(samps_name)
  samps <- samples[[samps_name]]
  
  #munge to array 
  corr_mats <- lapply(1:nrow(samps), function(ri) 
    matrix(c(as.matrix(samps[ri, grep("R", colnames(samps))])), ncol = D))
  corr_mats_array <- do.call(abind::abind, c(corr_mats, list(along = 3)))
  
  #compute posterior summaries
  mean_corr_mat <- apply(corr_mats_array, c(1,2), function(x) mean(x))
  var_corr_mat <- apply(corr_mats_array, c(1,2), function(x) var(x))
  q_corr_mat <- apply(combn(D, 2), 2, function(i){
    mean(corr_mats_array[i[1], i[2],] < R[i[1], i[2]])
  })
  
  return(list(mean_R = mean_corr_mat[upper.tri(R)], 
              var_R = var_corr_mat[upper.tri(R)], 
              q_R = q_corr_mat))
  
})

#recover means
true_R <- R[upper.tri(R)]
means <- lapply(post, function(x) x$mean_R)
vars <- lapply(post, function(x) x$var_R)
qs <- lapply(post, function(x) x$q_R)


# panel function to show Pearson correlation
panel.cor <- function(x, y, digits = 4, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0("r = ", prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.5 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

# panel function to draw scatterplot with abline
panel.smooth.abline <- function(x, y, ...) {
  points(x, y, ...)
  abline(0, 1, col = 2, lty = 2, lwd = 2)
}

# Call pairs with custom panel functions
pairs(list(true = true_R, means), 
      lower.panel = panel.cor, 
      upper.panel = panel.smooth.abline,
      diag.panel = NULL)

#now look at posterior variances 
pairs(vars, 
      lower.panel = panel.cor, 
      upper.panel = panel.smooth.abline,
      diag.panel = NULL)


#now look at calibration
par(mfrow = c(length(qs), 1))
for(i in seq_along(qs)){
  hist(qs[[i]], main = names(qs)[i], xlab = "P(R < true R)", 
       breaks = 0:10/10, xlim = c(0,1), freq = F)
}
