library(cmdstanr)
library(posterior)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#does iterative EM-flavored EB converge on Full Bayes MCMC?
#or maybe I can call this "maximum marginal likelihood"

#### specify high level parameters ####

#properties of iteration
accommodate_scale <- T
use_pathfinder <- F
use_errprop_mcmc <- T
propagate_scale_error <- F
n_iter <- 10
eps <- 1E-3

#simulate data
n_obs <- 30 #per group
n_groups <- 100
sd_a <- 0.5
sd_b <- 1
sd_e <- 2
a <- rnorm(n_groups, sd = sd_a)
b <- rnorm(n_groups, sd = sd_b)
# a <- rep(0, n_rep)
# b <- rep(0, n_rep)
xmat <- replicate(n_groups, rnorm(n_obs))
ymat <- sapply(1:n_groups, function(i) a[i] + b[i] * xmat[,i] + rnorm(n_obs, sd = sd_e))

#melt to vector
x <- c(xmat)
y <- c(ymat)
g <- rep(1:n_groups, each = n_obs)

#specify Stan model
stan_program <- "
data {
  int<lower=1> n_obs;
  int<lower=1> n_groups;
  array[n_obs*n_groups] int<lower=1, upper=n_groups> g;
  vector[n_obs*n_groups] x;
  vector[n_obs*n_groups] y;
}
parameters {
  vector[n_groups] a;
  vector[n_groups] b;
  real<lower=0> sd_a;
  real<lower=0> sd_b;
  real<lower=0> sd_e;
}
model {
  sd_a ~ std_normal();
  sd_b ~ std_normal();
  sd_e ~ std_normal();
  a ~ normal(0, sd_a);
  b ~ normal(0, sd_b);
  y ~ normal(a[g] + x .* b[g], sd_e);
}
"

stan_program_vb <- "
data {
  int<lower=1> n_obs;
  int<lower=1> n_groups;
  array[n_obs*n_groups] int<lower=1, upper=n_groups> g;
  vector[n_obs*n_groups] x;
  vector[n_obs*n_groups] y;
  real<lower=0> sd_a;
  real<lower=0> sd_b;
}
parameters {
  vector[n_groups] a;
  vector[n_groups] b;
  real<lower=0> sd_e;
}
model {
  a ~ normal(0, sd_a);
  b ~ normal(0, sd_b);
  y ~ normal(a[g] + x .* b[g], sd_e);
}
"

stan_program_error_scale <- "
data {
  int<lower=1> n_obs;
  int<lower=1> n_groups;
  array[n_obs*n_groups] int<lower=1, upper=n_groups> g;
  vector[n_obs*n_groups] x;
  vector[n_obs*n_groups] y;
  real log_sd_a;
  real log_sd_b;
  real<lower=0> sd_log_sd_a;
  real<lower=0> sd_log_sd_b;
}
parameters {
  vector[n_groups] true_log_sd_a;
  vector[n_groups] true_log_sd_b;
  vector[n_groups] a;
  vector[n_groups] b;
  real<lower=0> sd_e;
}
transformed parameters {
  vector[n_groups] sd_a = exp(true_log_sd_a);
  vector[n_groups] sd_b = exp(true_log_sd_b);
}
model {
  true_log_sd_a ~ normal(log_sd_a, sd_log_sd_a);
  true_log_sd_b ~ normal(log_sd_b, sd_log_sd_b);
  a ~ normal(0, sd_a);
  b ~ normal(0, sd_b);
  y ~ normal(a[g] + x .* b[g], sd_e);
}
"

mod <- cmdstan_model(write_stan_file(stan_program))
mod_vb <- cmdstan_model(write_stan_file(stan_program_vb))
mod_es <- cmdstan_model(write_stan_file(stan_program_error_scale))

#first fit w/ full bayes
d <- list(n_obs = n_obs, n_groups = n_groups, g = g, x = x, y = y)
fit <- mod$sample(d, chains = 4, parallel_chains = 4)
summ <- fit$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]
samps <- data.frame(as_draws_df(fit$draws()))
apply(samps[,grepl("sd_", colnames(samps))], 2, mean)

#now fit w/ VB once using wide priors
if(!use_errprop_mcmc){
  
  d_vb <- c(d, sd_a = 1E3, sd_b = 1E3)
  if(use_pathfinder){
    fit_approx <- mod_vb$pathfinder(data = d_vb, draws = 5E3)
    fit_approx_multilevel <- mod$pathfinder(data = d, draws = 5E3)
  } else {
    fit_approx <- mod_vb$variational(data = d_vb, algorithm="meanfield", draws = 5E3)
    fit_approx_multilevel <- mod$variational(data = d, algorithm="meanfield", draws = 5E3)
  }
  
  # fit_frvb <- mod_vb$variational(data = d_vb, algorithm="fullrank", draws = 5E3, iter = 5E4)
  samps_approx <- as.data.table(data.frame(as_draws_df(fit_approx$draws())))
  # samps_frvb <- as.data.table(data.frame(as_draws_df(fit_frvb$draws())))
  
  #compute the average sample variance across parameters
  scale_param_targets <- c("a", "b")
  scale_param_targets <- setNames(scale_param_targets, paste0("sd_", scale_param_targets))
  source("~/repos/Stan2R/R/functions.R")
  
  #calculate mean variance across parameters
  param_mean_var_approx <- sapply(scale_param_targets, function(scale_param){
    param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_approx)))
    param_vars <- apply(param_samps, 1, var)
    mean(param_vars)
  })
  
  # param_mean_var_frvb <- sapply(scale_param_targets, function(scale_param){
  #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
  #   param_vars <- apply(param_samps, 1, var)
  #   mean(param_vars)
  # })
  
  #and across iterations
  param_var_mean_approx <- sapply(scale_param_targets, function(scale_param){
    param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_approx)))
    param_vars <- apply(param_samps, 2, var)
    mean(param_vars)
  })
  
  # param_var_mean_frvb <- sapply(scale_param_targets, function(scale_param){
  #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
  #   param_vars <- apply(param_samps, 2, var)
  #   mean(param_vars)
  # })
  
  #see file '~/scripts/minor_scripts/postdoc/scale_mixture_variances.R' for derivation + demonstration
  if(accommodate_scale){
    param_sds_approx <- sqrt(param_mean_var_approx - param_var_mean_approx)  
  } else {
    param_sds_approx <- sqrt(param_mean_var_approx)
  }
  # param_sds_frvb <- sqrt(param_mean_var_frvb - param_var_mean_frvb)
  
  if(use_errprop_mcmc){
    param_sds_approx <- c(sd_a = mean(a_errp_samps$x_sigma),
                          sd_b = mean(b_errp_samps$x_sigma))
  }
  
  apply(samps[,names(scale_param_targets)], 2, mean) #"true" values
  param_sds_approx
  # param_sds_frvb
}

#check error propagation model?
d_flatmcmc <- c(d, sd_a = 1E3, sd_b = 1E3)
fit_flatmcmc <- mod_vb$sample(data = d_flatmcmc, chains = 4, parallel_chains = 4)
a_samps <- data.frame(as_draws_df(fit_flatmcmc$draws("a")))[,-((1:3)+n_groups)]
b_samps <- data.frame(as_draws_df(fit_flatmcmc$draws("b")))[,-((1:3)+n_groups)]
d_err_a <- list(n = n_groups, 
                x_est = apply(a_samps, 2, mean), 
                x_err_sd = apply(a_samps, 2, sd)
)
d_err_b <- list(n = n_groups, 
                x_est = apply(b_samps, 2, mean), 
                x_err_sd = apply(b_samps, 2, sd)
)

stan_program_errprop <- "
data {
  int<lower=1> n;
  vector[n] x_est;
  vector<lower=0>[n] x_err_sd;
}
parameters {
  real x_mu;
  real<lower=0> x_sigma;
  vector[n] x_raw;
}
transformed parameters {
  vector[n] x = x_raw * x_sigma + x_mu;
}
model {
  x_mu ~ std_normal();
  x_sigma ~ std_normal();
  x_raw ~ std_normal();
  x_est ~ normal(x, x_err_sd);
}
"

stan_program_errprop_fixed_mean <- "
data {
  int<lower=1> n;
  vector[n] x_est;
  vector<lower=0>[n] x_err_sd;
}
parameters {
  real<lower=0> x_sigma;
  vector[n] x_raw;
}
transformed parameters {
  vector[n] x = x_raw * x_sigma;
}
model {
  x_sigma ~ std_normal();
  x_raw ~ std_normal();
  x_est ~ normal(x, x_err_sd);
}
"

#fit
mod_errprop <- cmdstan_model(write_stan_file(stan_program_errprop_fixed_mean))
fit_errprop_a <- mod_errprop$sample(d_err_a, chains = 4, parallel_chains = 4)
fit_errprop_b <- mod_errprop$sample(d_err_b, chains = 4, parallel_chains = 4)

#extract samples and inspect results
a_errp_samps <- data.frame(as_draws_df(fit_errprop_a$draws()))
b_errp_samps <- data.frame(as_draws_df(fit_errprop_b$draws()))
mean(a_errp_samps$x_sigma)
mean(samps$sd_a)
mean(b_errp_samps$x_sigma)
mean(samps$sd_b)

#plot comparison
a_range <- range(c(a_errp_samps$x_sigma, samps$sd_a))
a_breaks <- seq(a_range[1], a_range[2], length.out = 20)
hist(a_errp_samps$x_sigma, freq = F, col = adjustcolor("orange", 0.5), breaks = a_breaks,
     ylim = c(0,9), xlab = latex2exp::TeX("$\\sigma_\\alpha$ (intercept scale)"),
     main = latex2exp::TeX("Posterior Samples for $\\sigma_\\alpha$ (intercept scale)"))
hist(samps$sd_a, add = T, freq = F, col = adjustcolor("blue",0.5), breaks = a_breaks)
legend("topright", pch = 15, col = adjustcolor(c("blue", "orange"), 0.5), bty = "n",
      legend = c("target distribution (full Bayes)", "meta-analytic approximation"))

b_range <- range(c(b_errp_samps$x_sigma, samps$sd_b))
b_breaks <- seq(b_range[1], b_range[2], length.out = 20)
hist(b_errp_samps$x_sigma, freq = F, col = adjustcolor("orange", 0.5), breaks = b_breaks,
     ylim = c(0,9), xlab = latex2exp::TeX("$\\sigma_\\beta$ (slope scale)"),
     main = latex2exp::TeX("Posterior Samples for $\\sigma_\\beta$ (slope scale)"))
hist(samps$sd_b, add = T, freq = F, col = adjustcolor("blue",0.5), breaks = b_breaks)
legend("topright", pch = 15, col = adjustcolor(c("blue", "orange"), 0.5), bty = "n",
       legend = c("target distribution (full Bayes)", "meta-analytic approximation"))

#compare to fitting multilevel model with VB
# samps_approx_multilevel <- data.frame(as_draws_df(fit_approx_multilevel$draws()))
# param_sds_approx_multilevel <- apply(samps_approx_multilevel[,names(scale_param_targets)], 2, mean)

#now try fitting the model with the fixed scale parameter, 
#either to the posterior mean or to the MoM normal, and see which does better
param_sds_pmeans <- c(sd_a = mean(a_errp_samps$x_sigma),
                      sd_b = mean(b_errp_samps$x_sigma))
#same as MoM estimator of the hyperprior scale mixture of normals, 
#bc of the law of total variance
d_fixed_hyperprior <- c(d, param_sds_pmeans)
fit_fixed_hyperprior <- mod_vb$sample(data = d_fixed_hyperprior, chains = 4, 
                                      parallel_chains = 4)
samps_fixed_hyperprior <- data.frame(as_draws_df(fit_fixed_hyperprior$draws()))

fixed_a_samps <- samps_fixed_hyperprior[,paste0("a.", 1:n_groups, ".")]
fixed_b_samps <- samps_fixed_hyperprior[,paste0("b.", 1:n_groups, ".")]
full_a_samps <- samps[,paste0("a.", 1:n_groups, ".")]
full_b_samps <- samps[,paste0("b.", 1:n_groups, ".")]

#compare means
plot(apply(full_a_samps, 2, mean), apply(fixed_a_samps, 2, mean)); abline(0,1,col=2)
plot(apply(full_b_samps, 2, mean), apply(fixed_b_samps, 2, mean)); abline(0,1,col=2)

#compare whole distributions
probs <- 0:100/100
plot(apply(full_a_samps, 2, quantile, probs = probs), 
     apply(fixed_a_samps, 2, quantile, probs = probs), type = "l"); abline(0,1,col=2,lwd=2)
plot(apply(full_b_samps, 2, quantile, probs = probs), 
     apply(fixed_b_samps, 2, quantile, probs = probs), type = "l"); abline(0,1,col=2, lwd=2)


#### try iteration? ####

#now iterate over this procedure
i <- 1
mean_diff <- Inf
approx_progress <- t(data.frame(c(param_sds_approx, mlp = mean(mean(fit_approx$lp())))))
# frvb_progress <- t(data.frame(param_sds_frvb))

while(i <= n_iter & mean_diff > eps){
  
  i <- i + 1
  print(i)
  
  #compile data
  d_approx <- c(d, param_sds_approx["sd_a"], param_sds_approx["sd_b"])
  # d_frvb <- c(d, param_sds_frvb["sd_a"], param_sds_frvb["sd_b"])
  
  #fit the model
  sink(tempfile())
  
  if(use_errprop_mcmc){
  
    if(propagate_scale_error){
      
      #compile data
      d_approx <- c(d,
                    log_sd_a = mean(log(a_errp_samps$x_sigma)),
                    log_sd_b = mean(log(b_errp_samps$x_sigma)),
                    sd_log_sd_a = sd(log(a_errp_samps$x_sigma)),
                    sd_log_sd_b = sd(log(b_errp_samps$x_sigma))
      )
      
      #fit model with previously determined scale
      fit_approx <- mod_es$sample(data = d_approx, chains = 4, parallel_chains = 4)
      
    } else {
      
      #compile data
      d_approx <- c(d, param_sds_approx["sd_a"], param_sds_approx["sd_b"])
      
      #fit model with previously determined scale
      fit_approx <- mod_vb$sample(data = d_approx, chains = 4, parallel_chains = 4)  
    }
    
    
    #extract samples from full model fit
    a_samps <- data.frame(as_draws_df(fit_approx$draws("a")))[,-((1:3)+n_groups)]
    b_samps <- data.frame(as_draws_df(fit_approx$draws("b")))[,-((1:3)+n_groups)]
    d_err_a <- list(n = n_groups, 
                    x_est = apply(a_samps, 2, mean), 
                    x_err_sd = apply(a_samps, 2, sd)
    )
    d_err_b <- list(n = n_groups, 
                    x_est = apply(b_samps, 2, mean), 
                    x_err_sd = apply(b_samps, 2, sd)
    )
    
    #fit error model
    fit_errprop_a <- mod_errprop$sample(d_err_a, chains = 4, parallel_chains = 4)
    fit_errprop_b <- mod_errprop$sample(d_err_b, chains = 4, parallel_chains = 4)
    
    #extract samples and inspect results
    a_errp_samps <- data.frame(as_draws_df(fit_errprop_a$draws("x_sigma")))
    b_errp_samps <- data.frame(as_draws_df(fit_errprop_b$draws("x_sigma")))
    
    #take sample means as new estimates
    param_sds_approx <- c(sd_a = mean(a_errp_samps$x_sigma),
                          sd_b = mean(b_errp_samps$x_sigma))
    
  } else {
    
    #compile data
    d_approx <- c(d, param_sds_approx["sd_a"], param_sds_approx["sd_b"])
    # d_frvb <- c(d, param_sds_frvb["sd_a"], param_sds_frvb["sd_b"])
    
    
    if(use_pathfinder){
      fit_approx <- mod_vb$pathfinder(data = d_approx, draws = 5E3)
    } else {
      fit_approx <- mod_vb$variational(data = d_approx, algorithm="meanfield", draws = 5E3)
    }
    
    # fit_frvb <- mod_vb$variational(data = d_frvb, algorithm="fullrank", draws = 5E3, iter = 5E4)
    samps_approx <- as.data.table(data.frame(as_draws_df(fit_approx$draws())))
    # samps_frvb <- as.data.table(data.frame(as_draws_df(fit_frvb$draws())))
    
    #compute the average sample variance across parameters
    scale_param_targets <- c("a", "b")
    scale_param_targets <- setNames(scale_param_targets, paste0("sd_", scale_param_targets))
    source("~/repos/Stan2R/R/functions.R")
    
    #calculate mean variance across parameters
    param_mean_var_approx <- sapply(scale_param_targets, function(scale_param){
      param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_approx)))
      param_vars <- apply(param_samps, 1, var)
      mean(param_vars)
    })
    
    # param_mean_var_frvb <- sapply(scale_param_targets, function(scale_param){
    #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
    #   param_vars <- apply(param_samps, 1, var)
    #   mean(param_vars)
    # })
    
    #and across iterations
    param_var_mean_approx <- sapply(scale_param_targets, function(scale_param){
      param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_approx)))
      param_vars <- apply(param_samps, 2, var)
      mean(param_vars)
    })
    
    # param_var_mean_frvb <- sapply(scale_param_targets, function(scale_param){
    #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
    #   param_vars <- apply(param_samps, 2, var)
    #   mean(param_vars)
    # })
    
    #see file '~/scripts/minor_scripts/postdoc/scale_mixture_variances.R' for derivation + demonstration
    if(accommodate_scale){
      param_sds_approx <- sqrt(param_mean_var_approx - param_var_mean_approx)
    } else {
      param_sds_approx <- sqrt(param_mean_var_approx)  
    }
    # param_sds_frvb <- sqrt(param_mean_var_frvb - param_var_mean_frvb)
    
  }
  
  #remove the sink
  sink()
  
  
  diff_approx <- param_sds_approx - approx_progress[i-1, names(param_sds_approx)]
  # diff_frvb <- param_sds_frvb - unlist(d_frvb[names(param_sds_frvb)])
  
  cat(paste0("diff_approx (", names(param_sds_approx), ") = ", diff_approx, collapse = ", "))
  # cat("\n")
  # cat(paste0("diff_frvb = ", diff_frvb, collapse = ", "))
  
  # mean_diff <- mean(c(abs(diff_approx), abs(diff_frvb)))
  mean_diff <- mean(c(abs(diff_approx)))
  approx_progress <- rbind(approx_progress, c(param_sds_approx, mlp = mean(mean(fit_approx$lp()))))
  # frvb_progress <- rbind(frvb_progress, param_sds_frvb)
  
}

cols <- viridisLite::inferno(101, end = 0.8)
cols_mlp_pos <- approx_progress[,"mlp"] - min(approx_progress[,"mlp"])
cols_mlp_pos <- cols_mlp_pos / max(cols_mlp_pos) * 100 + 1
cols_mlp_pos <- log10(cols_mlp_pos) / log10(max(cols_mlp_pos))
cols_mlp_pos <- ceiling(cols_mlp_pos * 100 + 1E-4)
cols_mlp <- cols[cols_mlp_pos]
best_MLP <- which.max(approx_progress[,"mlp"])
par(mfrow = c(length(scale_param_targets),1))
for(i in 1:length(scale_param_targets)){
  
  #retrieve useful vars
  sp_name <- names(scale_param_targets)[i]
  post_samps <- samps[,sp_name]
  post_mean <- mean(post_samps)
  post_kde <- density(post_samps)
  post_mode <- post_kde$x[which.max(post_kde$y)]
  post_95q <- quantile(post_samps, c(0.025, 0.975))
  
  #find nice breaks
  varlims <- range(c(post_samps, get(sp_name)))
  breaks_interv <- c(diff(post_95q)) / 12
  breaks <- seq(post_95q[1], varlims[1], by = -breaks_interv)
  breaks <- c(rev(breaks[-1]), seq(post_95q[1], varlims[2], by = breaks_interv))
  breaks <- c(breaks[1] - breaks_interv, breaks, breaks[length(breaks)] + breaks_interv)
  
  #plot histogram
  histd <- hist(post_samps, plot = F, breaks = breaks)
  ylims <- c(0, max(histd$counts) * 2.1)
  
  hist(post_samps, main = sp_name, border = "darkslategray", col = "darkslategrey", axes = F,
       breaks = breaks, xlim = c(0, varlims[2]), ylim = ylims)
  axis(2, at = pretty(c(0,max(histd$counts))), labels = pretty(c(0,max(histd$counts))))
  
  #plot posterior summary with true value and multilevel approx valye
  displ_y <- diff(par("usr")[3:4]) / 40
  axis(1, pos = -displ_y * 5)
  mtext(side = 1, padj = displ_y * 3, text = paste0("posterior ", sp_name), xpd = NA)
  displ_y <- diff(par("usr")[3:4]) / 50
  segments(x1 = post_95q[1], 
           x0 = post_95q[2], 
           y0 = par("usr")[3] - displ_y, 
           y1 = par("usr")[3] - displ_y, col = "deepskyblue4", lwd = 2, xpd = NA)
  points(x = post_mode, y = par("usr")[3] - displ_y,
         pch = 19, col = "deepskyblue4", cex = 2, xpd = NA)
  points(x = get(sp_name), y = par("usr")[3] - displ_y, 
         pch = 19, col = 2, cex = 1, xpd = NA)
  abline(v = post_mode, col = adjustcolor("deepskyblue4", 0.5), lty = 3)
  hist(post_samps, main = sp_name, border = "darkslategray", col = "darkslategrey", axes = F,
       breaks = breaks, add = T)
  points(x = param_sds_approx_multilevel[sp_name], y = par("usr")[3] - displ_y,
         pch = 19, col = "darkslateblue", cex = 1, xpd = NA)

  #plot legend
  legend_tr <- mean(approx_progress[,sp_name] < par("usr")[2]/2) > 0.5
  legend_info <- legend(x = ifelse(legend_tr, "topr", "topl"), legend = c("full Bayes MAP", "95% credible interval", "true value", "iter. MF estimate", "multil. MF estimate"), 
         pch = c(19,NA,19,19,19), lty = c(NA,1,NA,1,NA), col = c("deepskyblue4","deepskyblue4",2,1, "darkslateblue"), pt.cex = c(2,NA,1,0.5,1), inset = 0.01)
  
  #plot variation approx
  progress_range <- max(histd$counts) * c(1.1, 2)
  pts_y <- seq(progress_range[2], progress_range[1], length.out = nrow(approx_progress))
  points(approx_progress[,sp_name], pts_y, pch = 19, cex = 0.5)
  lines(approx_progress[,sp_name], pts_y, pch = 19, cex = 0.5)
  text(x = approx_progress[,sp_name], pts_y, pch = 19, cex = 0.5, pos = 4, labels = round(approx_progress[,"mlp"], 1), col = cols_mlp)
  contleg_x <- legend_info$rect$left + legend_info$rect$w * ifelse(legend_tr, 0, 1) + diff(par("usr")[1:2])/30 * ifelse(legend_tr, -1, 1)
  add_continuous_legend(x = contleg_x, 
                        y = legend_info$rect$top, cols_mlp, w = diff(par("usr")[1:2])/50, h = legend_info$rect$h,
                        labels = approx_progress[,"mlp"], positions = 1-(cols_mlp_pos-1)/100, left_below = legend_tr, main = "mlp")
  
  #single out best log prob
  arrowtip_x <- approx_progress[best_MLP,sp_name] + strwidth(round(approx_progress[best_MLP,"mlp"], 1))
  arrows(x1 = arrowtip_x, x0 = arrowtip_x + diff(par("usr")[1:2])/40, y0 = pts_y[best_MLP], y1 = pts_y[best_MLP], length = par("pin")[1]/100)
  text(x = arrowtip_x + diff(par("usr")[1:2])/50, y = pts_y[best_MLP], pos = 4, labels = "best mlp", cex = 0.75)
  
  
}


#hmm, perhaps not quite...
