library(cmdstanr)
library(posterior)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#does iterative EM-flavored EB converge on Full Bayes MCMC?
#or maybe I can call this "maximum marginal likelihood"

#simulate data
n_obs <- 10 #per group
n_groups <- 20
sd_a <- 1.75
sd_b <- 4
sd_e <- 5
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

mod <- cmdstan_model(write_stan_file(stan_program))
mod_vb <- cmdstan_model(write_stan_file(stan_program_vb))

#first fit w/ full bayes
d <- list(n_obs = n_obs, n_groups = n_groups, g = g, x = x, y = y)
fit <- mod$sample(d, chains = 4, parallel_chains = 4)
summ <- fit$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]
samps <- data.frame(as_draws_df(fit$draws()))
apply(samps[,grepl("sd_", colnames(samps))], 2, mean)

#now fit w/ VB once using wide priors
d_vb <- c(d, sd_a = 1E3, sd_b = 1E3)
fit_mfvb <- mod_vb$variational(data = d_vb, algorithm="meanfield", draws = 5E3)
fit_mfvb_multilevel <- mod$variational(data = d, algorithm="meanfield", draws = 5E3)

# fit_frvb <- mod_vb$variational(data = d_vb, algorithm="fullrank", draws = 5E3, iter = 5E4)
samps_mfvb <- as.data.table(data.frame(as_draws_df(fit_mfvb$draws())))
# samps_frvb <- as.data.table(data.frame(as_draws_df(fit_frvb$draws())))

#compute the average sample variance across parameters
scale_param_targets <- c("a", "b")
scale_param_targets <- setNames(scale_param_targets, paste0("sd_", scale_param_targets))
source("~/repos/Stan2R/R/functions.R")

#calculate mean variance across parameters
param_mean_var_mfvb <- sapply(scale_param_targets, function(scale_param){
  param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
  param_vars <- apply(param_samps, 1, var)
  mean(param_vars)
})

# param_mean_var_frvb <- sapply(scale_param_targets, function(scale_param){
#   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
#   param_vars <- apply(param_samps, 1, var)
#   mean(param_vars)
# })

#and across iterations
param_var_mean_mfvb <- sapply(scale_param_targets, function(scale_param){
  param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
  param_vars <- apply(param_samps, 2, var)
  mean(param_vars)
})

# param_var_mean_frvb <- sapply(scale_param_targets, function(scale_param){
#   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
#   param_vars <- apply(param_samps, 2, var)
#   mean(param_vars)
# })

#see file '~/scripts/minor_scripts/postdoc/scale_mixture_variances.R' for derivation + demonstration
param_sds_mfvb <- sqrt(param_mean_var_mfvb - param_var_mean_mfvb)
# param_sds_mfvb <- sqrt(param_mean_var_mfvb)
# param_sds_frvb <- sqrt(param_mean_var_frvb - param_var_mean_frvb)

apply(samps[,names(scale_param_targets)], 2, mean) #"true" values
param_sds_mfvb
# param_sds_frvb

#compare to fitting multilevel model with VB
samps_mfvb_multilevel <- (data.frame(as_draws_df(fit_mfvb_multilevel$draws())))
param_sds_mfvb_multilevel <- apply(samps_mfvb_multilevel[,names(scale_param_targets)], 2, mean)
#now iterate over this procedure
n_iter <- 10
eps <- 1E-3
i <- 1
mean_diff <- Inf
mfvb_progress <- t(data.frame(c(param_sds_mfvb, mlp = mean(mean(fit_mfvb$lp())))))
# frvb_progress <- t(data.frame(param_sds_frvb))

while(i <= n_iter & mean_diff > eps){
  
  i <- i + 1
  print(i)
  
  #compile data
  d_mfvb <- c(d, param_sds_mfvb["sd_a"], param_sds_mfvb["sd_b"])
  # d_frvb <- c(d, param_sds_frvb["sd_a"], param_sds_frvb["sd_b"])
  
  #fit the model
  sink(tempfile())
  fit_mfvb <- mod_vb$variational(data = d_mfvb, algorithm="meanfield", draws = 5E3)
  sink()
  # fit_frvb <- mod_vb$variational(data = d_frvb, algorithm="fullrank", draws = 5E3, iter = 5E4)
  samps_mfvb <- as.data.table(data.frame(as_draws_df(fit_mfvb$draws())))
  # samps_frvb <- as.data.table(data.frame(as_draws_df(fit_frvb$draws())))
  
  #compute the average sample variance across parameters
  scale_param_targets <- c("a", "b")
  scale_param_targets <- setNames(scale_param_targets, paste0("sd_", scale_param_targets))
  source("~/repos/Stan2R/R/functions.R")
  
  #calculate mean variance across parameters
  param_mean_var_mfvb <- sapply(scale_param_targets, function(scale_param){
    param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
    param_vars <- apply(param_samps, 1, var)
    mean(param_vars)
  })
  
  # param_mean_var_frvb <- sapply(scale_param_targets, function(scale_param){
  #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
  #   param_vars <- apply(param_samps, 1, var)
  #   mean(param_vars)
  # })
  
  #and across iterations
  param_var_mean_mfvb <- sapply(scale_param_targets, function(scale_param){
    param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
    param_vars <- apply(param_samps, 2, var)
    mean(param_vars)
  })
  
  # param_var_mean_frvb <- sapply(scale_param_targets, function(scale_param){
  #   param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_frvb)))
  #   param_vars <- apply(param_samps, 2, var)
  #   mean(param_vars)
  # })
  
  #see file '~/scripts/minor_scripts/postdoc/scale_mixture_variances.R' for derivation + demonstration
  param_sds_mfvb <- sqrt(param_mean_var_mfvb - param_var_mean_mfvb)
  param_sds_mfvb <- sqrt(param_mean_var_mfvb)
  # param_sds_frvb <- sqrt(param_mean_var_frvb - param_var_mean_frvb)
  
  diff_mfvb <- param_sds_mfvb - unlist(d_mfvb[names(param_sds_mfvb)])
  # diff_frvb <- param_sds_frvb - unlist(d_frvb[names(param_sds_frvb)])
  
  cat(paste0("diff_mfvb = ", diff_mfvb, collapse = ", "))
  # cat("\n")
  # cat(paste0("diff_frvb = ", diff_frvb, collapse = ", "))
  
  # mean_diff <- mean(c(abs(diff_mfvb), abs(diff_frvb)))
  mean_diff <- mean(c(abs(diff_mfvb)))
  mfvb_progress <- rbind(mfvb_progress, c(param_sds_mfvb, mlp = mean(mean(fit_mfvb$lp()))))
  # frvb_progress <- rbind(frvb_progress, param_sds_frvb)
  
}


cols <- viridisLite::inferno(101, end = 0.8)
cols_mlp_pos <- mfvb_progress[,"mlp"] - min(mfvb_progress[,"mlp"])
cols_mlp_pos <- cols_mlp_pos / max(cols_mlp_pos) * 100 + 1
cols_mlp_pos <- log10(cols_mlp_pos) / log10(max(cols_mlp_pos))
cols_mlp_pos <- ceiling(cols_mlp_pos * 100 + 1E-4)
cols_mlp <- cols[cols_mlp_pos]
best_MLP <- which.max(mfvb_progress[,"mlp"])
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
  
  #plot posterior summary with true value and multilevel MFVB valye
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
  points(x = param_sds_mfvb_multilevel[sp_name], y = par("usr")[3] - displ_y,
         pch = 19, col = "darkslateblue", cex = 1, xpd = NA)

  #plot legend
  legend_tr <- mean(mfvb_progress[,sp_name] < par("usr")[2]/2) > 0.5
  legend_info <- legend(x = ifelse(legend_tr, "topr", "topl"), legend = c("full Bayes MAP", "95% credible interval", "true value", "iter. MF estimate", "multil. MF estimate"), 
         pch = c(19,NA,19,19,19), lty = c(NA,1,NA,1,NA), col = c("deepskyblue4","deepskyblue4",2,1, "darkslateblue"), pt.cex = c(2,NA,1,0.5,1), inset = 0.01)
  
  #plot variation approx
  progress_range <- max(histd$counts) * c(1.1, 2)
  pts_y <- seq(progress_range[2], progress_range[1], length.out = nrow(mfvb_progress))
  points(mfvb_progress[,sp_name], pts_y, pch = 19, cex = 0.5)
  lines(mfvb_progress[,sp_name], pts_y, pch = 19, cex = 0.5)
  text(x = mfvb_progress[,sp_name], pts_y, pch = 19, cex = 0.5, pos = 4, labels = round(mfvb_progress[,"mlp"], 1), col = cols_mlp)
  contleg_x <- legend_info$rect$left + legend_info$rect$w * ifelse(legend_tr, 0, 1) + diff(par("usr")[1:2])/30 * ifelse(legend_tr, -1, 1)
  add_continuous_legend(x = contleg_x, 
                        y = legend_info$rect$top, cols_mlp, w = diff(par("usr")[1:2])/50, h = legend_info$rect$h,
                        labels = mfvb_progress[,"mlp"], positions = 1-(cols_mlp_pos-1)/100, left_below = legend_tr, main = "mlp")
  
  #single out best log prob
  arrowtip_x <- mfvb_progress[best_MLP,sp_name] + strwidth(round(mfvb_progress[best_MLP,"mlp"], 1))
  arrows(x1 = arrowtip_x, x0 = arrowtip_x + diff(par("usr")[1:2])/40, y0 = pts_y[best_MLP], y1 = pts_y[best_MLP], length = par("pin")[1]/100)
  text(x = arrowtip_x + diff(par("usr")[1:2])/50, y = pts_y[best_MLP], pos = 4, labels = "best mlp", cex = 0.75)
  
  
}


#hmm, perhaps not quite...
