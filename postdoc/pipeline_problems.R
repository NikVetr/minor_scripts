library(cmdstanr)
library(posterior)
library(data.table)

#### functions ####
softmax <- function(x){exp(x) / sum(exp(x))}
ifelse2 <- function(test, out1, out2){if(test){return(out1)}else{return(out2)}}
invlogit <- function(x, par) {
  1 / (1 + exp(-(x + par)))
}
obj <- function(x, par, target) {
  abs(sum(invlogit(x, par)) - target)
}
mean_center_vec <- function(x){x - mean(x)}
mean_center_mat <- function(x){
  x <- x - rowMeans(x) #center by rows
  x <- t(t(x) - colMeans(x)) #center by columns
  return(x)
}
mean_center <- function(x, na.rm = FALSE) {
  for (d in seq_along(dim(x))) {
    x <- sweep(x, d, apply(x, d, mean, na.rm = na.rm), FUN = "-")
  }
  x
}

ZerosumBasis <- function(n) {
  v <- rep(1, n)
  M <- diag(n) - (1/n) * (v %*% t(v))
  Q_full <- qr.Q(qr(M), complete = FALSE)
  Q_part <- Q_full[, 1:(n - 1)]
  Q_part
}

plot_posterior_diagnostics <- function(truth,
                                       samples,
                                       CI_prop = 0.9,
                                       name     = NULL,
                                       col      = NULL,
                                       col_names = NULL,
                                       alpha    = 0.8,
                                       new_plot = T)
{
  ## — input checks 
  stopifnot(is.numeric(truth),
            length(dim(samples)) >= 2)
  
  if (is.null(name)) name <- "parameter"
  p <- length(truth)
  if (NCOL(samples) != p)
    stop("ncol(samples) must equal length(truth)")
  
  ## colour handling -
  if (is.null(col)) {
    base_cols <- c("#4477AA", "#EE6677", "#228833",
                   "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    col       <- grDevices::adjustcolor(rep(base_cols, length.out = p),
                                        alpha.f = alpha)
  } else {
    col <- grDevices::adjustcolor(rep(col, length.out = p),
                                  alpha.f = alpha)
  }
  
  ## summaries -
  post_mean <- colMeans(samples)
  ci_bounds <- apply(samples, 2, stats::quantile,
                     probs = c(0.5 - CI_prop/2,
                               0.5 + CI_prop/2))
  
  ## histogram of posterior quantiles -
  post_q <- vapply(seq_len(p),
                   \(j) mean(samples[, j] > truth[j]),
                   numeric(1))
  
  ## calibration curve -
  calib_probs   <- c(0.001, seq(0.01, 0.99, 0.01), 0.999)
  widths        <- calib_probs[calib_probs > .5] -
    rev(calib_probs[calib_probs < .5])
  cover_mat <- vapply(seq_len(p), \(j) {
    qpost <- stats::quantile(samples[, j], probs = calib_probs)
    lower <- rev(qpost[calib_probs < 0.5])
    upper <-       qpost[calib_probs > 0.5]
    truth[j] > lower & truth[j] < upper
  }, logical(length(widths)))
  cover_prop <- rowMeans(cover_mat)
  
  ## — plotting 
  if(new_plot){
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE),
                     heights = c(2, 1, 1))  
  }
  
  ## 1. scatter: posterior mean vs. truth + CI bars -
  graphics::par(mar = c(5, 5, 2, 7))
  graphics::plot(truth, post_mean,
                 pch = 19, col = col,
                 xlab = bquote(bold(true)~italic(.(name))),
                 ylab = bquote(bold(inferred)~italic(.(name))),
                 ylim = range(ci_bounds))
  segments(x0 = truth, x1 = truth,
           y0 = ci_bounds[1, ], y1 = ci_bounds[2, ],
           col = col, lwd = 2)
  abline(h = 0,           lty = 2, col = "grey50", lwd = 2)
  abline(a = 0, b = 1,    lty = 2, col = "red2",   lwd = 2)
  
  graphics::legend("topleft",
                   legend = c("posterior mean",
                              sprintf("%.1f%% CI", 100*CI_prop),
                              "y = 0 line",
                              "1-to-1 line"),
                   pch    = c(19, NA, NA, NA),
                   lty    = c(NA, 1, 2, 2),
                   lwd    = c(NA, 2, 2, 2),
                   col    = c("grey30", "grey30", "red2", "grey50"),
                   box.lty = 2, bg = "grey95", cex = 0.8)
  
  ## optional colour legend (right margin) 
  if(!is.null(col_names)){
    usr <- graphics::par("usr")
    graphics::legend(x = usr[2] + diff(usr[1:2])*0.02,
                     y = usr[4],
                     legend = names(col_names),
                     pch    = 19,
                     col    = col_names,
                     xpd    = NA, bty = "n", cex = 0.7)  
  }
  
  ## 2. posterior-quantile histogram 
  graphics::par(mar = c(5, 5, 1, 2))
  graphics::hist(post_q,
                 breaks = seq(0, 1, 0.05),
                 freq   = FALSE,
                 xlab   = "Posterior quantile",
                 main   = "")
  
  ## 3. calibration curve -
  graphics::plot(widths, cover_prop, type = "l",
                 xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "Credible-interval width",
                 ylab = "Coverage proportion")
  abline(0, 1, col = "red2", lty = 2)
  graphics::legend("topleft", lty = 2, col = "red2",
                   legend = "1-to-1 line", bty = "n")
  
  invisible(list(mean   = post_mean,
                 ci     = ci_bounds,
                 q_hist = post_q,
                 calib  = data.frame(width = widths,
                                     coverage = cover_prop)))
}

#### simulation ####

#specify population and filtration parameters
n_rounds <- 3
n_indiv_start <- 1E3
n_indiv_end <- 10
filtration_factors <- exp((log(n_indiv_start) - log(n_indiv_end)) * 
                            softmax(runif(n_rounds, min = 0.1, max = 0.9)))
n_indiv <- c(n_indiv_start, round(n_indiv_start / cumprod(filtration_factors)))
n_dim <- 3
n_groups <- sample(3:5, n_dim, T)
group_membership_probs <- lapply(n_groups, function(ng) softmax(rnorm(ng, mean = 2)))

individual_sd <- 0
individual_x_round_sd <- 0.25
mean_group_sd <- 0.5
group_x_round_sd <- 0.25
group_x_group_sd <- 0.5
group_x_group_x_round_sd <- 0.25

# specify mcmc parameters
nchains <- 4
niter <- 1E3
adapt_delta <- 0.9
max_treedepth <- 10
thin <- 1
init_mcmc <- 0.2
refresh <- 10

#simulate individuals
indiv_data <- data.frame(do.call(cbind, lapply(1:n_dim, function(gi) 
  sample(size = n_indiv_start, 
         x = 1:n_groups[gi], 
         replace = T, 
         prob = group_membership_probs[[gi]])
)))
colnames(indiv_data) <- paste0("g_", 1:n_dim)
indiv_liab <- rnorm(n_indiv_start) * individual_sd
indiv_liab_x_round <- do.call(cbind, replicate(n = n_rounds, 
  rnorm(n_indiv_start) * individual_x_round_sd, simplify = F))

#specify main filtration effects
mean_group_effects <- lapply(1:n_dim, function(gi) mean_center(rnorm(n_groups[gi]) * mean_group_sd))
group_x_round_deviations <- lapply(1:n_dim, function(gi) 
  mean_center(do.call(rbind, replicate(n = n_rounds, rnorm(n_groups[gi]) * group_x_round_sd, simplify = F)))
)
names(group_x_round_deviations) <- paste0("g_", 1:n_dim)

#specify interaction filtration effects
npairs <- choose(n_dim, 2)
group_pairs_key <- expand.grid(1:n_dim, 1:n_dim)
group_pairs_key <- group_pairs_key[!apply(group_pairs_key, 1, function(x) x[1] == x[2]),]
group_pairs_key <- t(apply(group_pairs_key, 1, sort))
group_pairs_key <- group_pairs_key[!duplicated(group_pairs_key),]
group_x_group_effects <- lapply(1:nrow(group_pairs_key), function(gpi) {
  pairwise_n_groups <- n_groups[unlist(group_pairs_key[gpi,])]
  mean_center(matrix(rnorm(prod(pairwise_n_groups)) * group_x_group_sd, 
         nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2]))
})
names(group_x_group_effects) <- apply(group_pairs_key, 1, paste0, collapse = "x")
group_x_group_x_round_deviations <- lapply(1:nrow(group_pairs_key), function(gpi) {
  group_sets <- unlist(group_pairs_key[gpi,])
  pairwise_n_groups <- n_groups[group_sets]
  round_effects <- replicate(n = n_rounds, expr = {
    mat <- matrix(rnorm(prod(pairwise_n_groups)) * group_x_group_x_round_sd, 
           nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2])
    rownames(mat) <- paste0("gr", ".", group_sets[1], ".", 1:pairwise_n_groups[1])
    colnames(mat) <- paste0("gr", ".", group_sets[2], ".", 1:pairwise_n_groups[2])
    return(mat)
  }, simplify = F)
  names(round_effects) <- paste0("r_", 1:n_rounds)
  round_effects <- mean_center(abind::abind(round_effects, along = 3))
  return(round_effects)
})
names(group_x_group_x_round_deviations) <- apply(group_pairs_key, 1, paste0, collapse = "x")

#compute progress probabilities
progress_liabs <- do.call(cbind, lapply(1:n_rounds, function(ri){
  
  #get individual effects for round
  x_round_vals <- indiv_liab + indiv_liab_x_round[,ri]
  
  #get group effects for round
  x_round_vals <- x_round_vals + Reduce("+", lapply(1:n_dim, function(gi1){
    group_x_roundvals <- mean_group_effects[[gi1]][indiv_data[,paste0("g_", gi1)]] + 
      group_x_round_deviations[[gi1]][ri, indiv_data[,paste0("g_", gi1)]]
  }))
  
  #get group interaction effects for round
  x_round_vals <- x_round_vals + Reduce("+", lapply(1:nrow(group_pairs_key), function(gpi) {
    curr_gis <- paste0("g_", group_pairs_key[gpi,])
    group_x_group_effects[[gpi]][cbind(indiv_data[,curr_gis[1]], indiv_data[,curr_gis[2]])] + 
      group_x_group_x_round_deviations[[gpi]][cbind(indiv_data[,curr_gis[1]], 
                                                    indiv_data[,curr_gis[2]],
                                                    ri)]
  }))
    
  
}))

#now transform so propto probabilities of progress?
progress_probs <- apply(progress_liabs, 2, softmax)

#simulate selection process
hypergeom <- F #sample filtration from a binomial or a hypergeometric?-
indiv_data <- cbind(indiv_data, matrix(F, nrow = n_indiv_start, ncol = n_rounds, 
                                       dimnames = list(NULL, paste0("r_", 1:n_rounds))))
round_disps <- rep(NA, n_rounds)
for(ri in 1:n_rounds){
  still_here <- ifelse2(ri==1, rep(T, n_indiv_start), indiv_data[,paste0("r_", ri-1)])
  still_here_inds <- which(still_here)
  if(hypergeom){
    passing <- sample(which(still_here), size = n_indiv[ri+1], replace = F, 
                        prob = progress_probs[still_here,ri])  
  } else {
    n_still_here <- sum(still_here)
    disp <- optim(par = 0, fn = obj, 
                  x = progress_liabs[still_here, ri], 
                  target = n_indiv[ri+1], 
                  method = "SANN")$par
    round_disps[ri] <- disp
    indiv_prob_pass <- invlogit(progress_liabs[still_here, ri], disp)
    passing <- still_here_inds[which(rbinom(n = n_still_here, 
                                            size = 1, 
                                            prob = indiv_prob_pass) == 1)]
  }
  
  indiv_data[passing, paste0("r_", ri)] <- T #oops, the passing inds are for everyone, not just the still here subset
}
rounds_progressed <- apply(indiv_data[,paste0("r_", 1:n_rounds)], 1, sum)
n_indiv <- c(n_indiv_start, apply(indiv_data[,paste0("r_", 1:n_rounds)], 2, sum))

#do some quick visualization
par(mar = c(5,5,4,8.5))
plot(indiv_liab, rounds_progressed + rnorm(n = n_indiv_start) * 0.1, pch = 19, 
     col = sapply((1:(n_rounds+1)/(n_rounds+1)*0.9)^1.5, adjustcolor, col = 1)[rounds_progressed+1],
     xlab = "individual liability (log-odds deviation to overall success)",
     ylab = "final round of filtration (jittered)")
text(x = par("usr")[2], y = 0:n_rounds, xpd = NA, pos = 4, labels = paste(n_indiv, "individuals left"))

#construct Q matrices
Q_R  <- ZerosumBasis(n_rounds)
Qg1  <- ZerosumBasis(n_groups[1])
Qg2  <- ZerosumBasis(n_groups[2])
Qg3  <- ZerosumBasis(n_groups[3])

#### inference ####

#specify the data object to pass to Stan
dat <- list(
  n0 = n_indiv_start,
  R = n_rounds,
  n = n_indiv,
  G = n_dim,
  g = n_groups,
  gi = as.matrix(indiv_data[,paste0("g_", 1:n_dim)]),
  gpk = as.matrix(group_pairs_key),
  fri = rounds_progressed + 1,
  Q_R = Q_R,
  Q_g1 = Qg1,
  Q_g2 = Qg2,
  Q_g3 = Qg3
)

# compile model

model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#fit model

#s2z implementation
model_path_s2z <- paste0("~/scripts/minor_scripts/postdoc/iterative_filtration-s2z.stan")
mod_s2z <- cmdstan_model(model_path_s2z, cpp_options = list(stan_threads = TRUE))
fit_s2z <- mod$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                      data = dat, parallel_chains = nchains, adapt_delta = adapt_delta,
                      refresh = refresh, max_treedepth = max_treedepth,
                      thin = thin, threads_per_chain = 1,
                      init = init_mcmc)

#pop implementation
model_path_pop <- paste0("~/scripts/minor_scripts/postdoc/iterative_filtration.stan")
mod_pop <- cmdstan_model(model_path_pop, cpp_options = list(stan_threads = TRUE))
fit_pop <- mod$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                      data = dat, parallel_chains = nchains, adapt_delta = adapt_delta,
                      refresh = refresh, max_treedepth = max_treedepth,
                      thin = thin, threads_per_chain = 1,
                      init = init_mcmc)

#check mcmc diags
summ_s2z <- fit_s2z$summary()
summ_s2z[order(summ_s2z$ess_bulk), c("variable", 
                                     "rhat", "ess_bulk", "ess_tail", 
                                     "mean", "sd")]
summ_s2z[order(summ_s2z$rhat, decreasing = T), c("variable", 
                                                 "rhat", "ess_bulk", "ess_tail", 
                                                 "mean", "sd")]

summ_pop <- fit_pop$summary()
summ_pop[order(summ_pop$ess_bulk), c("variable", 
                                     "rhat", "ess_bulk", "ess_tail", 
                                     "mean", "sd")]
summ_pop[order(summ_pop$rhat, decreasing = T), c("variable", 
                                                 "rhat", "ess_bulk", "ess_tail", 
                                                 "mean", "sd")]

#recover samples
params_of_interest <- intersect(fit_s2z$metadata()$stan_variables,
                                fit_pop$metadata()$stan_variables)
params_of_interest <- params_of_interest[startsWith(params_of_interest, "B_")]

samps_s2z <- as.data.table(as_draws_df(fit_s2z$draws(params_of_interest)))
samps_s2z <- lapply(params_of_interest, function(var_name) 
  subset_samps(var_name = var_name, samps = samps_s2z)
)
samps_s2z <- as.data.frame(do.call(cbind, samps_s2z))

samps_pop <- as.data.table(as_draws_df(fit_pop$draws(params_of_interest)))
samps_pop <- lapply(params_of_interest, function(var_name) 
  subset_samps(var_name = var_name, samps = samps_pop)
)
samps_pop <- as.data.frame(do.call(cbind, samps_pop))

all_samps <- list(s2z = samps_s2z,
                  pop = samps_pop)

#### evaluation ####

source("~/repos/Stan2R/R/functions.R")
params_of_interest <- fit$metadata()$stan_variables

#set up main plot params
graphics::layout(matrix(1:6, 3, 2, byrow = F),
                 heights = c(2, 1, 1))
types <- names(all_samps)

#main group effects
col_alpha <- 0.8
dim_cols <- adjustcolor(c('#4477AA', '#EE6677', '#228833', 
                          '#CCBB44', '#66CCEE', '#AA3377', 
                          '#BBBBBB')[1:n_dim], col_alpha)
group_inds <- rep(1:n_dim, times = n_groups)
group_names <- paste0("B_g", group_inds, ".", unlist(sapply(n_groups, function(i) 1:i)))
param_names <- do.call(rbind, strsplit(group_names, "\\."))[,1]
g_cols <- setNames(dim_cols, unique(param_names))

par(oma = c(0,0,2,0))
pcors <- list()
for(type in types){
  samps <- all_samps[[type]]
  plot_posterior_diagnostics(truth = setNames(unlist(mean_group_effects), group_names), 
                             samples = samps[,group_names], 
                             col = g_cols[param_names], 
                             name = "group effect",
                             col_names = g_cols, 
                             new_plot = F)
  pcors[[type]] <- cor(unlist(mean_group_effects), colMeans(samps[,group_names]))
}
title_string <- paste0(types, " (corr = ", round(unlist(pcors[types]), 2), ")", collapse = ", ")
mtext(text = title_string, side = 3, outer = T, cex = 1.25, line = -0.5)

#interaction (group x group) effects
B_group_x_group <- lapply(group_x_group_effects, function(gxg){
  inds <- as.matrix(expand.grid(1:nrow(gxg), 1:ncol(gxg)))
  out <- gxg[inds]
  names(out) <- apply(inds, 1, paste0, collapse = ".")
  out
})
setNm <- function(x) setNames(x,x)
B_group_x_group <- lapply(setNm(names(B_group_x_group)), function(gxgn){
  out <- B_group_x_group[[gxgn]]
  names(out) <- paste0("B_g", gxgn, ".", names(out))
  out
})
B_group_x_group <- unlist(setNames(B_group_x_group, NULL))

gxg_cols <- adjustcolor(c('#CC6677', '#332288', '#DDCC77', 
                          '#117733', '#88CCEE', '#882255', 
                          '#44AA99', '#999933', '#AA4499', 'darkgrey')[1:npairs], col_alpha)
param_names <- do.call(rbind, strsplit(names(B_group_x_group), "\\."))[,1]
gxg_cols <- setNames(gxg_cols, unique(param_names))

for(type in types){
  samps <- all_samps[[type]]
  plot_posterior_diagnostics(truth = B_group_x_group, 
                             samples = samps[,names(B_group_x_group)], 
                             col = gxg_cols[param_names], 
                             name = "group x group effect",
                             col_names = gxg_cols, 
                             new_plot = F)
  
  pcors[[type]] <- cor(B_group_x_group, colMeans(samps[,names(B_group_x_group)]))
}
title_string <- paste0(types, " (corr = ", round(unlist(pcors[types]), 2), ")", collapse = ", ")
mtext(text = title_string, side = 3, outer = T, cex = 1.25, line = -0.5)

#round effects?
for(type in types){
  samps <- all_samps[[type]]
  plot_posterior_diagnostics(truth = round_disps, 
                             samples = samps[,paste0("B_round.", 1:n_rounds)], 
                             name = "round effect", 
                             new_plot = F)
}
title_string <- paste0(types, collapse = ", ")
mtext(text = title_string, side = 3, outer = T, cex = 1.25, line = -0.5)


#interaction (group x group x round) effects
B_group_x_group_x_round <- unlist(lapply(names(group_x_group_x_round_deviations), function(gxgn){
  gxg <- group_x_group_x_round_deviations[[gxgn]]
  gxg_out <- lapply(1:n_rounds, function(ri){
    inds <- cbind(as.matrix(expand.grid(1:nrow(gxg), 1:ncol(gxg))), ri)
    out <- gxg[inds]
    names(out) <- apply(inds[,c(3,1,2)], 1, paste0, collapse = ".")
    out
  })
  gxg_out <- unlist(gxg_out)
  names(gxg_out) <- paste0("B_g", gxgn, "_x_round.", names(gxg_out))
  gxg_out
}))
B_group_x_group_x_round
gxgxr_cols <- adjustcolor(c('#CC6677', '#332288', '#DDCC77', 
                          '#117733', '#88CCEE', '#882255', 
                          '#44AA99', '#999933', '#AA4499', 'darkgrey')[1:n_rounds], col_alpha)
param_names <- do.call(rbind, strsplit(names(B_group_x_group_x_round), "\\."))[,2]
gxgxr_cols <- setNames(gxg_cols, unique(param_names))

for(type in types){
  samps <- all_samps[[type]]
  plot_posterior_diagnostics(truth = B_group_x_group_x_round, 
                             samples = samps[,names(B_group_x_group_x_round)], 
                             col = gxgxr_cols[param_names], 
                             name = "group x group x round effect",
                             col_names = gxgxr_cols, 
                             new_plot = F)
  pcors[[type]] <- cor(B_group_x_group_x_round, colMeans(samps[,names(B_group_x_group_x_round)]))
}
title_string <- paste0(types, " (corr = ", round(unlist(pcors[types]), 2), ")", collapse = ", ")
mtext(text = title_string, side = 3, outer = T, cex = 1.25, line = -0.5)

#compare posterior variances and ess / s?
specific_params_of_interest <- intersect(colnames(samps_pop), colnames(samps_s2z))
specific_params_of_interest_roots <- sapply(strsplit(specific_params_of_interest, "\\."), head, 1)
samps_s2z


#### end ####