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

#### simulation ####

#specify population and filtration parameters
n_rounds <- 6
n_indiv_start <- 1E3
n_indiv_end <- 50
filtration_factors <- exp((log(n_indiv_start) - log(n_indiv_end)) * 
                            softmax(runif(n_rounds, min = 0.1, max = 0.9)))
n_indiv <- c(n_indiv_start, round(n_indiv_start / cumprod(filtration_factors)))
n_dim <- 5
n_groups <- sample(2:4, n_dim, T)
group_membership_probs <- lapply(n_groups, function(ng) softmax(rnorm(ng, mean = 2)))

individual_sd <- 1
individual_x_round_sd <- 0.25
mean_group_sd <- 1
group_x_round_sd <- 0.25
group_x_group_sd <- 1
group_x_group_x_round_sd <- 0.25

# specify mcmc parameters
nchains <- 4
niter <- 1E3
adapt_delta <- 0.9
max_treedepth <- 10
thin <- 1
init_mcmc <- 0.2
refresh <- 1

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
mean_group_effects <- lapply(1:n_dim, function(gi) rnorm(n_groups[gi]) * mean_group_sd)
group_x_round_deviations <- lapply(1:n_dim, function(gi) 
  do.call(rbind, replicate(n = n_rounds, rnorm(n_groups[gi]) * group_x_round_sd, simplify = F))
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
  matrix(rnorm(prod(pairwise_n_groups)) * group_x_group_sd, 
         nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2])
})
names(group_x_group_effects) <- apply(group_pairs_key, 1, paste0, collapse = "x")
group_x_group_x_round_deviations <- lapply(1:nrow(group_pairs_key), function(gpi) {
  pairwise_n_groups <- n_groups[unlist(group_pairs_key[gpi,])]
  round_effects <- replicate(n = n_rounds, expr = {
    matrix(rnorm(prod(pairwise_n_groups)) * group_x_group_x_round_sd, 
           nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2])
  }, simplify = F)
  names(round_effects) <- paste0("r_", 1:n_rounds)
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
      group_x_group_x_round_deviations[[gpi]][[ri]][cbind(indiv_data[,curr_gis[1]], indiv_data[,curr_gis[2]])]
  }))
    
  
}))

#now transform so propto probabilities of progress
progress_probs <- apply(progress_liabs, 2, softmax)

#simulate selection process
hypergeom <- F #sample filtration from a binomial or a hypergeometric?
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
    passing <- still_here_inds[which(rbinom(n = n_still_here, size = 1, prob = indiv_prob_pass) == 1)]
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

#### inference ####
#now fit the Stan model
dat <- list(
  n0 = n_indiv_start,
  R = n_rounds,
  n = n_indiv,
  G = n_dim,
  g = n_groups,
  gi = as.matrix(indiv_data[,paste0("g_", 1:n_dim)]),
  gpk = as.matrix(group_pairs_key),
  fri = rounds_progressed + 1
)

# compile model
model_path <- paste0("~/scripts/minor_scripts/postdoc/pipeline_problems.stan")
mod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#fit model
fit <- mod$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                  data = dat, parallel_chains = nchains, adapt_delta = adapt_delta,
                  refresh = refresh, max_treedepth = max_treedepth,
                  thin = thin, threads_per_chain = 1,
                  init = init_mcmc)


summ <- fit$summary()
summ[order(summ$ess_bulk), c("variable", 
                             "rhat", "ess_bulk", "ess_tail", 
                             "mean", "sd")]
summ[order(summ$rhat, decreasing = T), c("variable", 
                                         "rhat", "ess_bulk", "ess_tail", 
                                         "mean", "sd")]


#### evaluation ####
source("~/repos/Stan2R/R/functions.R")
params_of_interest <- fit$metadata()$stan_variables
samps <- as.data.table(as_draws_df(fit$draws(params_of_interest)))
samps_by_param <- lapply(setNames(params_of_interest, params_of_interest), function(var_name) 
  subset_samps(var_name = var_name, samps = samps)
)

#main group effects
col_alpha <- 0.8
dim_cols <- adjustcolor(c('#4477AA', '#EE6677', '#228833', 
                          '#CCBB44', '#66CCEE', '#AA3377', 
                          '#BBBBBB')[1:n_dim], col_alpha)
group_inds <- rep(1:n_dim, times = n_groups)
group_names <- paste0("B_group.", group_inds, ".", unlist(sapply(n_groups, function(i) 1:i)))
posterior_mean_group_effects <- apply(samps_by_param$B_group[,..group_names] * 
                                        matrix(samps$B_group_sd, 
                                               nrow = length(samps$B_group_sd), 
                                               ncol = sum(n_groups)), 2, mean)

CI_prop <- 0.5
posterior_CI_group_effects <- apply(samps_by_param$B_group[,..group_names] * 
                                        matrix(samps$B_group_sd, 
                                               nrow = length(samps$B_group_sd), 
                                               ncol = sum(n_groups)), 2, quantile, probs = c(0.5 - CI_prop/2, 
                                                                                             0.5 + CI_prop/2))

#adequacy
plot(unlist(mean_group_effects), posterior_mean_group_effects, 
     ylim = range(posterior_CI_group_effects), pch = 19, col = dim_cols[group_inds], 
     xlab = latex2exp::TeX("\\textbf{true} \\textit{group} effect"),
     ylab = latex2exp::TeX("\\textbf{inferred} \\textit{group} effect"))
for(i in 1:length(unlist(mean_group_effects))){
  segments(x0 = unlist(mean_group_effects)[i], x1 = unlist(mean_group_effects)[i], 
           y0 = posterior_CI_group_effects[1,i], y1 = posterior_CI_group_effects[2,i],
           col = dim_cols[group_inds][i], lwd = 2)
}
abline(h = 0, lty = 2, col = "grey", lwd = 2)
abline(a = 0, b = 1, lwd = 2, col = 2)
legend(x = par("usr")[1] + diff(par("usr")[1:2])/100, y = par("usr")[4] - diff(par("usr")[3:4])/100, pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 1), 
       legend = c("posterior mean", 
                  paste0(round(CI_prop*100, 1), "% CI"), 
                  "y = 0 line", 
                  "1-to-1 line"), 
       lwd = c(NA, 2, 2, 2), box.lty = 2, col = c(adjustcolor(c("grey20", "grey20"), col_alpha), 2, "grey"), bg = "grey99")
legend(x = par("usr")[2] + diff(par("usr")[1:2])/100, y = par("usr")[4] - diff(par("usr")[3:4])/100, pch = 19, 
       legend = paste0("Feature ", 1:n_dim), bty = "n", col = group_cols, xpd = NA)


#calibration
scaled_liabs <- samps_by_param$B_group[,..group_names] * 
  matrix(samps$B_group_sd, 
         nrow = length(samps$B_group_sd), 
         ncol = sum(n_groups))
calibration_probs <- c(0.001, 1:99/100, 0.999)
posterior_quantile <- sapply(1:length(unlist(mean_group_effects)), function(i) 
                                mean(scaled_liabs[,..i] > unlist(mean_group_effects)[i]))
posterior_calibration <- lapply(1:length(unlist(mean_group_effects)), function(i){
  qpost <- quantile(unlist(scaled_liabs[,..i]), probs = calibration_probs)
  bounds <- data.frame(lower = rev(qpost[calibration_probs < 0.5]), 
                  upper = qpost[calibration_probs > 0.5],
                  width = calibration_probs[calibration_probs > 0.5] - rev(calibration_probs[calibration_probs < 0.5]))
  bounds$in_bounds <- unlist(mean_group_effects)[i] > bounds$lower & unlist(mean_group_effects)[i] < bounds$upper
  return(setNames(nm = bounds$width, object = bounds$in_bounds))
})
posterior_calibration <- apply(do.call(rbind, posterior_calibration), 2, mean)

par(mfrow = c(2,1), mar = c(5,5,0,2))
hist(posterior_quantile, breaks = 0:20/20, xlab = "Posterior Quantile", freq = F, main = "")
plot(as.numeric(names(posterior_calibration)), posterior_calibration, type = "l", 
     xlim = c(0,1), ylim = c(0,1), xlab = "Width of CI", ylab = "Coverage Proportion")
abline(0,1,col=2,lty=2)
legend("topleft", lty = 2, col = 2, legend = "1-to-1 line", bty = "n")

hist(unlist(samps_by_param$B_round_sd * 100))
apply(samps_by_param$B_round, 2, function(x) mean(x * unlist(samps_by_param$B_round) * 100))




#interaction (group x group) effects
pairwise_interaction_dims <- lapply(group_x_group_effects, dim)
group_inds <- rep(1:npairs, times = sapply(pairwise_interaction_dims, prod))
gxg_cols <- adjustcolor(c('#CC6677', '#332288', '#DDCC77', 
                          '#117733', '#88CCEE', '#882255', 
                          '#44AA99', '#999933', '#AA4499', 'darkgrey')[1:npairs], col_alpha)

B_group_x_group <- abind::abind(munge_samps("B_group_x_group", samps_by_param$B_group_x_group), along = 4)
B_group_x_group_list <- lapply(1:npairs, function(di){
  B_group_x_group[di,1:pairwise_interaction_dims[[di]][1], 1:pairwise_interaction_dims[[di]][2],]
})

#compute posterior summaries
B_group_x_group_posterior_means <- lapply(1:npairs, function(di) 
  apply(B_group_x_group_list[[di]], c(1, 2), function(x) mean(x * samps$B_group_x_group_sd)))

B_group_x_group_posterior_lowerCI <- lapply(1:npairs, function(di) 
  apply(B_group_x_group_list[[di]], c(1, 2), function(x) 
    quantile(x * samps$B_group_x_group_sd, prob = 0.5 - CI_prop/2)))

B_group_x_group_posterior_upperCI <- lapply(1:npairs, function(di) 
  apply(B_group_x_group_list[[di]], c(1, 2), function(x) 
    quantile(x * samps$B_group_x_group_sd, prob = 0.5 + CI_prop/2)))

par(mar = c(5,5,4,8.5), mfrow = c(1,1))
plot(unlist(group_x_group_effects), unlist(B_group_x_group_posterior_means), 
     ylim = range(c(unlist(B_group_x_group_posterior_lowerCI), 
                    unlist(B_group_x_group_posterior_upperCI))), 
     pch = 19, col = gxg_cols[group_inds], 
     xlab = latex2exp::TeX("\\textbf{true} \\textit{group x group} interaction effect"),
     ylab = latex2exp::TeX("\\textbf{inferred} \\textit{group x group} interaction effect"))
for(i in 1:length(unlist(group_x_group_effects))){
  segments(x0 = unlist(group_x_group_effects)[i], x1 = unlist(group_x_group_effects)[i], 
           y0 = unlist(B_group_x_group_posterior_lowerCI)[i], y1 = unlist(B_group_x_group_posterior_upperCI)[i],
           col = gxg_cols[group_inds][i], lwd = 2)
}
abline(h = 0, lty = 2, col = "grey", lwd = 2)
abline(a = 0, b = 1, lwd = 2, col = 2)
legend(x = par("usr")[1] + diff(par("usr")[1:2])/100, y = par("usr")[4] - diff(par("usr")[3:4])/100, pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 1), 
       legend = c("posterior mean", 
                  paste0(round(CI_prop*100, 1), "% CI"), 
                  "y = 0 line", 
                  "1-to-1 line"), bty = "n",
       lwd = c(NA, 2, 2, 2), box.lty = 2, col = c(adjustcolor(c("grey20", "grey20"), col_alpha), 2, "grey"), bg = "grey99")
legend(x = par("usr")[2] + diff(par("usr")[1:2])/100, y = par("usr")[4] - diff(par("usr")[3:4])/100, pch = 19, 
       legend = sapply(strsplit(names(group_x_group_effects), split = ""), paste0, collapse = " "), 
       bty = "n", col = gxg_cols, xpd = NA)



#calibration
scaled_liabs <- lapply(1:npairs, function(di) 
  apply(B_group_x_group_list[[di]], c(1, 2), function(x) x * samps$B_group_x_group_sd))
calibration_probs <- c(0.001, 1:99/100, 0.999)
posterior_quantile <- lapply(1:length(group_x_group_effects), function(ggi){
  gge <- group_x_group_effects[[ggi]]
  sapply(1:nrow(gge), function(ri){
    sapply(1:ncol(gge), function(ci){
      mean(scaled_liabs[[ggi]][,ri,ci] > gge[ri,ci])
    })
  })
})
  
posterior_calibration <- do.call(rbind, lapply(1:length(group_x_group_effects), function(ggi){
  gge <- group_x_group_effects[[ggi]]
  do.call(rbind, lapply(1:nrow(gge), function(ri){
    do.call(rbind, lapply(1:ncol(gge), function(ci){
      qpost <- quantile(scaled_liabs[[ggi]][,ri,ci], probs = calibration_probs)
      bounds <- data.frame(lower = rev(qpost[calibration_probs < 0.5]), 
                           upper = qpost[calibration_probs > 0.5],
                           width = calibration_probs[calibration_probs > 0.5] - rev(calibration_probs[calibration_probs < 0.5]))
      bounds$in_bounds <- gge[ri,ci] > bounds$lower & gge[ri,ci] < bounds$upper
      return(setNames(nm = bounds$width, object = bounds$in_bounds))
    }))
  }))
}))
posterior_calibration <- apply(posterior_calibration, 2, mean)

par(mfrow = c(2,1), mar = c(5,5,0.5,2))
hist(unlist(posterior_quantile), breaks = 0:20/20, xlab = "Posterior Quantile", freq = F, main = "")
plot(as.numeric(names(posterior_calibration)), posterior_calibration, type = "l", 
     xlim = c(0,1), ylim = c(0,1), xlab = "Width of CI", ylab = "Coverage Proportion")
abline(0,1,col=2,lty=2)
legend("topleft", lty = 2, col = 2, legend = "1-to-1 line", bty = "n")





# individuals
indiv_names <- paste0("B_indiv.", 1:n_indiv_start)
posterior_mean_indiv_effects <- apply(samps_by_param$B_indiv[,..indiv_names] * 
        matrix(samps$B_indiv_sd, 
               nrow = length(samps$B_indiv_sd), 
               ncol = sum(n_indiv_start)), 2, mean)
posterior_CI_indiv_effects <- apply(samps_by_param$B_indiv[,..indiv_names] * 
                                      matrix(samps$B_indiv_sd, 
                                             nrow = length(samps$B_indiv_sd), 
                                             ncol = sum(n_indiv_start)), 2, quantile, probs = c(0.5 - CI_prop/2, 
                                                                                           0.5 + CI_prop/2))

plot(indiv_liab, posterior_mean_indiv_effects)





par(mar = c(5,5,4,8.5), mfrow = c(1,1))
plot(indiv_liab, posterior_mean_indiv_effects, 
     ylim = range(posterior_CI_indiv_effects), 
     pch = 19, col = adjustcolor(1, 0.5), 
     xlab = latex2exp::TeX("\\textbf{true} \\textit{individual} effect"),
     ylab = latex2exp::TeX("\\textbf{inferred} \\textit{individual} effect"))
for(i in 1:length(indiv_liab)){
  segments(x0 = indiv_liab[i], x1 = indiv_liab[i], 
           y0 = posterior_CI_indiv_effects[1,i], y1 = posterior_CI_indiv_effects[2,i],
           col = adjustcolor(1, 0.5), lwd = 2)
}
abline(h = 0, lty = 2, col = "grey", lwd = 2)
abline(a = 0, b = 1, lwd = 2, col = 2)
legend(x = par("usr")[1] + diff(par("usr")[1:2])/100, y = par("usr")[4] - diff(par("usr")[3:4])/100, pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 1), 
       legend = c("posterior mean", 
                  paste0(round(CI_prop*100, 1), "% CI"), 
                  "y = 0 line", 
                  "1-to-1 line"), bty = "n",
       lwd = c(NA, 2, 2, 2), box.lty = 2, col = c(adjustcolor(c("grey20", "grey20"), col_alpha), 2, "grey"), bg = "grey99")


#calibration
scaled_liabs <- samps_by_param$B_indiv[,..indiv_names] * 
  matrix(samps$B_indiv_sd, 
         nrow = length(samps$B_indiv_sd), 
         ncol = n_indiv_start)
calibration_probs <- c(0.001, 1:99/100, 0.999)
posterior_quantile <- sapply(1:n_indiv_start, function(i) 
  mean(scaled_liabs[,..i] > unlist(indiv_liab)[i]))
posterior_calibration <- lapply(1:n_indiv_start, function(i){
  qpost <- quantile(unlist(scaled_liabs[,..i]), probs = calibration_probs)
  bounds <- data.frame(lower = rev(qpost[calibration_probs < 0.5]), 
                       upper = qpost[calibration_probs > 0.5],
                       width = calibration_probs[calibration_probs > 0.5] - rev(calibration_probs[calibration_probs < 0.5]))
  bounds$in_bounds <- unlist(indiv_liab)[i] > bounds$lower & unlist(indiv_liab)[i] < bounds$upper
  return(setNames(nm = bounds$width, object = bounds$in_bounds))
})
posterior_calibration <- apply(do.call(rbind, posterior_calibration), 2, mean)

par(mfrow = c(2,1), mar = c(5,5,0,2))
hist(posterior_quantile, breaks = 0:20/20, xlab = "Posterior Quantile", freq = F, main = "")
plot(as.numeric(names(posterior_calibration)), posterior_calibration, type = "l", 
     xlim = c(0,1), ylim = c(0,1), xlab = "Width of CI", ylab = "Coverage Proportion")
abline(0,1,col=2,lty=2)
legend("topleft", lty = 2, col = 2, legend = "1-to-1 line", bty = "n")


cor(indiv_liab, posterior_mean_indiv_effects)

apply(samps_by_param$B_round, 2, mean)
filtration_factors
