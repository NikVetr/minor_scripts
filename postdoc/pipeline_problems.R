#functions
softmax <- function(x) exp(x) / sum(exp(x))

#specify high-level parameters
n_rounds <- 5
n_indiv_start <- 5E3
n_indiv_end <- 50
filtration_factors <- exp((log(n_indiv_start) - log(n_indiv_end)) * 
                            softmax(runif(n_rounds, min = 0.1, max = 0.9)))
filtration_probabilities <- 1 / filtration_factors
n_indiv <- c(n_indiv_start, round(n_indiv_start / cumprod(filtration_factors)))
n_dim <- 3
n_groups <- sample(2:6, n_dim, T)
group_membership_probs <- lapply(n_groups, function(ng) softmax(rnorm(ng, mean = 2)))
individual_sd <- 2
individual_roundwise_sd <- 1
mean_group_sd <- 1
roundwise_group_sd <- 0.5
paired_mean_group_sd <- 0.5
paired_roundwise_group_sd <- 0.25

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
  rnorm(n_indiv_start) * individual_roundwise_sd, simplify = F))

#specify main filtration effects
mean_group_effects <- lapply(1:n_dim, function(gi) rnorm(n_groups[gi]) * mean_group_sd)
roundwise_group_deviations <- lapply(1:n_dim, function(gi) 
  do.call(rbind, replicate(n = n_rounds, rnorm(n_groups[gi]) * roundwise_group_sd, simplify = F))
)
names(roundwise_group_deviations) <- paste0("g_", 1:n_dim)

#specify interaction filtration effects
group_pairs_key <- expand.grid(1:n_dim, 1:n_dim)
group_pairs_key <- group_pairs_key[!apply(group_pairs_key, 1, function(x) x[1] == x[2]),]
paired_mean_group_effects <- lapply(1:nrow(group_pairs_key), function(gpi) {
  pairwise_n_groups <- n_groups[unlist(group_pairs_key[gpi,])]
  matrix(rnorm(prod(pairwise_n_groups)) * paired_mean_group_sd, 
         nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2])
})
names(paired_mean_group_effects) <- apply(group_pairs_key, 1, paste0, collapse = "x")
paired_roundwise_group_deviations <- lapply(1:nrow(group_pairs_key), function(gpi) {
  pairwise_n_groups <- n_groups[unlist(group_pairs_key[gpi,])]
  round_effects <- replicate(n = n_rounds, expr = {
    matrix(rnorm(prod(pairwise_n_groups)) * paired_roundwise_group_sd, 
           nrow = pairwise_n_groups[1], ncol = pairwise_n_groups[2])
  }, simplify = F)
  names(round_effects) <- paste0("r_", 1:n_rounds)
  return(round_effects)
})
names(paired_roundwise_group_deviations) <- apply(group_pairs_key, 1, paste0, collapse = "x")

#compute progress probabilities
progress_liabs <- do.call(cbind, lapply(1:n_rounds, function(ri){
  
  #get individual effects for round
  roundwise_vals <- indiv_liab + indiv_liab_x_round[,ri]
  
  #get group effects for round
  roundwise_vals <- roundwise_vals + Reduce("+", lapply(1:n_dim, function(gi1){
    roundwise_groupvals <- mean_group_effects[[gi1]][indiv_data[,paste0("g_", gi1)]] + 
      roundwise_group_deviations[[gi1]][ri, indiv_data[,paste0("g_", gi1)]]
  }))
  
  #get group interaction effects for round
  roundwise_vals <- roundwise_vals + Reduce("+", lapply(1:nrow(group_pairs_key), function(gpi) {
    curr_gis <- paste0("g_", group_pairs_key[gpi,])
    paired_mean_group_effects[[gpi]][cbind(indiv_data[,curr_gis[1]], indiv_data[,curr_gis[2]])] + 
      paired_roundwise_group_deviations[[gpi]][[ri]][cbind(indiv_data[,curr_gis[1]], indiv_data[,curr_gis[2]])]
  }))
    
  
}))

#now transform so propto probabilities of progress
progress_probs <- apply(progress_liabs, 2, softmax)

#simulate selection process
hypergeom <- F #sample filtration from a binomial or a hypergeometric?
ifelse2 <- function(test, out1, out2){if(test){return(out1)}else{return(out2)}}
invlogit <- function(x, par) {
  1 / (1 + exp(-(x + par)))
}
obj <- function(x, par, target) {
  abs(sum(invlogit(x, par)) - target)
}

indiv_data <- cbind(indiv_data, matrix(F, nrow = n_indiv_start, ncol = n_rounds, 
                                       dimnames = list(NULL, paste0("r_", 1:n_rounds))))
for(ri in 1:n_rounds){
  still_here <- ifelse2(ri==1, rep(T, n_indiv_start), indiv_data[,paste0("r_", ri-1)])
  if(hypergeom){
    passing <- sample(which(still_here), size = n_indiv[ri+1], replace = F, 
                        prob = progress_probs[still_here,ri])  
  } else {
    n_still_here <- sum(still_here)
    disp <- optim(par = 0, fn = obj, 
                  x = progress_liabs[still_here, ri], 
                  target = n_indiv[ri+1], 
                  method = "SANN")$par
    # disp
    # plot(-100:100/10, sapply(-100:100/10, function(i) obj(x, i, target)))
    indiv_prob_pass <- invlogit(progress_liabs[still_here, ri], disp)
    passing <- which(rbinom(n = n_still_here, size = 1, prob = indiv_prob_pass) == 1)
  }
  
  indiv_data[passing, paste0("r_", ri)] <- T
}
rounds_progressed <- apply(indiv_data[,paste0("r_", 1:n_rounds)], 1, sum)
n_indiv <- c(n_indiv_start, apply(indiv_data[,paste0("r_", 1:n_rounds)], 2, sum))

#do some quick visualization
plot(indiv_liab, rounds_progressed + rnorm(n = n_indiv_start) * 0.1, pch = 19, col = adjustcolor(1, 0.1))

#now fit the Stan model
dat <- list(
  n0 = n_indiv_start,
  R = n_rounds,
  n = n_indiv,
  G <- n_dim,
  g <- n_groups,
  gi <- indiv_data[,paste0("g_", 1:n_dim)],
  gpk <- group_pairs_key,
  fri <- rounds_progressed + 1
)
