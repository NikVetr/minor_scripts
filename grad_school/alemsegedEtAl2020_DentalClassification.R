library(mvtnorm)
library(coda)

#read in the data
data <- read.csv("data/alemsegedEtAl_DentalClassificationData.csv", header = T)
d <- data[1:102, 3:5]
unknown_spec <- data[103:104,3:5]

#log transform -- lognormals are imo cleaner than idk gammas or w/e
d[,1:2] <- log(d[,1:2])
unknown_spec[,1:2] <- log(unknown_spec[,1:2])

#extract some useful info
spi <- unique(d[,3]) #species indices
nsp <- length(spi) #number of species

#write short convenience functions
cor_to_cor <- function(r){
  cor <- diag(2)
  cor[1,2] <- cor[2,1] <- r
  return(cor)
}

#itialize parameters to reasonable starting values
mus <- sapply(1:2, function(trait) sapply(spi, function(sp) mean(d[d[,3]==sp,trait])))
sds <- sapply(1:2, function(trait) sapply(spi, function(sp) sd(d[d[,3]==sp,trait])))
cors <- sapply(spi, function(sp) cor(d[d[,3]==sp,1:2])[1,2])
unknown_sp_member <- sample(spi, size = 2)

mus_mus <- apply(mus, 2, mean)
mus_sd <- apply(mus, 2, sd)
mus_cor <- cor(mus)
beta_shape = c(1,1)

#specify priors for parameters
mus_mus_exp_rates <- c(1,1) #the exponential distribution's rates for the means of the means
mus_sds_exp_rates <- c(1,1) #the exponential distribution's rates for the variances of the means
cors_beta_shape_exp_rate <- c(0.25,0.25) #you guessed it, another exponential distribution, now for the beta shape parameters 
unknown_sp_member_multinom <- matrix(rep(1/nsp, length(unknown_sp_member) * nsp), ncol = 2, nrow = nsp)

#run MCMC
n_iter <- 1E5
thin <- 1E2
thin_print <- floor(n_iter / 100)
params <- list(mus = array(0, dim = c(nrow(mus), ncol(mus), floor(n_iter/thin))),
               sds = array(0, dim = c(nrow(sds), ncol(sds), floor(n_iter/thin))),
               cors = array(0, dim = c(length(cors), 1, floor(n_iter/thin))),
               sp = array(0, dim = c(2, 1, floor(n_iter/thin))),
               mus_mus = array(0, dim = c(length(mus_mus), 1, floor(n_iter/thin))),
               mus_sd = array(0, dim = c(length(mus_sd), 1, floor(n_iter/thin))),
               mus_cor = array(0, dim = c(nrow(mus_cor), ncol(mus_cor), floor(n_iter/thin))),
               beta_shape = array(0, dim = c(2, 1, floor(n_iter/thin))))

for(i in 1:n_iter){
  if(i %% thin_print == 0){cat(paste0(" ", i / thin_print))}
  param_type_to_poke <- sample(1:length(params), 1)
  if(param_type_to_poke == 1){
    which_mu_to_poke <- sample(1:nrow(mus), 1)
    prop_mu <- mus[which_mu_to_poke,] + rnorm(2, 0, 0.1)
    probs_obs <- sum(dmvnorm(x = d[d[,3] == which_mu_to_poke, 1:2], 
                             mean = prop_mu,
                             sigma = diag(sds[which_mu_to_poke,]) %*% cor_to_cor(cors[which_mu_to_poke]) %*% diag(sds[which_mu_to_poke,]),
                             log = T))
    probs_mu <- dmvnorm(x = prop_mu, mean = mus_mus, sigma = diag(mus_sd) %*% mus_cor %*% diag(mus_sd), log = T)
    prop_probs_p1[which_mu_to_poke] <- probs_mu + probs_obs
    if(exp(prop_probs_p1[which_mu_to_poke] - curr_probs_p1[which_mu_to_poke]) > runif(1,0,1)){
      curr_probs_p1[which_mu_to_poke] <- prop_probs_p1[which_mu_to_poke]
      mus[which_mu_to_poke,] <- prop_mu
    }
  }
  
  if(param_type_to_poke == 2){
    which_sd_to_poke <- sample(1:nrow(sds), 1)
    prop_sd <- sds[which_sd_to_poke,] + rnorm(2, 0, 0.025)
    probs_obs <- sum(dmvnorm(x = d[d[,3] == which_sd_to_poke, 1:2], 
                             mean = mus[which_sd_to_poke,],
                             sigma = diag(prop_sd) %*% cor_to_cor(cors[which_mu_to_poke]) %*% diag(prop_sd),
                             log = T))
    probs_sd <- dmvnorm(x = prop_mu, mean = mus_mus, sigma = diag(mus_sd) %*% mus_cor %*% diag(mus_sd), log = T)
    prop_probs_p2[which_sd_to_poke] <- probs_sd + probs_obs
    if(exp(prop_probs_p2[which_mu_to_poke] - curr_probs_p2[which_mu_to_poke]) > runif(1,0,1)){
      curr_probs_p2[which_mu_to_poke] <- prop_probs_p2[which_mu_to_poke]
      mus[which_mu_to_poke,] <- prop_mu
    }
  }
  
  if(param_type_to_poke == 3){
    
  }
  
  if(param_type_to_poke == 4){
    
  }
  
  if(param_type_to_poke == 5){
    
  }
  
  if(param_type_to_poke == 6){
    
  }
  
  if(param_type_to_poke == 7){
    
  }
  
  if(param_type_to_poke == 8){
    
  }
}

# library(rethinking)
# d <- as.data.frame(d)
# colnames(d)[3] <- "SP"
# pop_model <- ulam(
#   alist(
#     c(MD, BL)[SP] ~ multi_normal(c(mu_MD[SP], mu_BL[SP]), rho[SP], sd[SP]),
#     
#     rho[SP] ~ dbeta(s1, s2),
#     c(s1, s2) ~ dexp(0.25),
#     
#     c(mu_MD, mu_BL)[SP] ~ multi_normal(c(mu_mu_MD, mu_mu_BL), rho_mu, sd_mu),
#     rho_mu ~ dbeta(s1, s2),
#     sd_mu ~ dexp(1),
#     c(mu_mu_MD, mu_mu_BL) ~ dnorm(2, 2),
#     
#     sd[SP] ~ dexp(sd_r),
#     sd_r ~ dhalfnorm(0,1)
#     
#   ) ,
#   data = d, log_lik = T,
#   iter = 1E3, warmup = 1E3, chains = 1, cores = 1
# )
