#### narrow-sense heritability calculations ####
set.seed(123)

#SNP-based heritability
N <- 1E7
minor_allele_freq <- 0.2
genotype_probs <- c((1-minor_allele_freq)^2, 2*minor_allele_freq*(1-minor_allele_freq), minor_allele_freq^2) #at hw
snps <- sample(c(0,1,2), size = N, replace = T, prob = genotype_probs)
snps_std <- (snps - mean(snps)) / sd(snps)
beta <- 3 #effect of allele on trait
alpha <- 7 #intercept when no alleles present
env_sd <- 1.5 #standard deviation of environmental effect
error <- rnorm(N, mean = 0, sd = env_sd)
traits <- alpha + beta * snps + error
traits_std <- (traits - mean(traits)) / sd(traits) #mean center and rescale to unit variance
# plot(snps, traits)
summary(lm(traits ~ snps))$r.squared
summary(lm(traits_std ~ snps_std))$coefficients["snps_std","Estimate"]^2
# var(traits) - (var(snps * beta) + env_sd^2) #confirm additivity of variances

#plain ol' heritability, per the definition
var(snps * beta) / var(traits)

#midparent regression on dizygotic twins / full sibs raised in same environment
parent_pairs <- data.frame(p1 = sample(1:N, size = N/2))
parent_pairs$p2 <- sample(setdiff(1:N, parent_pairs$p1))
n_child_per_parents <- 2
snps_child <- replicate(n_child_per_parents, 
                        rbinom(n = N/2, size = 1, prob = snps[parent_pairs$p1]/2) + #allele from p1
                          rbinom(n = N/2, size = 1, prob = snps[parent_pairs$p2]/2)) #allele from p2
times_beta_plus_intercept_noise <- function(x){return(x*beta + alpha + rnorm(1, mean = 0, sd = env_sd))} 
traits_child <- t(apply(snps_child, 1, times_beta_plus_intercept_noise)) #assumes indep. env. effects across parents / children
midparent_traits <- (traits[parent_pairs$p1] + traits[parent_pairs$p2])/2
summary(lm(as.vector(traits_child) ~ rep(midparent_traits, n_child_per_parents)))$
  coefficients["rep(midparent_traits, n_child_per_parents)","Estimate"]

#twins separated at birth & randomly distributed correlation (w/ independent env. effects)
error_twin <- rnorm(N, sd = env_sd)
traits_twin <- alpha + beta * snps + error_twin
cor(traits, traits_twin)

#also midparent correlation
cor(as.vector(traits_child), rep(midparent_traits, n_child_per_parents)) * sqrt(2)

#parent-offspring correlation
sapply(1:2, function(pi) sapply(1:n_child_per_parents, function(ci) 
  cor(traits[parent_pairs[,pi]], traits_child[,ci]) * 2))

#falconer's formula
#separated at birth? i.e. indep. env. effects?
traits_child_sep <- snps_child * beta + alpha + matrix(rnorm(N/2*n_child_per_parents , 0, sd = env_sd), ncol = 2)
(cor(traits, traits_twin) - cor(traits_child_sep)[1,2])*2
#same env. effects?
(1 - cor(traits_child)[1,2])*2

