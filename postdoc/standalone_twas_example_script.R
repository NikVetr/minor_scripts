#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

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

#### OK, now let's have it do multi-tissues, one gene, same loci affecting everything ####

## specify simulation parameters
n_y <- 10 #number of genes (or tissues)
bs_wanted_var <- 20 #variance of tissue-specific effect sizes
bs_wanted = round(rnorm(n = n_y,0,sqrt(bs_wanted_var)),2) #effect of each gene expression on z, the phenotype; if high, most of genotype signal should come from here
# bs_0s <- sample(1:n_y, size = round(n_y/2), replace = F) #simulate spike and slab
# bs_wanted[bs_0s] <- 0
h2_y <- rbeta(n_y,15,5) #heritability of gene expression
prop_resid_var_z <- rbeta(1,10,10) #residual variance in trait, conditional on gene expression and SNPs
weight_with_perfect_correlation <- 0.5
corr_mat_ys <- rlkj(n_y, 1) * (1-weight_with_perfect_correlation) + 
  weight_with_perfect_correlation #correlation structure between SNP effects on trait
n_obs = c(800, 1.5E3) #total number of individuals in samples 1 and 2
n_rep <- 50 #total number of loci influencing each gene
x_freq <- rbeta(n_rep,1,2) * 0.99 + 0.01 #allele frequency across loci


# simulate uncorrelated alleles
# x1 <- sapply(x_freq, function(xf) rbinom(n_obs[1], 2, prob = xf)) #genotypes in pop 1
# x2 <- sapply(x_freq, function(xf) rbinom(n_obs[2], 2, prob = xf)) #genotypes in comparable pop2

# simulate autocorrelated alleles by thresholding a MVN random variable 
# rbf_d <- function(dist, params = list(s = 1, l = 1)){params$s^2*exp(-(dist)^2/(2*params$l^2))}
ou_d <- function(dist, params = list(s = 1, l = 1)){params$s^2*exp(-abs(dist)/params$l)}
# brownian <- function(x1, x2, params = list(s=1)){return(params$s^2*min(x1, x2))}
# fastsim_MVN <- function(corr_step_multiplier = 0, n){
#   x <- cumsum(rnorm(n, sd = sqrt(1 / n))) + rev(cumsum(rnorm(n, sd = sqrt(1 / n)))) / sqrt(1+1/n)
#   x <- (x + rnorm(n, sd = sqrt(corr_step_multiplier / n))) / sqrt(1 + corr_step_multiplier / n)
#   x
# }
std_normal_thresholds <- qnorm(x_freq)
x_locations <- sort(runif(n_rep, 0, 2))
x_dists <- as.matrix(dist(x_locations))
x_cov <- ou_d(x_dists, params = list(s = 1, l = 2)) #covariance in underlying MVN 

#genotypes in population 1
all.same <- function(x) length(unique(x)) == 1
any_same = T
while(any_same){
  x1 <- replicate(2, rmvnorm(n = n_obs[1], mean = rep(0, n_rep), sigma = x_cov))
  x1 <- sapply(1:n_rep, function(li) (x1[,li,1] > std_normal_thresholds[li]) + (x1[,li,2] > std_normal_thresholds[li]))
  any_same = any(apply(x1, 2, all.same))
}
#genotypes in population 2
any_same = T
while(any_same){
  x2 <- replicate(2, rmvnorm(n = n_obs[2], mean = rep(0, n_rep), sigma = x_cov))
  x2 <- sapply(1:n_rep, function(li) (x2[,li,1] > std_normal_thresholds[li]) + (x2[,li,2] > std_normal_thresholds[li]))
  any_same = any(apply(x2, 2, all.same))
}
#check for invariant sites

apply(x2, 2, all.same)

#effect of each marginal variant at each locus on gene expression
bs_expr <- t(sapply(1:n_rep, function(li) rmvnorm(1, sigma = corr_mat_ys))) #correlated effects on expression

#expected values for gene expression in pops 1 and 2
y1_exp <- sapply(1:n_y, function(gi) x1 %*% t(t(bs_expr[,gi])))
y2_exp <- sapply(1:n_y, function(gi) x2 %*% t(t(bs_expr[,gi])))

#residual variance
sd_ys_obs <- sapply(1:n_y, function(gi) sd(c(y1_exp[,gi], y2_exp[,gi]))) #can work out exactly later
sd_ys_remaining <- sqrt(sd_ys_obs^2 * (1-h2_y) / h2_y)

#expression levels in genes for pop1
# y1 <- sapply(1:n_y, function(gi) y1_exp[,gi] + rnorm(n = n_obs[1], 0, sd_ys_remaining[gi]))
resids_ys_1 <- rmvnorm(n = n_obs[1], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y1 <- sapply(1:n_y, function(gi) y1_exp[,gi] + resids_ys_1[,gi]) #correlated residuals

#expression levels in genes for pop2
# y2 <- sapply(1:n_y, function(gi) y2_exp[,gi] + rnorm(n = n_obs[2], 0, sd_ys_remaining[gi]))
resids_ys_2 <- rmvnorm(n = n_obs[2], mean = rep(0, n_y), sigma = diag(sd_ys_remaining) %*% corr_mat_ys %*% diag(sd_ys_remaining))
y2 <- sapply(1:n_y, function(gi) y2_exp[,gi] + resids_ys_2[,gi]) #correlated residuals

#direct SNP effects on phenotype z for each population
bs_pheno_mean <- 2
bs_pheno <- rnorm(n_rep, mean = bs_pheno_mean, sd = 1)

z1_exp <- sapply(1:n_obs[1], function(indiv)
  sum(x1[indiv,] * bs_pheno) + #direct genetic effect 
    sum(y1[indiv,] * bs_wanted) #effect of gene expression
)

z2_exp <- sapply(1:n_obs[2], function(indiv)
  sum(x2[indiv,] * bs_pheno) + #direct genetic effect 
    sum(y2[indiv,] * bs_wanted) #effect of gene expression
)

#... can also work out exactly later
sd_z <- sqrt(var(c(z1_exp, z2_exp)) * prop_resid_var_z / (1-prop_resid_var_z))
z1 <- z1_exp + rnorm(n_obs[1], 0, sd_z)
z2 <- z2_exp + rnorm(n_obs[2], 0, sd_z)

#fit single regressions
eQTL_coef <- do.call(rbind, lapply(1:n_y, function(gi) sapply(1:n_rep, function(li) lm(c(y1[,gi]) ~ x1[,li])$coefficients[2])))
GWAS_coef <- sapply(1:n_rep, function(li) lm(z2 ~ x2[,li])$coefficients[2])

#coerce single regressions to multiple regressions
# Cov_x1 <- cov(x1) * (n_obs[1] - 1) / n_obs[1] + diag(n_rep) / n_obs[1]
Cov_x1 <- corpcor::cov.shrink(x1)
Cov_x2 <- corpcor::cov.shrink(x2)

Cov_x1_y <- eQTL_coef %*% diag(apply(x1, 2, var))
eQTL_coef <- t(pracma::pinv(Cov_x1) %*% t(Cov_x1_y))

Cov_x2_z <- GWAS_coef * apply(x2, 2, var)
GWAS_coef <- c(pracma::pinv(Cov_x2) %*% Cov_x2_z)

cov_x2_z_full <- cbind(rbind(Cov_x2, Cov_x2_z), c(Cov_x2_z, var(z2)))
r2_x2_z <- t(cov2cor(cov_x2_z_full)[n_rep+1,-(n_rep+1)]) %*% pracma::pinv(cov2cor(Cov_x2)) %*% t(t(cov2cor(cov_x2_z_full)[n_rep+1,-(n_rep+1)]))
adj_r2_x2_z <- 1-((1- r2_x2_z) * (n_obs[2]-1) / (n_obs[2]-n_rep-1))
cov_GWAS_coef <- pracma::pinv(Cov_x2) / c(((n_obs[2]-1) / (var(z2) * (1-adj_r2_x2_z))))
GWAS_coef_SE <- sqrt(diag(cov_GWAS_coef))

cov_x1_y_full <- lapply(1:n_y, function(gi) cbind(rbind(Cov_x1, Cov_x1_y[gi,]), c(Cov_x1_y[gi,], var(y1[,gi]))))
r2_x1_y <- sapply(1:n_y, function(gi) t(cov2cor(cov_x1_y_full[[gi]])[n_rep+1,-(n_rep+1)]) %*% pracma::pinv(cov2cor(Cov_x1)) %*% t(t(cov2cor(cov_x1_y_full[[gi]])[n_rep+1,-(n_rep+1)])))
adj_r2_x1_y <- 1-((1- r2_x1_y) * (n_obs[1]-1) / (n_obs[1]-n_rep-1))
cov_eQTL_coef <- lapply(1:n_y, function(gi) pracma::pinv(Cov_x1) / c(((n_obs[1]-1) / (var(y1[,gi]) * (1-adj_r2_x1_y[gi])))))
eQTL_coef_SE <- t(sapply(1:n_y, function(gi) sqrt(diag(cov_eQTL_coef[[gi]]))))

#alternatively, fit multiple regressions, and pretend we got their coefficients from the single regressions
eQTL_coef <- do.call(rbind, lapply(1:n_y, function(gi) lm(y1[,gi] ~ x1)$coefficients[-1])) #fitting the OLS eQTL model
cov_eQTL_coef <- lapply(1:n_y, function(gi) vcov(lm(y1[,gi] ~ x1)))
eQTL_coef_SE <- do.call(rbind, lapply(1:n_y, function(gi) summary(lm(y1[,gi] ~ x1))$coefficients[-1,2])) #not very efficient but w/e
GWAS_coef <- lm(z2 ~ x2)$coefficients[-1] #fitting the OLS GWAS model
cov_GWAS_coef <- vcov(lm(z2 ~ x2))[-1,-1]
GWAS_coef_SE <- summary(lm(z2 ~ x2))$coefficients[-1,2]


#fit coefficients model
# standardize <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
par(mfrow = c(4,2), mar = c(4,7,3,3), xpd = NA)
fit3 <- rlm(t(t(GWAS_coef)) ~ t(eQTL_coef), na.action = na.exclude)
# plot(standardize(resid(fit3)), standardize(unlist(bs_pheno)))
# fit3 <- lm(t(t(GWAS_coef)) ~ t(eQTL_coef), na.action = na.exclude)
fit3_summary <- summary(fit3)
plot(bs_wanted, y = fit3$coefficients[-1], cex.lab = 1.25, pch = 19, col = adjustcolor(1, 0.5), cex = 1.5, main = "Robust IWLS ('M-Estimator')",
     ylim = range(c(fit3$coefficients[-1] + 2*fit3_summary$coefficients[,"Std. Error"][-1], 
                    fit3$coefficients[-1] - 2*fit3_summary$coefficients[,"Std. Error"][-1])),
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i} + \\epsilon$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = fit3$coefficients[-1][i] + 2*fit3_summary$coefficients[,"Std. Error"][-1][i], 
           x1 = bs_wanted[i], y1 = fit3$coefficients[-1][i] - 2*fit3_summary$coefficients[,"Std. Error"][-1][i])}
abline(0,1, col = 2, lty = 2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.03*(par("usr")[4] - par("usr")[3]), labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, fit3$coefficients[-1]), 3))), pos = 2)

#95% CI
in95CI <- abs(fit3$coefficients[-1] - bs_wanted) < qt(0.975, df = n_y-1)*fit3_summary$coefficients[,"Std. Error"][-1]
text(x = par("usr")[2], y = par("usr")[3] + 0.12*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("Prop. in CI_{95} = ", round(mean(in95CI), 3))), pos = 2)

#additional coverage check
tscore <- (fit3$coefficients[-1] - bs_wanted) / fit3_summary$coefficients[,"Std. Error"][-1]
hist(pt(tscore, df = n_y-1), breaks = seq(0,1,length.out = 20), xlab = "Quantile of True Values in T-Distribution")

#plot raw coefficient estimates
plot(bs_expr, t(eQTL_coef), pch = 19, col = adjustcolor(1, 0.5), 
     ylim = range(c(t(eQTL_coef) + 2*t(eQTL_coef_SE), t(eQTL_coef) - 2*t(eQTL_coef_SE))),
     xlab = "true eQTL effect sizes", ylab = "estimated (OLS) eQTL effect sizes"); 
for(i in 1:length(eQTL_coef_SE)){
  segments(x0 = bs_expr[i], y0 = t(eQTL_coef)[i] + 2*t(eQTL_coef_SE)[i], 
           x1 = bs_expr[i], y1 = t(eQTL_coef)[i] - 2*t(eQTL_coef_SE)[i],
           col = adjustcolor(1, 0.5))}
abline(0,1,col=2,lty=2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))

plot(bs_pheno, GWAS_coef, pch = 19, col = adjustcolor(1, 0.5), 
     ylim = range(GWAS_coef + 2*GWAS_coef_SE, GWAS_coef - 2*GWAS_coef_SE),
     xlab = "true horizontal pleiotropic (i.e. direct) effect sizes", ylab = "estimated (OLS) total SNP effect sizes"); 
for(i in 1:length(bs_pheno)){
  segments(x0 = bs_pheno[i], y0 = (GWAS_coef)[i] + 2*(GWAS_coef_SE)[i], 
           x1 = bs_pheno[i], y1 = (GWAS_coef)[i] - 2*(GWAS_coef_SE)[i],
           col = adjustcolor(1, 0.5))}
abline(0,1,col=2,lty=2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SE"), lty = c(2,1), col = c(2,1))

bs_pheno_mean
fit3$coefficients[1]

# try out the Bayesian model #
d <- list(eQTL_coef = t(eQTL_coef),
          eQTL_coef_SE = t(eQTL_coef_SE),
          GWAS_coef = c(GWAS_coef),
          GWAS_coef_SE = c(GWAS_coef_SE),
          n_y = n_y,
          n_rep = n_rep)

stan_program <- "
data {
  int<lower=0> n_y;
  int<lower=0> n_rep;
  real GWAS_coef[n_rep];
  real GWAS_coef_SE[n_rep];
  matrix[n_rep, n_y] eQTL_coef;
  matrix[n_rep, n_y] eQTL_coef_SE;
}
parameters {
  vector[n_y] beta;
  real<lower=0> beta_var;
  real<lower=0> sigma2;
  real alpha;
  matrix[n_rep, n_y] eQTL_coef_AUG;
  real GWAS_coef_AUG[n_rep];
}
model {
  // intermediate variables
  vector[n_y * n_rep] mu;
  real sigma;

  // priors
  alpha ~ normal(0,5);
  sigma2 ~ exponential(0.1);
  beta_var ~ exponential(0.1);
  beta ~ normal(0, sqrt(beta_var));
  
  // upstream inferential uncertainty, since normals are symmetric
  // normal() doesn't seem to accept matrix arguments but does accept vector args
  for(ri in 1:n_rep){
    eQTL_coef[ri,] ~ normal(eQTL_coef_AUG[ri,], eQTL_coef_SE[ri,]);  
  }
  GWAS_coef ~ normal(GWAS_coef_AUG, GWAS_coef_SE);

  // actual model and likelihood
  mu = alpha + eQTL_coef_AUG * beta;
  sigma = sqrt(sigma2);
  // GWAS_coef_AUG ~ student_t(5, mu, sigma);
  GWAS_coef_AUG ~ normal(mu, sigma);
  
}
"

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)  
}
mod <- cmdstan_model(f)

#fit model
out <- mod$sample(chains = 4, iter_sampling = 4E3, iter_warmup = 4E3, data = d, parallel_chains = 4, adapt_delta = 0.95)
samps <- data.frame(as_draws_df(out$draws()))

#print out some hyperparameter estimates
bs_pheno_mean
mean(samps[,grep("alpha", colnames(samps))])

bs_wanted_var
mean(samps[,grep("beta_var", colnames(samps))])

#visualize posterior means and sds
plot(bs_wanted, y = apply(samps[,grep("beta\\.", colnames(samps))],2,mean), cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Hierarchical Bayesian Analysis",
     ylim = range(c(apply(samps[,grep("beta\\.", colnames(samps))],2,mean) + 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd),
                    apply(samps[,grep("beta\\.", colnames(samps))],2,mean) - 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))),
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i} + \\epsilon$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = (apply(samps[,grep("beta\\.", colnames(samps))],2,mean) + 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))[i],
           x1 = bs_wanted[i], y1 = (apply(samps[,grep("beta\\.", colnames(samps))],2,mean) - 2*apply(samps[,grep("beta\\.", colnames(samps))],2,sd))[i])}
abline(0,1, col = 2, lty = 2, lwd = 2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SD"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.02*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, apply(samps[,grep("beta\\.", colnames(samps))],2,mean)), 3))), pos = 2)
#check 95% coverage
in95CI_Bayes <- apply(samps[,grep("beta\\.", colnames(samps))],2,quantile, probs = c(0.025,0.975))
in95CI_Bayes <- bs_wanted > in95CI_Bayes[1,] & bs_wanted < in95CI_Bayes[2,]  
text(x = par("usr")[2], y = par("usr")[3] + 0.12*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("Prop. in CI_{95} = ", round(mean(in95CI_Bayes), 3))), pos = 2)


#compare bayesian vs OLS estimates
plot(fit3$coefficients[-1], apply(samps[,grep("beta\\.", colnames(samps))],2,mean),
     ylab = "Posterior Means", xlab = "Robust IWLS Estimates", cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Shrinkage Comparison") 
abline(0,1, lty = 2, col = 2, xpd = F)
abline(lm(apply(samps[,grep("beta\\.", colnames(samps))],2,mean) ~ fit3$coefficients[-1]), col = 4, lty = 3, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "OLS fit to Points"), lty = c(2,1), col = c(2,4))


# try out the multivariate Bayesian model #

#buffer covariance matrices so they fit in an array?
cov_eQTL_coef_array <- aperm(abind::abind(cov_eQTL_coef, along = 3), c(3,1,2))
cholcov_eQTL_coef_array <- aperm(abind::abind(lapply(cov_eQTL_coef, function(x) t(chol(x))), along = 3), c(3,1,2))

# d_mv <- list(eQTL_coef = t(eQTL_coef),
#              cov_eQTL_coef_array = cov_eQTL_coef_array,
#              GWAS_coef = c(GWAS_coef),
#              cov_GWAS_coef = cov_GWAS_coef,
#              n_y = n_y,
#              n_rep = n_rep)

d_mv <- list(eQTL_coef = t(eQTL_coef),
             cholcov_eQTL_coef_array = cholcov_eQTL_coef_array,
             GWAS_coef = c(GWAS_coef),
             cholcov_GWAS_coef = t(chol(cov_GWAS_coef)),
             n_y = n_y,
             n_rep = n_rep)

stan_program_mv <- "
data {
  int<lower=0> n_y;
  int<lower=0> n_rep;
  vector[n_rep] GWAS_coef;
  matrix[n_rep, n_rep] cholcov_GWAS_coef;
  matrix[n_rep, n_y] eQTL_coef;
  matrix[n_rep, n_rep] cholcov_eQTL_coef_array[n_y];
}
parameters {
  vector[n_y] beta;
  real<lower=0> beta_var;
  real<lower=0> sigma2;
  real alpha;
  matrix[n_rep, n_y] eQTL_coef_AUG;
  vector[n_rep] GWAS_coef_AUG;
}
model {
  // intermediate variables
  vector[n_y * n_rep] mu;
  real sigma;

  // priors
  alpha ~ normal(0,5);
  sigma2 ~ exponential(0.1);
  beta_var ~ exponential(0.1);
  beta ~ normal(0, sqrt(beta_var));
  
  // upstream inferential uncertainty, since normals are symmetric
  // normal() doesn't seem to accept matrix arguments but does accept vector args
  for(gi in 1:n_y){
    // eQTL_coef[,gi] ~ multi_normal(eQTL_coef_AUG[,gi], cov_eQTL_coef_array[gi]);  
    eQTL_coef[,gi] ~ multi_normal_cholesky(eQTL_coef_AUG[,gi], cholcov_eQTL_coef_array[gi]);  
  }
  // GWAS_coef ~ multi_normal(GWAS_coef_AUG, cov_GWAS_coef);
  GWAS_coef ~ multi_normal_cholesky(GWAS_coef_AUG, cholcov_GWAS_coef);

  // actual model and likelihood
  mu = alpha + eQTL_coef_AUG * beta;
  sigma = sqrt(sigma2);
  // GWAS_coef_AUG ~ student_t(5, mu, sigma);
  GWAS_coef_AUG ~ normal(mu, sigma);
  
}
"

if(!exists("curr_stan_program_mv") || stan_program_mv != curr_stan_program_mv){
  curr_stan_program_mv <- stan_program_mv
  f_mv <- write_stan_file(stan_program_mv)  
}
mod_mv <- cmdstan_model(f_mv)

#fit model
out_mv <- mod_mv$sample(chains = 4, iter_sampling = 4E3, iter_warmup = 4E3, data = d_mv, parallel_chains = 4, adapt_delta = 0.95)
samps_mv <- data.frame(as_draws_df(out_mv$draws()))

#print out some hyperparameter estimates
bs_pheno_mean
mean(samps_mv[,grep("alpha", colnames(samps_mv))])

bs_wanted_var
mean(samps_mv[,grep("beta_var", colnames(samps_mv))])

#visualize posterior means and sds
plot(bs_wanted, y = apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean), cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Mv. Error Hierarchical Bayesian Analysis",
     ylim = range(c(apply(samps[,grep("beta\\.", colnames(samps_mv))],2,mean) + 2*apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,sd),
                    apply(samps[,grep("beta\\.", colnames(samps_mv))],2,mean) - 2*apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,sd))),
     xlab = latex2exp::TeX("true value $\\beta_i$"), ylab = latex2exp::TeX("estimated value $\\hat{\\beta_i}$ from $\\Beta_z ~ \\Sigma \\beta_i\\Beta_{y_i} + \\epsilon$"))
for(i in 1:length(fit3$coefficients[-1])){
  segments(x0 = bs_wanted[i], y0 = (apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean) + 2*apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,sd))[i],
           x1 = bs_wanted[i], y1 = (apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean) - 2*apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,sd))[i])}
abline(0,1, col = 2, lty = 2, lwd = 2, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "±2SD"), lty = c(2,1), col = c(2,1))
text(x = par("usr")[2], y = par("usr")[3] + 0.02*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("r_{ij} = ", round(cor(bs_wanted, apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean)), 3))), pos = 2)
#check 95% coverage
in95CI_Bayes_mv <- apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,quantile, probs = c(0.025,0.975))
in95CI_Bayes_mv <- bs_wanted > in95CI_Bayes_mv[1,] & bs_wanted < in95CI_Bayes_mv[2,]  
text(x = par("usr")[2], y = par("usr")[3] + 0.12*(par("usr")[4] - par("usr")[3]), 
     labels = latex2exp::TeX(paste0("Prop. in CI_{95} = ", round(mean(in95CI_Bayes_mv), 3))), pos = 2)


#compare bayesian vs OLS estimates
plot(apply(samps[,grep("beta\\.", colnames(samps))],2,mean), 
     apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean),
     ylab = "Posterior Means MV. Errors", xlab = "Posterior Means UV. Errors", cex.lab = 1.25, pch = 19, 
     col = adjustcolor(1, 0.5), cex = 1.5, main = "Shrinkage Comparison") 
abline(0,1, lty = 2, col = 2, xpd = F)
abline(lm(apply(samps_mv[,grep("beta\\.", colnames(samps_mv))],2,mean) ~ apply(samps[,grep("beta\\.", colnames(samps))],2,mean)), 
       col = 4, lty = 3, xpd = F)
legend(lwd = 1, x = "topleft", legend = c("1-to-1 line", "OLS fit to Points"), lty = c(2,1), col = c(2,4))


