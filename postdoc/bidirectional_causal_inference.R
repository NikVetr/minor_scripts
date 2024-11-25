library(cmdstanr)
library(posterior)
#### sim data and fit freq mod ####

#specify high level parameters
standardize_vars <- F
p = 50
n = 200
prop_yvar = 0.5 #proportion non-aleatoric variance (ie, not in error term)
prop_zvar = 0.4
sd_b_xz_scale <- 0.2
sd_b_xz = sd_b_xz_scale / sqrt(p)
b_yz = -4
sd_u_scale <- 0 #relative confounding effect of unobserved_variable
sd_u <- sd_u_scale * sqrt(p)

#set stage to simulate correlated counts in {0,1,2}
rop <- 0.8
rs <- c(1, rop^(1:p / 3))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)
cholR <- chol(R)
p_freqs <- rbeta(p, 2, 2)
p_liab <- qnorm(p_freqs)

#actually simulate these
x1 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs <- t(t(liabs) + p_liab)
  probs <- pnorm(liabs)
  apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
}))

x2 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs <- t(t(liabs) + p_liab)
  probs <- pnorm(liabs)
  apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
}))

x3 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs <- t(t(liabs) + p_liab)
  probs <- pnorm(liabs)
  apply(probs, c(1,2), function(prob) rbinom(n = 1, prob = prob, size = 1))  
}))

#further randomize effects of x?
# x1 <- x1 %*% matrix(rnorm(p^2), p, p)
# x2 <- x2 %*% matrix(rnorm(p^2), p, p)

#specify mediator, outcome, and exposure relationships
b_xy <- rnorm(p)
b_xz <- rnorm(p) * sd_b_xz #when sd_b_xz = 0, no horizontal pleiotropy 

### simulate data for ###
### z <- x -> y -> z  ###

#simulate y for pop 1
mu_y1 <- c(x1 %*% t(t(b_xy)))
u1 <- rnorm(n) * sd_u
e_y1 <- rnorm(n) * c(sqrt(var(mu_y1 + u1) / prop_yvar * (1-prop_yvar)))
y1 <- mu_y1 + u1 + e_y1

#simulate y and z for pop 2
mu_y2 <- c(x2 %*% t(t(b_xy)))
u2 <- rnorm(n) * sd_u
e_y2 <- rnorm(n) * c(sqrt(var(mu_y2 + u2) / prop_yvar * (1-prop_yvar)))
y2 <- mu_y2 + u2 + e_y2
mu_z2 <- y2 * b_yz + c(x2 %*% t(t(b_xz)))
e_z2 <- rnorm(n) * c(sqrt(var(mu_z2 + u2) / prop_zvar * (1-prop_zvar)))
z2 <- mu_z2 + u2 + e_z2

#standardize variables?
if(standardize_vars){
  x1 <- matrix(c(scale(x1)), n, p)
  x2 <- matrix(c(scale(x2)), n, p)
  y1 <- c(scale(y1))
  z2 <- c(scale(z2))
  # y <- y - mean(y)
  # z <- z - mean(z)
}

#perform univariate regressions
fits_xy <- do.call(rbind, lapply(1:p, function(i) summary(lm(y1 ~ 1 + x1[,i]))$coefficients["x1[, i]",1:2]))
fits_xz <- do.call(rbind, lapply(1:p, function(i) summary(lm(z2 ~ 1 + x2[,i]))$coefficients["x2[, i]",1:2]))

#### recover multiple regression ###

covx1 <- cov(x1)
covx2 <- cov(x2)
covx1 <- covx2 <- (covx1 + covx2)/2

#first for y
covx1i <- solve(covx1)
varx1 <- diag(covx1)
covxy1 <- fits_xy[,"Estimate"] * varx1
mfit_xy <- covx1i %*% fits_xy[,"Estimate"] * varx1
r2xy <- c(t(fits_xy[,"Estimate"]) %*% covx1i %*% fits_xy[,"Estimate"])
r2adj_xy <- 1 - ((1 - r2xy) * (n - 1) / (n - p - 1))
covxy1_multi <- solve(covx1) / c(((n-1) / (1 * (1-r2adj_xy))))
fit_xy_coerced <- cbind(coef = mfit_xy, SE = sqrt(diag(covxy1_multi)))

#then for z
covx2i <- solve(covx2)
varx2 <- diag(covx2)
covxz2 <- fits_xz[,"Estimate"] * varx2
mfit_xz <- covx2i %*% fits_xz[,"Estimate"] * varx2
r2xz <- c(t(fits_xz[,"Estimate"]) %*% covx2i %*% fits_xz[,"Estimate"])
r2adj_xz <- 1 - ((1 - r2xz) * (n - 1) / (n - p - 1))
covxz2_multi <- solve(covx2) / c(((n-1) / (1 * (1-r2adj_xz))))
fit_xz_coerced <- cbind(coef = mfit_xz, SE = sqrt(diag(covxz2_multi)))

#compare to straightforward multiple regression
fit_xy <- summary(lm(y1 ~ 1 + x1))$coefficients[paste0("x1", 1:p), 1:2]
fit_xz <- summary(lm(z2 ~ 1 + x2))$coefficients[paste0("x2", 1:p), 1:2]
fit_xy - fit_xy_coerced
fit_xz - fit_xz_coerced

#### recover PCR ###
# eig_x1 <- eigen(covx1)
# eig_x2 <- eigen(covx2)
#can't use independent eigendecompositions for x2 and x2 transformations
#eigenvector sign gets flipped! which yields expected slope of 0
#need to use single eigendecomposition
eig_x2 <- eig_x1 <<- eigen((covx1 + covx2) / 2)

x1_eval <- eig_x1$values
x1_evec <- eig_x1$vectors

x2_eval <- eig_x2$values
x2_evec <- eig_x2$vectors

#coefficients
mfit_xy_pcr <- c(t(mfit_xy) %*% x1_evec) * sqrt(x1_eval)
mfit_xz_pcr <- c(t(mfit_xz) %*% x2_evec) * sqrt(x2_eval)

#can also try to do something like a cholesky transformation?
chol_x <- chol((covx1 + covx2) / 2)
chol_xi <- solve(chol_x)
x1_ind <- x1 %*% solve(chol_xi)
x2_ind <- x2 %*% solve(chol_xi)
fit_xy_chol <- data.frame(summary(lm(y1 ~ 1 + x1_ind))$coefficients[paste0("x1_ind", 1:p), 1:2])
fit_xz_chol <- data.frame(summary(lm(z2 ~ 1 + x2_ind))$coefficients[paste0("x2_ind", 1:p), 1:2])

#vcov matrices
covxy1_multi_pcr <- t(x1_evec %*% diag(sqrt(x1_eval))) %*% covxy1_multi %*% (x1_evec %*% diag(sqrt(x1_eval)))
covxz2_multi_pcr <- t(x2_evec %*% diag(sqrt(x2_eval))) %*% covxz2_multi %*% (x2_evec %*% diag(sqrt(x2_eval)))

#put together
fit_xy_pcr_coerced <- cbind(coef = mfit_xy_pcr, SE = sqrt(diag(covxy1_multi_pcr)))
fit_xz_pcr_coerced <- cbind(coef = mfit_xz_pcr, SE = sqrt(diag(covxz2_multi_pcr)))

#compare to straightforward PCR
x1_PCs <- x1 %*% x1_evec %*% diag(1/sqrt(x1_eval))
x2_PCs <- x2 %*% x2_evec %*% diag(1/sqrt(x2_eval))
fit_xy_pcr <- summary(lm(y1 ~ 1 + x1_PCs))$coefficients[paste0("x1_PCs", 1:p), 1:2]
fit_xz_pcr <- summary(lm(z2 ~ 1 + x2_PCs))$coefficients[paste0("x2_PCs", 1:p), 1:2]
# head(sort(abs(cov2cor(vcov(lm(y1 ~ 1 + x1_PCs))[-1,-1])[upper.tri(diag(p))]), T))
# head(sort(abs(cov2cor(vcov(lm(z2 ~ 1 + x2_PCs))[-1,-1])[upper.tri(diag(p))]), T))
fit_xy_pcr - fit_xy_pcr_coerced
fit_xz_pcr - fit_xz_pcr_coerced

#### run simple causal regression ####
par(mfrow = c(1, 3), mar = c(5,5,2,2))
plot(mfit_xy, mfit_xz, col = NA, 
     xlab = "Multiple regression coefficients for X → Y", 
     ylab = "Multiple regression coefficients for X → Z")
text(mfit_xy, mfit_xz, labels = 1:p)
text(x = par("usr")[1], y = par("usr")[4] - diff(par("usr")[3:4])/30, 
     label = "# = Instrument Index", xpd = NA, pos = 4, col = 2)

plot(mfit_xy_pcr, mfit_xz_pcr, col = NA, 
     xlab = "PCR coefficients for X → Y", 
     ylab = "PCR coefficients for X → Z")
text(mfit_xy_pcr, mfit_xz_pcr, labels = 1:p, 
     col = viridis::magma(p, end = 0.8))
text(x = par("usr")[1], y = par("usr")[4] - diff(par("usr")[3:4])/30, 
     label = "# = PC Axis Index", xpd = NA, pos = 4, col = 2)

nPCs_fits <- data.frame(do.call(rbind, lapply(3:p, function(i) summary(lm(mfit_xz_pcr[1:i] ~ 1 + mfit_xy_pcr[1:i]))$coefficients[2,1:4])))
# nPCs_fits <- data.frame(do.call(rbind, lapply(3:p, function(i) summary(lm(mfit_xz_pcr[1:i] ~ 0 + mfit_xy_pcr[1:i]))$coefficients[1,1:4])))
# plot(-log10(nPCs_fits$Pr...t..), type = "l")
CI95 <- data.frame(lb = nPCs_fits$Estimate - nPCs_fits$Std..Error * qnorm(0.975),
                   ub = nPCs_fits$Estimate + nPCs_fits$Std..Error * qnorm(0.975)
)
plot(3:p, nPCs_fits$Estimate, type = "l", ylim = range(c(CI95, 0, b_yz)), xlab = "number of PCs", ylab = "Estimate b_yz")
polygon(x = c(3:p, p:3, 3), y = c(CI95$lb, rev(CI95$ub), CI95$lb[1]), col = adjustcolor(1,0.25), border = NA)
abline(h = b_yz, col = 2, lty = 2)
abline(h = 0, lty = 3)
legend(x = 15, y = -15, 
       legend = c("Estimate b_yz", "95% CI (shaded)", "True b_yz", "Zero line"), 
       col = c(1, adjustcolor(1, 0.25), 2, 1), 
       lty = c(1, NA, 2, 3), 
       pch = c(NA, 15, NA, NA),
       fill = c(NA, adjustcolor(1, 0), NA, NA), 
       border = c(NA, NA, NA, NA), 
       bty = "n", 
       cex = 0.8, 
       pt.cex = 3, 
       y.intersp = 1.5)

#####

summary(lm(mfit_xz ~ 1 + mfit_xy))$coefficients[2,1:2]
summary(lm(mfit_xz_pcr ~ 1 + mfit_xy_pcr))$coefficients[2,1:2]
summary(lm(fit_xz_chol$Estimate ~ 1 + fit_xy_chol$Estimate))$coefficients[2,1:2]


#eventually, PC coefs get to be random noise?
#as you add random noise, slope gets increasingly confident at 0?
#but not if you propagate uncertainty 
#(the true coefs get centered on 0, but that is fine! enrichment relationship still holds)


#### check true PCR coefs? ####
b_xy_PCs <- solve(t(x1_PCs) %*% x1_PCs) %*% t(x1_PCs) %*% mu_y1
b_xz_PCs <- (solve(t(cbind(y2, x2_PCs)) %*% cbind(y2, x2_PCs)) %*% t(cbind(y2, x2_PCs)) %*% mu_z2)[-1]
b_yz_PCs <- (solve(t(cbind(y2, x2_PCs)) %*% cbind(y2, x2_PCs)) %*% t(cbind(y2, x2_PCs)) %*% mu_z2)[1]
b_xz_noy_pcs <- solve(t(x2_PCs) %*% x2_PCs) %*% t(x2_PCs) %*% mu_z2
solve(t(x2_PCs) %*% x2_PCs) %*% t(x2_PCs) %*% mu_z2
plot(1:p, b_xy_PCs)
plot(b_xz_noy_pcs, b_xy_PCs)
summary(lm(b_xz_noy_pcs ~ b_xy_PCs))$coefficients[2,]

plot(summary(lm(mu_z2 ~ x2))$coefficients[-1,1], summary(lm(mu_y1 ~ x1))$coefficients[-1,1])
plot(summary(lm(mu_z2 ~ x2_PCs))$coefficients[-1,1], summary(lm(mu_y1 ~ x1_PCs))$coefficients[-1,1], col = NA)
text(summary(lm(mu_z2 ~ x2_PCs))$coefficients[-1,1], summary(lm(mu_y1 ~ x1_PCs))$coefficients[-1,1], labels = 1:p)
hist(summary(lm(mu_z2 ~ x2_PCs))$coefficients[-1,1], breaks = 100)

#### bidirectional inference conditioning! ####
fishers_method <- function(pvalues) {
  test_statistic <- -2 * sum(log(pvalues))
  degrees_of_freedom <- 2 * length(pvalues)
  combined_pvalue <- pchisq(test_statistic, degrees_of_freedom, lower.tail = FALSE)
  return(combined_pvalue)
}

#first in the same population
par(mfrow = c(2,1))
yz_pvals <- summary(lm(z2 ~ x2 + y2))$coefficients[2:(p+1),4]
zy_pvals <- summary(lm(y2 ~ x2 + z2))$coefficients[2:(p+1),4]
hist(yz_pvals)
hist(zy_pvals)

#then in a new population where we predict stuff
pred_y3 <- predict(lm(y1 ~ 1 + x1), newdata = data.frame(x3))
pred_z3 <- predict(lm(z2 ~ 1 + x2), newdata = data.frame(x3))
pred_yz_pvals <- summary(lm(pred_y3 ~ x3 + pred_z3))$coefficients[2:(p+1),4]
pred_zy_pvals <- summary(lm(pred_z3 ~ x3 + pred_y3))$coefficients[2:(p+1),4]
# hist(pred_yz_pvals)
# hist(pred_zy_pvals)

#check pvals
fishers_method(yz_pvals)
fishers_method(zy_pvals)
fishers_method(pred_yz_pvals)
fishers_method(pred_zy_pvals)


#what if we just predict one of these?
pred_y2 <- predict(lm(y1 ~ 1 + x1), newdata = data.frame(x2))
pred_z1 <- predict(lm(z2 ~ 1 + x2), newdata = data.frame(x1))

yz_pz_pvals <- summary(lm(pred_z1 ~ x1 + y1))$coefficients[2:(p+1),4]
zy_pz_pvals <- summary(lm(y1 ~ x1 + pred_z1))$coefficients[2:(p+1),4]
yz_py_pvals <- summary(lm(pred_y2 ~ x2 + z2))$coefficients[2:(p+1),4]
zy_py_pvals <- summary(lm(z2 ~ x2 + pred_y2))$coefficients[2:(p+1),4]

fishers_method(yz_pz_pvals)
fishers_method(zy_pz_pvals)
fishers_method(yz_py_pvals)
fishers_method(zy_py_pvals)


#### run Bayes ####

#specify uncertainty propagation model
dat_1 <- list(p = p,
              b_xy_est = fit_xy_pcr[,"Estimate"],
              b_xz_est = fit_xz_pcr[,"Estimate"],
              b_xy_se = fit_xy_pcr[,"Std. Error"],
              b_xz_se = fit_xz_pcr[,"Std. Error"])

dat_2 <- list(p = p,
              b_xy_est = fit_xz_pcr[,"Estimate"],
              b_xz_est = fit_xy_pcr[,"Estimate"],
              b_xy_se = fit_xz_pcr[,"Std. Error"],
              b_xz_se = fit_xy_pcr[,"Std. Error"])

stan_model_1 <- 
"data{
  int p;
  vector[p] b_xy_est;
  vector<lower=0>[p] b_xy_se;
  vector[p] b_xz_est;
  vector<lower=0>[p] b_xz_se;
}

parameters{
  //true coefficients
  vector[p] b_xy_raw;
  vector[p] b_xz_raw;
  real b_yz;
  
  //intercepts and errors
  real a;
  real<lower=0> b_xy_sd;
  real<lower=0> b_xz_sd;
}

transformed parameters{
  //re-center true coefs
  vector[p] b_xy = b_xy_raw * b_xy_sd;
  vector[p] b_xz = b_xz_raw * b_xz_sd + a + b_yz * b_xy;
}

model{
  //observed coef model
  a ~ std_normal();
  b_yz ~ normal(0,5);
  b_xy_est ~ normal(b_xy, b_xy_se);
  b_xz_est ~ normal(b_xz, b_xz_se);
  
  //true coef model
  b_xy_raw ~ std_normal();
  b_xz_raw ~ std_normal();
}
"
#fit models
mod_1 <- cmdstan_model(write_stan_file(stan_model_1))
fit_1 <- mod_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.95, parallel_chains = 4,
                      refresh = 100, max_treedepth = 10, 
                      thin = 1, init = 0.1, data = dat_1)
fit_2 <- mod_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.99, parallel_chains = 4,
                      refresh = 100, max_treedepth = 10, 
                      thin = 1, init = 0.1, data = dat_2)

#check mcmc diagnostics
summ_1 <- fit_1$summary()
print(summ_1[order(summ_1$rhat, decreasing = T),])

summ_2 <- fit_2$summary()
print(summ_2[order(summ_2$rhat, decreasing = T),])

#extract samples
samps_1 <- data.frame(as_draws_df(fit_1$draws()))
samps_2 <- data.frame(as_draws_df(fit_2$draws()))

#inspect focal effect estimate
plot(samps_1$b_yz, type = "l"); abline(h=b_yz, col=2, lty=2, lwd=4)
hist(samps_1$b_yz, breaks = 100); abline(v=b_yz, col=2, lty=2, lwd=4)

plot(samps_2$b_yz, type = "l")
hist(samps_2$b_yz, breaks = 100); abline(v=b_yz, col=2, lty=2, lwd=4)

hist(samps_1$lp__)
hist(samps_2$lp__)
mean((samps_1$lp__-samps_2$lp__) > 0)

hist(samps_1$b_xy_sd)
hist(samps_1$b_xz_sd)
hist(samps_2$b_xy_sd)
hist(samps_2$b_xz_sd)

#inspect associations statistics
b_xy_samps_1 <- samps_1[,paste0("b_xy.", 1:p, ".")]
b_xz_samps_1 <- samps_1[,paste0("b_xz.", 1:p, ".")]
b_xy_samps_2 <- samps_2[,paste0("b_xy.", 1:p, ".")]
b_xz_samps_2 <- samps_2[,paste0("b_xz.", 1:p, ".")]
par(mfrow = c(2,1))
plot(apply(b_xy_samps_1, 2, mean), apply(b_xz_samps_1, 2, mean)); abline(0,b_yz,col=2)
plot(apply(b_xz_samps_2, 2, mean), apply(b_xy_samps_2, 2, mean)); abline(0,b_yz,col=2)

plot(mfit_xy_pcr, mfit_xz_pcr, col = NA)
text(mfit_xy_pcr, mfit_xz_pcr, labels = 1:p, 
     col = viridis::magma(p, end = 0.8))

#try bridge sampling? 
#specify a mixture model over the two alternatives?
#and integrate the mixture parameter out? 

stan_model_2 <- "
data {
  int<lower=0> p;
  vector[p] b_xy_est;
  vector<lower=0>[p] b_xy_se;
  vector[p] b_xz_est;
  vector<lower=0>[p] b_xz_se;
}

parameters {
  // True coefficients
  vector[p] b_xy_raw;
  vector[p] b_xz_raw;
  real b_yz;
  
  // Intercepts and errors
  real a;
  real<lower=0> b_xy_sd;
  real<lower=0> b_xz_sd;
  
  // Mixture proportion
  real<lower=0, upper=1> theta;
}

transformed parameters {
  vector[p] b_xy_1 = b_xy_raw * b_xy_sd;
  vector[p] b_xz_1 = b_xz_raw * b_xz_sd + a + b_yz * b_xy_1;
  
  vector[p] b_xy_2 = b_xz_raw * b_xz_sd;
  vector[p] b_xz_2 = b_xy_raw * b_xy_sd + a + b_yz * b_xy_2;
}

model {
  real log_lik_1;
  real log_lik_2;
  
  // Priors
  theta ~ beta(1, 1);  // Uniform prior for the mixture proportion
  a ~ std_normal();
  b_yz ~ normal(0, 5);
  b_xy_raw ~ std_normal();
  b_xz_raw ~ std_normal();
  b_xy_sd ~ std_normal();
  b_xz_sd ~ std_normal();
  
  // Likelihood for Alternative 1
  log_lik_1 = normal_lpdf(b_xy_est | b_xy_1, b_xy_se) + 
              normal_lpdf(b_xz_est | b_xz_1, b_xz_se);
              
  // Likelihood for Alternative 2 (flipped)
  log_lik_2 = normal_lpdf(b_xy_est | b_xy_2, b_xy_se) + 
              normal_lpdf(b_xz_est | b_xz_2, b_xz_se);
              
  // Combine the log-likelihoods using log_sum_exp
  target += log_sum_exp(log(theta) + log_lik_1, log1m(theta) + log_lik_2);
}
"

mod_mix_1 <- cmdstan_model(write_stan_file(stan_model_2))
fit_mix_1 <- mod_mix_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.95, parallel_chains = 4,
                      refresh = 100, max_treedepth = 10, 
                      thin = 1, init = 2, data = dat_1)

#check output
summ_mix_1 <- fit_mix_1$summary()
print(summ_mix_1[order(summ_mix_1$rhat, decreasing = T),])
samps_mix_1 <- data.frame(as_draws_df(fit_mix_1$draws()))
hist(samps_mix_1$theta, breaks = 0:40/40)

stan_model_3 <- "
data {
  int<lower=0> p;
  vector[p] b_xy_est;
  vector<lower=0>[p] b_xy_se;
  vector[p] b_xz_est;
  vector<lower=0>[p] b_xz_se;
}

parameters {
  // Mixture proportion
  real<lower=0, upper=1> theta;  // Probability of choosing Alternative 1

  // Parameters for Alternative 1
  vector<lower=-10, upper=10>[p] b_xy_raw_1;
  vector<lower=-10, upper=10>[p] b_xz_raw_1;
  real<lower=-10, upper=10> b_yz_1;
  real<lower=-10, upper=10> a_1;
  real<lower=0, upper=10> b_xy_sd_1;
  real<lower=0, upper=10> b_xz_sd_1;

  // Parameters for Alternative 2
  vector<lower=-10, upper=10>[p] b_xy_raw_2;
  vector<lower=-10, upper=10>[p] b_xz_raw_2;
  real<lower=-10, upper=10> b_yz_2;
  real<lower=-10, upper=10> a_2;
  real<lower=0, upper=10> b_xy_sd_2;
  real<lower=0, upper=10> b_xz_sd_2;
}

transformed parameters {
  // Transformed parameters for Alternative 1
  vector[p] b_xy_1 = b_xy_raw_1 * b_xy_sd_1;
  vector[p] b_xz_1 = b_xz_raw_1 * b_xz_sd_1 + a_1 + b_yz_1 * b_xy_1;

  // Transformed parameters for Alternative 2 (flipped)
  vector[p] b_xy_2 = b_xz_raw_2 * b_xz_sd_2; 
  vector[p] b_xz_2 = b_xy_raw_2 * b_xy_sd_2 + a_2 + b_yz_2 * b_xy_2;
}

model {
  real log_lik_1;
  real log_lik_2;
  real log_prior_1;
  real log_prior_2;
  
  // Priors for the mixture proportion
  theta ~ beta(1, 1);  // Uniform prior for mixture proportion
  
  // Log-priors for Alternative 1
  log_prior_1 = std_normal_lpdf(a_1) +
                std_normal_lpdf(b_xy_sd_1) +
                std_normal_lpdf(b_xz_sd_1) +
                std_normal_lpdf(b_yz_1 / 5) +
                std_normal_lpdf(b_xy_raw_1) +
                std_normal_lpdf(b_xz_raw_1);
  
  // Log-priors for Alternative 2
  log_prior_2 = std_normal_lpdf(a_2) +
                std_normal_lpdf(b_xy_sd_2) +
                std_normal_lpdf(b_xz_sd_2) +
                std_normal_lpdf(b_yz_2 / 5) +
                std_normal_lpdf(b_xy_raw_2) +
                std_normal_lpdf(b_xz_raw_2);
  
  // Log-likelihood for Alternative 1
  log_lik_1 = normal_lpdf(b_xy_est | b_xy_1, b_xy_se) + 
              normal_lpdf(b_xz_est | b_xz_1, b_xz_se);
              
  // Log-likelihood for Alternative 2 (flipped)
  log_lik_2 = normal_lpdf(b_xy_est | b_xy_2, b_xy_se) + 
              normal_lpdf(b_xz_est | b_xz_2, b_xz_se);
              
  Mixture model using log_sum_exp to combine log-likelihoods and log-priors
  target += log_sum_exp(log(theta) + log_prior_1 + log_lik_1,
                        log1m(theta) + log_prior_2 + log_lik_2);
}
"


# print(\"target = \", target());
#   print(\"theta = \", theta);
#   print(\"log_lik_1 = \", log_lik_1);
#   print(\"log_lik_2 = \", log_lik_2);
#   print(\"log_prior_1 = \", log_prior_1);
#   print(\"log_prior_2 = \", log_prior_2);
#   
# print(\"log_lik_1 = \", log_lik_1);
# print(\"log_lik_2 = \", log_lik_2);
# print(\"log_prior_1 = \", log_prior_1);
# print(\"log_prior_2 = \", log_prior_2);
# print(\"theta = \", theta);
# print(\"b_xy_sd_2 = \", b_xy_sd_2);
# print(\"b_xy_sd_1 = \", b_xy_sd_1);
# print(\"b_xz_sd_2 = \", b_xz_sd_2);
# print(\"b_xz_sd_1 = \", b_xz_sd_1);

mod_mix_2 <- cmdstan_model(write_stan_file(stan_model_3))
fit_mix_2 <- mod_mix_2$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                              adapt_delta = 0.95, parallel_chains = 4,
                              refresh = 100, max_treedepth = 10, 
                              thin = 1, init = 0.1, data = dat_1)

#check output
summ_mix_2 <- fit_mix_2$summary()
print(summ_mix_2[order(summ_mix_2$rhat, decreasing = T),])
samps_mix_2 <- data.frame(as_draws_df(fit_mix_2$draws()))
hist(samps_mix_2$theta, breaks = 0:40/40)
hist(samps_mix_2$b_yz_1)
hist(samps_mix_2$b_yz_2)
plot(samps_mix_2$b_yz_1, samps_mix_2$b_yz_2)

#### alternative joint spec ####
#specify z <- x -> y -> z inference model
dat <- list(n = n,
            p = p,
            x = x,
            y = y,
            z = z)

stan_model_1 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n,p] x;
  vector[n] y;
  vector[n] z;
}
parameters {
  
  //parameters
  vector[p] b_xy;
  vector[p] b_xz;
  real b_yz;
  
  //scale terms
  real<lower=0> b_xy_sd;
  real<lower=0> b_xz_sd;
  vector<lower=0>[p] x_sd;
  real<lower=0> y_sd;
  real<lower=0> z_sd;
  
}
transformed parameters {
  cov_matrix[2] S_yz;
  
  S_yz[1,1] = y_sd^2;
  S_yz[2,2] = z_sd^2;
  S_yz[1,2] = sum(x_sd^2 .* (b_xy .* b_xz + b_xy^2 * b_yz));
  S_yz[2,1] = S_yz[1,2];
  
}
model {
  // priors //

  //scale parameters
  x_sd ~ std_normal();
  y_sd ~ std_normal();
  z_sd ~ std_normal();
  b_xy_sd ~ normal(0,3);
  b_xz_sd ~ normal(0,3);
  
  //coefficients
  b_yz ~ std_normal();
  b_xy ~ normal(0, b_xy_sd);
  b_xz ~ normal(0, b_xz_sd);
  
  // model //
  for(j in 1:n){
    [y[j],z[j]] ~ multi_normal([0,0], S_yz);
  }
  
  matrix[2,2] Si_yz = inverse(S_yz);
  for(i in 1:p){
    
    //calculate conditional distribution for ith entry of x
    row_vector[2] Ss = [x_sd[i]^2 * b_xy[i], 
                       x_sd[i]^2 * (b_xz[i] + b_xy[i] * b_yz)];
    vector[n] xi_mu; 
    for(j in 1:n){
      xi_mu[j] = Ss * Si_yz * to_vector([y[j], z[j]]);
    }
    real xi_sd = x_sd[i]^2 - Ss * Si_yz * to_vector([Ss[2], Ss[1]]);
    
    // sample conditional distribution for x
    x[,i] ~ normal(xi_mu, xi_sd);
  }
  
}
"


mod_1 <- cmdstan_model(write_stan_file(stan_model_1))
fit_1 <- mod_1$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.9, parallel_chains = 4,
                      refresh = 10, max_treedepth = 10, 
                      thin = 1, init = 0.1, data = dat)

summ_1 <- fit_1$summary()
print(summ_1[order(summ_1$rhat, decreasing = T),])
samps_1 <- data.frame(as_draws_df(fit_1$draws()))
plot(samps_1$b_yz, type = "l")
hist(samps_1$b_yz)
