library(cmdstanr)
library(posterior)

#### normal ####

#simulate data
set.seed(1)
n  <- 1000
p1  <- 0.95
p <- c(p1, 1-p1)
mu <- c(0, 0)      # means
sd <- c(1, 5)      # standard deviations
comp <- sample(1:2, size = n, replace = TRUE, prob = p) #component indices
y <- rnorm(n, mean = mu[comp], sd = sd[comp]) #observations

#specify Stan model
stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
}4
parameters {
  positive_ordered[2] sigma;
  real<lower=0, upper=1> theta;
}
model {
  sigma ~ normal(0,10);
  theta ~ beta(5,1);

  // Likelihood
  for (i in 1:N) {
    target += log_mix(theta,
                      normal_lpdf(y[i] | 0, sigma[1]),
                      normal_lpdf(y[i] | 0, sigma[2]));
  }
}
generated quantities {
  // Probability that each y[i] belongs to the second component
  vector[N] p_comp2;
  for (i in 1:N) {
    real comp1 = theta * exp(normal_lpdf(y[i] | 0, sigma[1]));
    real comp2 = (1 - theta) * exp(normal_lpdf(y[i] | 0, sigma[2]));
    p_comp2[i] = comp2 / (comp1 + comp2);
  }
}
"
stan_file_path <- write_stan_file(stan_code)
mod <- cmdstan_model(stan_file_path)

#fit model
fit <- mod$sample(
  data = list(N = n, y = y),
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  thin = 1
)

#check MCMC diagnostics
summ <- fit$summary(c("sigma", "theta"))
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#inspect posterior samples
samps <- data.frame(as_draws_df(fit$draws()))

#scale parameters
c(mean(samps$sigma.1.), mean(samps$sigma.2.))
sd

#population mixture probability
mean(samps$theta)
p1

#individual mixture probabilities
p_comp2 <- samps[,paste0("p_comp2.", 1:n, ".")]
mean_p_comp2 <- colMeans(p_comp2)
true_comp1 <- mean_p_comp2[comp==1]
true_comp2 <- mean_p_comp2[comp==2]
hist(true_comp1, breaks = 0:50/50, col = adjustcolor("darkblue", 0.5), freq = F)
hist(true_comp2, breaks = 0:50/50, col = adjustcolor("darkorange", 0.5), add = T, freq = F)

#compare to single distribution z-score
z <- scale.default(y)
plot(z, colMeans(p_comp2), pch = 19, col = comp)
pvals <- 1 - abs(0.5 - pnorm(z)) * 2
pvals_adj <- p.adjust(pvals, method = "BY")
par(mar = c(4,4,8,2))
mean_p_comp2_jittered <- mean_p_comp2 + rnorm(n)/100
mean_p_comp2_jittered[mean_p_comp2_jittered > 1] <- 1
pvals_adj_jittered <- pvals_adj + rnorm(n)/100
pvals_adj_jittered[pvals_adj_jittered > 1] <- 1
plot(mean_p_comp2_jittered, pvals_adj_jittered, col = adjustcolor(c("darkblue", "darkorange"), 0.5)[comp], 
     xlab = "outlier posterior probability (mixture model, jittered)", 
     ylab = "outlier BY-adjusted p-value (jittered)")
abline(h = 0.05, lty = 2)
abline(v = 0.95, lty = 2)
pusr <- par("usr")
wh <- c(diff(pusr[1:2]), diff(pusr[3:4]))
leg <- legend(x = pusr[1], y = pusr[4], col = adjustcolor(c("darkblue", "darkorange", "black"), 0.5), pch = c(19,19,NA), 
       legend = c("true component 1 (common, low variance)", "true component 2 (rare, high variance)", "probability threshold"), 
       lty = c(NA,NA,2), xpd = NA, title = "normal outlier calling", plot = F)
legend(x = pusr[1], y = pusr[4] + leg$rect$h, col = adjustcolor(c("darkblue", "darkorange", "black"), 0.5), pch = c(19,19,NA), 
       legend = c("true component 1 (common, low variance)", "true component 2 (rare, high variance)", "probability threshold"), 
       lty = c(NA,NA,2), xpd = NA, title = "normal mixture outlier calling", title.font = 2)
comp2_thresh <- data.frame(pp = mean_p_comp2[comp==2] > 0.95, pval = pvals_adj[comp==2] < 0.05)
text(pusr[2] - wh[1]/40, pusr[3] + wh[2]/20, sum(comp2_thresh$pp & comp2_thresh$pval), font = 2, col = 1)
text(pusr[2] - wh[1]/40, pusr[4] - wh[2]/20, sum(comp2_thresh$pp & !comp2_thresh$pval), font = 2, col = 1)
text(pusr[1] + wh[1]/40, pusr[4] - wh[2]/20, sum(!comp2_thresh$pp & !comp2_thresh$pval), font = 2, col = 1)

#### gampois ####
set.seed(1)
n <- 100000
mu <- rlnorm(n, meanlog = 5, sdlog = 1)
p1 <- 0.95
p <- c(p1, 1-p1)
iphi <- c(0.1,20)
comp <- sample(1:2, size = n, replace = TRUE, prob = p) #component indices
y <- rnbinom(n, size = 1/iphi[comp], mu = mu)

stan_code_gpois <- '
data {
  int<lower=1> N;        // number of observations
  array[N] int<lower=0> y;     // observed counts
  vector<lower=0>[N] mu; // known means for each observation
}
parameters {
  positive_ordered[2] iphi;
  real<lower=0, upper=1> theta;
}
model {
  iphi ~ normal(0,1);
  theta ~ beta(5,1);

  // Mixture likelihood
  for (i in 1:N) {
    // NB_2 parameterization in Stan: neg_binomial_2_lpmf(y | mu, phi)
    // E[y]=mu, Var[y] = mu + mu^2/phi
    target += log_mix(theta,
      neg_binomial_2_lpmf(y[i] | mu[i], 1/iphi[1]),
      neg_binomial_2_lpmf(y[i] | mu[i], 1/iphi[2])
    );
  }
}
generated quantities {
  vector[N] p_comp2;  // posterior probability each y[i] from 2nd component
  for (i in 1:N) {
    real comp1 = theta * exp(neg_binomial_2_lpmf(y[i] | mu[i], 1/iphi[1]));
    real comp2 = (1 - theta) * exp(neg_binomial_2_lpmf(y[i] | mu[i], 1/iphi[2]));
    p_comp2[i] = comp2 / (comp1 + comp2);
  }
}
'

# Write the Stan code to a file, then compile
stan_file_path_gpois <- write_stan_file(stan_code_gpois)
mod_gpois <- cmdstan_model(stan_file_path_gpois)

#fit model
fit_gpois <- mod_gpois$sample(
  data = list(N = n, y = y, mu = mu),
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000
)

#check MCMC diagnostics
summ <- fit_gpois$summary(c("iphi", "theta"))
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#inspect posterior samples
samps <- data.frame(as_draws_df(fit_gpois$draws()))

#scale parameters
c(mean(samps$iphi.1.), mean(samps$iphi.2.))
iphi

#population mixture probability
mean(samps$theta)
p1

#individual mixture probabilities
p_comp2 <- samps[,paste0("p_comp2.", 1:n, ".")]
mean_p_comp2 <- colMeans(p_comp2)
true_comp1 <- mean_p_comp2[comp==1]
true_comp2 <- mean_p_comp2[comp==2]
hist(true_comp1, breaks = 0:50/50, col = adjustcolor("darkblue", 0.5), freq = F)
hist(true_comp2, breaks = 0:50/50, col = adjustcolor("darkorange", 0.5), add = T, freq = F)

#compare to single distribution "z-score"
sample_iphi <- as.numeric(MASS::theta.ml(y, mu))
pvals <- 1-abs(0.5-pnbinom(y, size = sample_iphi, mu = mu)) * 2
pvals_adj <- p.adjust(pvals, method = "BY")
par(mar = c(4,4,8,2))
mean_p_comp2_jittered <- mean_p_comp2 + rnorm(n)/100
mean_p_comp2_jittered[mean_p_comp2_jittered > 1] <- 1
pvals_adj_jittered <- pvals_adj + rnorm(n)/100
pvals_adj_jittered[pvals_adj_jittered > 1] <- 1
plot(mean_p_comp2_jittered, pvals_adj_jittered, col = adjustcolor(c("darkblue", "darkorange"), 0.5)[comp], 
     xlab = "outlier posterior probability (mixture model, jittered)", 
     ylab = "outlier BY-adjusted p-value (jittered)")
abline(h = 0.05, lty = 2)
abline(v = 0.95, lty = 2)
pusr <- par("usr")
wh <- c(diff(pusr[1:2]), diff(pusr[3:4]))
leg <- legend(x = pusr[1], y = pusr[4], col = adjustcolor(c("darkblue", "darkorange", "black"), 0.5), pch = c(19,19,NA), 
              legend = c("true component 1 (common, low variance)", "true component 2 (rare, high variance)", "probability threshold"), 
              lty = c(NA,NA,2), xpd = NA, title = "normal mixture outlier calling", plot = F)
legend(x = pusr[1], y = pusr[4] + leg$rect$h, col = adjustcolor(c("darkblue", "darkorange", "black"), 0.5), pch = c(19,19,NA), 
       legend = c("true component 1 (common, low variance)", "true component 2 (rare, high variance)", "probability threshold"), 
       lty = c(NA,NA,2), xpd = NA, title = "neg-binomial outlier calling", title.font = 2)
comp2_thresh <- data.frame(pp = mean_p_comp2[comp==2] > 0.95, pval = pvals_adj[comp==2] < 0.05)
text(pusr[2] - wh[1]/40, pusr[3] + wh[2]/20, sum(comp2_thresh$pp & comp2_thresh$pval), font = 2, col = 1)
text(pusr[2] - wh[1]/40, pusr[4] - wh[2]/20, sum(comp2_thresh$pp & !comp2_thresh$pval), font = 2, col = 1)
text(pusr[1] + wh[1]/40, pusr[4] - wh[2]/20, sum(!comp2_thresh$pp & !comp2_thresh$pval), font = 2, col = 1)

