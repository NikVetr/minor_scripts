#packages
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))


#### simple binomial model ####

n_large <- 50
n_small <- 1000
total <- c(rep(50, n_large), rep(5, n_small))
logodds <- c(rnorm(n_large, 0, 2), rnorm(n_small, 0, 0.1))
probs <- invlogit(logodds)
count <- rbinom(n_large + n_small, total, prob = probs)


d <- list(count = count, 
          total = total, 
          n = n_large + n_small)

stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> total[n];
    int<lower=0> count[n];
}
parameters {
    vector[n] raw_logodds;
    real<lower=0> logodds_sd;
}
transformed parameters {
    vector[n] logodds = raw_logodds * logodds_sd;
}
model {
    raw_logodds ~ std_normal();
    logodds_sd ~ normal(0,2);
    counts ~ binomial_logit(total, logodds);
}
generated quantities {
    vector<lower=0, upper=1>[n] probs = inv_logit(logodds);
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)


#fit model
base = "multilevel_logodds"
write_stan_file(stan_program, dir = "~/Desktop/", basename = base)
write_stan_json(d, paste0("~/Desktop/", base,".json"))
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                  adapt_delta = 0.85, refresh = 100, init = 0.1, max_treedepth = 15, thin = 5)
# summ <- out$summary()
# summ[order(summ$ess_bulk),]
# summ[order(summ$rhat, decreasing = T),]
samps <- data.frame(as_draws_df(out$draws()))

par(mfrow = c(2,2))
hist(samps$logodds_sd)
plot(probs, apply((samps[,grep("probs", colnames(samps))]), 2, mean), 
     pch = 19, col = adjustcolor(c(rep("blue", n_large), rep("orange", n_small)),0.5), cex = 2,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "true probabilities",
     ylab = "posterior mean probabilities (logit-normal)")
abline(0,1, lwd = 2, lty = 2, col = adjustcolor(2,0.25))
abline(h = 0.5, lty = 3)
abline(v = 0.5, lty = 3)
plot(probs, (count + 1) / (total + 2), 
     pch = 19, col = adjustcolor(c(rep("blue", n_large), rep("orange", n_small)),0.5), cex = 2,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "true probabilities",
     ylab = "posterior mean probabilities (flat beta)")
abline(0,1, lwd = 2, lty = 2, col = adjustcolor(2,0.25))
abline(h = 0.5, lty = 3)
abline(v = 0.5, lty = 3)



#### row_n x col_n contingency table ####
#specify high-level parameters
row_n <- 10
col_n <- 15
total <- 5E3
row_prob <- invlogit(rnorm(row_n, -2, 0.5))
col_prob <- invlogit(rnorm(col_n, -1, 1))
row_count <- rbinom(n = row_n, size = total, prob = row_prob)
col_count <- rbinom(n = col_n, size = total, prob = col_prob)

#find bounds on cells
cell_maxes <- sapply(1:col_n, function(ci) sapply(1:row_n, function(ri) min(row_count[ri], col_count[ci])))
cell_exp_probs <- sapply(1:col_n, function(ci) sapply(1:row_n, function(ri) c(row_prob[ri], col_prob[ci])[which.max(c(row_count[ri], col_count[ci]))]))
cell_exp_logodds <- logit(cell_exp_probs)

#simulate bias parameters
bias_sds <- c(row = 0.25, col = 0.5, cell = 0.75)
overall_bias <- -0.35
row_bias <- rnorm(row_n,0,bias_sds["row"])
col_bias <- rnorm(col_n,0,bias_sds["col"])
cell_bias <- matrix(rnorm(row_n * col_n, 0, bias_sds["cell"]), row_n, col_n)

#propagate it out to cell counts
cell_biased_logodds <- t(t(cell_exp_logodds + row_bias) + col_bias) + cell_bias + overall_bias
cell_biased_probs <- invlogit(cell_biased_logodds)
cell_count <- t(sapply(1:row_n, function(ri) rbinom(n = col_n, size = cell_maxes[ri,], prob = cell_biased_probs[ri,])))
cell_expected_count <- t(t(row_prob)) %*% t(col_prob) * total

#pass data to list
d <- list(cell_count = c(cell_count),
          total = total,
          row_count = row_count,
          col_count = col_count,
          row_index = rep(1:row_n, col_n),
          col_index = unlist(lapply(1:col_n, function(i) rep(i, row_n))),
          row_n = row_n,
          col_n = col_n)

#do quick EDA
layout(matrix(c(1,1,1,1,2,3,4,5,6,7), ncol=2, byrow = T))
par(mar = c(4,5,4,2))
plot(y = range(cell_count), x = range(cell_expected_count),
     col = "white", ylab = "observed counts", xlab = "expected counts")
text(y = cell_count, x = cell_expected_count,
     col = rep(1:row_n, col_n), 
     labels = as.character(unlist(lapply(1:col_n, function(i) rep(i, row_n)))))
text(mean(par("usr")[1:2]), par("usr")[4], labels = "numbers represent column indices", pos = 3, xpd = NA)
legend(x = "bottomright", legend = paste0("row ", 1:row_n), col = 1:row_n, 
       pch = "X", ncol = 1, pt.cex = 1, cex = 0.75)
legend(x = "topleft", legend = "1-to-1 line", lty = 2, col = adjustcolor(1,0.5), lwd = 2)
abline(0,1, lty = 2, col = adjustcolor(1,0.5), lwd = 2)

#add in section showing poor calibration across columns when updated from a flat beta

basic_posteriors_90CI <- rbind(qbeta(p = 0.05, shape1 = 1 + d$cell_count, shape2 = 1 + d$total - d$cell_count),
                               qbeta(p = 0.95, shape1 = 1 + d$cell_count, shape2 = 1 + d$total - d$cell_count))

#now fit the model
# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=0> total;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> row_count[row_n];
#     int<lower=0> col_count[col_n];
#     int<lower=0> cell_count[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
# 
#     real row_mean;
#     real<lower=0> row_sd;
#     vector[row_n] raw_row_logodds;
#     
#     vector[n] raw_cell_logodds;
#     real<lower=0> cell_sd;
#     
#     //biases in deviations terms
#     vector[row_n] raw_row_bias;
#     vector[row_n] row_bias_log_sd;
#     vector[col_n] raw_col_bias;
#     vector[col_n] col_bias_log_sd;
#     
#     //multilevel deviation term params
#     real<lower=0> row_bias_sd;
#     real<lower=0> col_bias_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[row_n] row_logodds = raw_row_logodds * row_sd + row_mean;
#     vector[n] cell_mean = logit(inv_logit(col_logodds[col_index]) .* inv_logit(row_logodds[row_index]));
#     vector[n] cell_logodds = raw_cell_logodds * cell_sd .* exp(row_bias_log_sd[row_index]) .* exp(col_bias_log_sd[col_index]) +
#                                   cell_mean + raw_row_bias[row_index] * row_bias_sd + 
#                                   raw_col_bias[col_index] * col_bias_sd;
# }
# model {
#     //priors and hyperpriors
#     
#     //marginal params
#     raw_col_logodds ~ std_normal();
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     
#     raw_row_logodds ~ std_normal();
#     row_mean ~ normal(0,2);
#     row_sd ~ std_normal();
#     
#     //bias params
#     raw_row_bias ~ std_normal();
#     row_bias_log_sd ~ std_normal();
#     raw_col_bias ~ std_normal();
#     col_bias_log_sd ~ std_normal();
#     row_bias_sd ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_cell_logodds ~ std_normal();
#     cell_sd ~ std_normal();
#     
#     //likelihood
#     col_count ~ binomial_logit(total, col_logodds);
#     row_count ~ binomial_logit(total, row_logodds);
#     cell_count ~ binomial_logit(total, cell_logodds);
# }
# generated quantities {
#     vector[n] cell_bias = cell_logodds - cell_mean;
#     vector[n] cell_prob_bias = inv_logit(cell_logodds) - inv_logit(cell_mean);
#     vector[row_n] row_bias = raw_row_bias * row_bias_sd;
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
# }
# '


# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=0> total;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> row_count[row_n];
#     int<lower=0> col_count[col_n];
#     int<lower=0> cell_count[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
# 
#     real row_mean;
#     real<lower=0> row_sd;
#     vector[row_n] raw_row_logodds;
#     
#     vector[n] raw_cell_logodds;
#     real<lower=0> cell_sd;
#     
#     //biases in deviations terms
#     vector[row_n] raw_row_bias;
#     vector[col_n] raw_col_bias;
# 
#     //multilevel deviation term params
#     real<lower=0> row_bias_sd;
#     real<lower=0> col_bias_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[row_n] row_logodds = raw_row_logodds * row_sd + row_mean;
#     vector[n] cell_mean = logit(inv_logit(col_logodds[col_index]) .* inv_logit(row_logodds[row_index]));
#     vector[n] cell_logodds = raw_cell_logodds * cell_sd + cell_mean + 
#                              raw_row_bias[row_index] * row_bias_sd + 
#                              raw_col_bias[col_index] * col_bias_sd;
# }
# model {
#     //priors and hyperpriors
#     
#     //marginal params
#     raw_col_logodds ~ std_normal();
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     
#     raw_row_logodds ~ std_normal();
#     row_mean ~ normal(0,2);
#     row_sd ~ std_normal();
#     
#     //bias params
#     raw_row_bias ~ std_normal();
#     raw_col_bias ~ std_normal();
#     row_bias_sd ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_cell_logodds ~ std_normal();
#     cell_sd ~ std_normal();
#     
#     //likelihood
#     col_count ~ binomial_logit(total, col_logodds);
#     row_count ~ binomial_logit(total, row_logodds);
#     cell_count ~ binomial_logit(total, cell_logodds);
# }
# generated quantities {
#     vector[n] cell_bias = cell_logodds - (cell_mean + 
#                              raw_row_bias[row_index] * row_bias_sd + 
#                              raw_col_bias[col_index] * col_bias_sd);
#     vector[n] cell_total_prob_bias = inv_logit(cell_logodds) - inv_logit(cell_mean);
#     vector[row_n] row_bias = raw_row_bias * row_bias_sd;
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
# }
# '

# let's try the conditional model?
stan_program <- '
data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=0> total;
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=0> row_count[row_n];
    int<lower=0> col_count[col_n];
    int<lower=0> cell_count[row_n * col_n];
}
transformed data {
    int<lower=1> n = row_n * col_n;
    int<lower=0,upper=1> smaller_margin[n]; //0 if row, 1 if col
    int<lower=0> marginal_total[n];
    for(i in 1:n){
      //smaller_margin[i] = 1;
      smaller_margin[i] = row_count[row_index[i]] > col_count[col_index[i]];
      marginal_total[i] = smaller_margin[i] * col_count[col_index[i]] + (1-smaller_margin[i]) * row_count[row_index[i]];
    }
    vector[n] smaller_margin_vec = to_vector(smaller_margin);
}
parameters {
    //main parameters
    real col_mean;
    real<lower=0> col_sd;
    vector[col_n] raw_col_logodds;

    real row_mean;
    real<lower=0> row_sd;
    vector[row_n] raw_row_logodds;

    vector[n] raw_cell_logodds;
    real<lower=0> cell_sd;

    //biases in deviations terms
    real overall_bias;
    vector[row_n] raw_row_bias;
    vector[col_n] raw_col_bias;
    real<lower=0> row_bias_sd;
    real<lower=0> col_bias_sd;
}
transformed parameters {
    //recenter params
    vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
    vector[row_n] row_logodds = raw_row_logodds * row_sd + row_mean;
    vector[n] cell_mean = smaller_margin_vec .* row_logodds[row_index] + (1-smaller_margin_vec) .* col_logodds[col_index];
    vector[n] cell_logodds = raw_cell_logodds * cell_sd + cell_mean + overall_bias +
                             raw_row_bias[row_index] * row_bias_sd +
                             raw_col_bias[col_index] * col_bias_sd;
}
model {
    //priors and hyperpriors

    //marginal params
    raw_col_logodds ~ std_normal();
    col_mean ~ normal(0,2);
    col_sd ~ std_normal();

    raw_row_logodds ~ std_normal();
    row_mean ~ normal(0,2);
    row_sd ~ std_normal();

    //bias params
    overall_bias ~ std_normal();
    raw_row_bias ~ std_normal();
    raw_col_bias ~ std_normal();
    row_bias_sd ~ std_normal();
    col_bias_sd ~ std_normal();

    //cell params
    raw_cell_logodds ~ std_normal();
    cell_sd ~ std_normal();

    //likelihood
    col_count ~ binomial_logit(total, col_logodds);
    row_count ~ binomial_logit(total, row_logodds);
    cell_count ~ binomial_logit(marginal_total, cell_logodds);
}
generated quantities {
    vector[n] cell_bias = cell_logodds - (cell_mean + overall_bias + 
                             raw_row_bias[row_index] * row_bias_sd +
                             raw_col_bias[col_index] * col_bias_sd);
    vector[n] cell_total_prob_bias = inv_logit(cell_logodds) - inv_logit(cell_mean);
    vector[row_n] row_bias = raw_row_bias * row_bias_sd;
    vector[col_n] col_bias = raw_col_bias * col_bias_sd;
}
'

# stan_program <- '
# data {
#     int<lower=1> row_n;
#     int<lower=1> col_n;
#     int<lower=0> total;
#     int<lower=1,upper=row_n> row_index[row_n * col_n];
#     int<lower=1,upper=col_n> col_index[row_n * col_n];
#     int<lower=0> row_count[row_n * col_n];
#     int<lower=0> col_count[row_n * col_n];
#     int<lower=0> cell_count[row_n * col_n];
# }
# transformed data {
#     int<lower=1> n = row_n * col_n;
#     int<lower=0,upper=1> smaller_margin[n]; //0 if row, 1 if col
#     int<lower=0> marginal_total[n];
#     for(i in 1:n){
#       //smaller_margin[i] = 1;
#       smaller_margin[i] = row_count[i] > col_count[i];
#       marginal_total[i] = smaller_margin[i] * col_count[i] + (1-smaller_margin[i]) * row_count[i];
#     }
#     vector[n] smaller_margin_vec = to_vector(smaller_margin);
# }
# parameters {
#     //main parameters
#     real col_mean;
#     real<lower=0> col_sd;
#     vector[col_n] raw_col_logodds;
#     real<lower=0> col_indiv_sd;
#     vector[n] raw_col_indiv_logodds;
# 
#     real row_mean;
#     real<lower=0> row_sd;
#     vector[row_n] raw_row_logodds;
#     real<lower=0> row_indiv_sd;
#     vector[n] raw_row_indiv_logodds;
#     
#     vector[n] raw_cell_logodds;
#     real<lower=0> cell_sd;
#     
#     //biases in deviations terms
#     real overall_bias;
#     vector[row_n] raw_row_bias;
#     vector[col_n] raw_col_bias;
#     real<lower=0> row_bias_sd;
#     real<lower=0> col_bias_sd;
# }
# transformed parameters {
#     //recenter params
#     vector[col_n] col_logodds = raw_col_logodds * col_sd + col_mean;
#     vector[n] col_indiv_logodds = raw_col_indiv_logodds * col_indiv_sd + col_logodds[col_index];
#     vector[row_n] row_logodds = raw_row_logodds * row_sd + row_mean;
#     vector[n] row_indiv_logodds = raw_row_indiv_logodds * row_indiv_sd + row_logodds[row_index];
#     vector[n] cell_mean = smaller_margin_vec .* row_indiv_logodds + (1-smaller_margin_vec) .* col_indiv_logodds;
#     vector[n] cell_logodds = raw_cell_logodds * cell_sd + cell_mean + 
#                              overall_bias +
#                              raw_row_bias[row_index] * row_bias_sd + 
#                              raw_col_bias[col_index] * col_bias_sd;
# }
# model {
#     //priors and hyperpriors
#     
#     //marginal params
#     raw_col_logodds ~ std_normal();
#     raw_col_indiv_logodds ~ std_normal();
#     col_mean ~ normal(0,2);
#     col_sd ~ std_normal();
#     col_indiv_sd ~ std_normal();
#     
#     raw_row_logodds ~ std_normal();
#     raw_row_indiv_logodds ~ std_normal();
#     row_mean ~ normal(0,2);
#     row_sd ~ std_normal();
#     row_indiv_sd ~ std_normal();
#     
#     //bias params
#     overall_bias ~ std_normal();
#     raw_row_bias ~ std_normal();
#     raw_col_bias ~ std_normal();
#     row_bias_sd ~ std_normal();
#     col_bias_sd ~ std_normal();
#     
#     //cell params
#     raw_cell_logodds ~ std_normal();
#     cell_sd ~ std_normal();
#     
#     //likelihood
#     col_count ~ binomial_logit(total, col_indiv_logodds);
#     row_count ~ binomial_logit(total, row_indiv_logodds);
#     cell_count ~ binomial_logit(marginal_total, cell_logodds);
# }
# generated quantities {
#     vector[n] cell_bias = cell_logodds - (cell_mean + overall_bias +
#                              raw_row_bias[row_index] * row_bias_sd + 
#                              raw_col_bias[col_index] * col_bias_sd);
#     vector[n] cell_total_prob_bias = inv_logit(cell_logodds) - inv_logit(cell_mean);
#     vector[n] cell_total_prob_bias_logodds = cell_logodds - cell_mean;
#     vector[row_n] row_bias = raw_row_bias * row_bias_sd;
#     vector[col_n] col_bias = raw_col_bias * col_bias_sd;
# }
# '

#compile model
if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
base = "multilevel_logodds"
write_stan_file(stan_program, dir = "~/Desktop/", basename = base)
write_stan_json(d, paste0("~/Desktop/", base,".json"))
out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = d, parallel_chains = 4, 
                  adapt_delta = 0.9, refresh = 50, init = 0.1, max_treedepth = 20, thin = 5)

#quick check diagnostics
summ <- out$summary()
summ[order(summ$ess_bulk),]
summ[order(summ$rhat, decreasing = T),]

#examine output
samps <- data.frame(as_draws_df(out$draws()))
subset_samps <- function(include, exclude, samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}
seq2 <- function(lims, lout){
  seq(lims[1], lims[2], length.out = lout)
}
coverage_prop <- 0.9
coverage <- apply(samps, 2, quantile, probs = c((1-coverage_prop)/2,coverage_prop+(1-coverage_prop)/2))

cell_coverage <- round(mean(subset_samps("cell_bias", "raw", coverage)[1,] < c(cell_bias) & subset_samps("cell_bias", "raw", coverage)[2,] > c(cell_bias)), 2)
row_coverage <- round(mean(subset_samps("row_bias", c("raw", "sd"), coverage)[1,] < c(row_bias) & subset_samps("row_bias", c("raw", "sd"), coverage)[2,] > c(row_bias)), 2)
col_coverage <- round(mean(subset_samps("col_bias", c("raw", "sd"), coverage)[1,] < c(col_bias) & subset_samps("col_bias", c("raw", "sd"), coverage)[2,] > c(col_bias)), 2)

#do some quick plotting
plot(cell_bias, apply(subset_samps("cell_bias", "raw", samps), 2, mean), pch = 19, col = adjustcolor(1, 0.35), cex = 1.5,
     xlab = "true cell bias (logodds)", ylab = "posterior mean cell bias (logodds)",
     main = paste0(coverage_prop*100, "% coverage = ", cell_coverage))
abline(0,1, lty = 2, col = adjustcolor(2,0.75), lwd = 2)
hist(samps$cell_sd, probability = T, main = "Cell Bias SD", breaks = seq2(range(c(bias_sds["cell"], samps$cell_sd)), 10))
abline(v = bias_sds["cell"], col = 2, lty = 2, lwd = 2)

plot(row_bias, apply(subset_samps("row_bias", c("raw", "sd"), samps), 2, mean), pch = 19, col = adjustcolor(1, 0.35), cex = 1.5,
     xlab = "true row bias (logodds)", ylab = "posterior mean row bias (logodds)",
     main = paste0(coverage_prop*100, "% coverage = ", row_coverage))
abline(0,1, lty = 2, col = adjustcolor(2,0.75), lwd = 2)
hist(samps$row_bias_sd, probability = T, main = "Row Bias SD", breaks = seq2(range(c(bias_sds["row"], samps$row_bias_sd)), 10))
abline(v = bias_sds["row"], col = 2, lty = 2, lwd = 2)

plot(col_bias, apply(subset_samps("col_bias", c("raw", "sd"), samps), 2, mean), pch = 19, col = adjustcolor(1, 0.35), cex = 1.5,
     xlab = "true col bias (logodds)", ylab = "posterior mean col bias (logodds)",
     main = paste0(coverage_prop*100, "% coverage = ", col_coverage))
abline(0,1, lty = 2, col = adjustcolor(2,0.75), lwd = 2)
hist(samps$col_bias_sd, probability = T, main = "Col Bias SD", breaks = seq2(range(c(bias_sds["col"], samps$col_bias_sd)), 10))
abline(v = bias_sds["col"], col = 2, lty = 2, lwd = 2)



