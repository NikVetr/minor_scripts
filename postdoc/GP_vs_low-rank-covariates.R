library(cmdstanr)
library(posterior)
library(data.table)
library(ape)
library(phangorn)

#read helper functions
source("~/repos/Stan2R/R/functions.R")

generate_kde_plot <- function(x, cols = c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD')) {
  
  # Check if x is a matrix or a list and process accordingly
  if (is.matrix(x)) {
    p <- nrow(x)
    x_list <- split(x, row(x))
  } else if (is.list(x)) {
    p <- length(x)
    x_list <- x
  } else {
    stop("Input must be a matrix or a list")
  }
  
  # Compute parameter ranges and densities
  param_ranges <- do.call(rbind, lapply(x_list, function(xi) c(min(xi), quantile(xi, 0.999))))
  params_range <- c(min(param_ranges[, 1]), max(param_ranges[, 2]))
  densities <- lapply(x_list, function(xi) density(xi, from = params_range[1], to = params_range[2]))
  density_maxima <- sapply(densities, function(dens) max(dens$y))
  density_maxima_locs <- sapply(densities, function(dens) dens$x[which.max(dens$y)])
  
  # Set up plotting parameters
  par(mar = c(4, 6, 5, 1))
  if(is.null(colnames(x))){
    kde_labels <- paste0("Var_", 1:ncol(x))
  } else {
    kde_labels <- colnames(x)  
  }
  
  label_locs_x <- seq(from = params_range[1], to = params_range[2], length.out = p)
  label_loc_y <- max(density_maxima) * 1.05
  label_order <- rank(density_maxima_locs, ties.method = "first")
  plot_cols <- cols[c(1:min(p, length(cols)), rep(length(cols), p - min(p, length(cols))))]
  
  # Plot setup
  plot(x = NULL, y = NULL, xlim = params_range, ylim = c(0, max(density_maxima)), ylab = "density", xlab = "Parameter", frame = FALSE)
  
  # Plot each density
  for (i in seq_along(densities)) {
    polygon(x = c(densities[[i]]$x, rev(densities[[i]]$x)), 
            y = c(densities[[i]]$y, rep(0, length(densities[[i]]$y))), 
            col = adjustcolor(plot_cols[label_order[i]], 0.05), border = plot_cols[label_order[i]], xpd = NA)
    
    text(labels = kde_labels[i], x = label_locs_x[label_order[i]] - strheight(paste("Var", i)) / 1.5 * (par()$usr[2] - par()$usr[1]) / (par()$usr[4] - par()$usr[3]), 
         y = label_loc_y, srt = 90, pos = 4, xpd = NA, col = plot_cols[label_order[i]])
    segments(x0 = density_maxima_locs[i], x1 = label_locs_x[label_order[i]], 
             y0 = density_maxima[i], y1 = label_loc_y, lty = 2, col = plot_cols[label_order[i]], xpd = NA)
  }
}

#generate correlation matrix for obs
rop <- 0
n <- 100
pow_scale <- 50
rs <- c(1, rop^((1:n)/pow_scale))
R <- outer(1:n, 1:n, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)

#specify expectation model
p <- 5
b <- 2
bs <- rep(b, p)
x <- t(chol(R)) %*% matrix(rnorm(n * p), n, p)
y_exp <- x %*% t(t(bs))

#specify noise model
e_sd <- sqrt(p^2 * 1.5) / 10
dist_matrix <- 1-R
nPCs <- 20

# eigsystem <- eigen(R)
# V <- eigsystem$vectors
# U <- eigsystem$values
# PCs <- matrix(rnorm(n * n), n, n) %*% eigsystem$vectors[,1:nPCs]
# PCs <- matrix(rnorm(n * n), n, n) %*% V[,1:nPCs]

n_PCA_obs <- 1E3
PCA_out <- prcomp((t(chol(R)) %*% matrix(rnorm(n * n_PCA_obs), n, n_PCA_obs)))
U <- PCA_out$sdev^2
V <- PCA_out$rotation
PCs <- apply(PCA_out$x[,1:nPCs], 2, scale)

par(mfrow = c(1,1))
plot(upgma(dist_matrix)) 
paste0("Proportion variance captured by first ", nPCs, " PCs: ", round(sum(U[1:nPCs]) / sum(U), 3))
e <- t(chol(R)) %*% matrix(rnorm(n) * e_sd, n, 1)
y <- y_exp + e

#fit conventional regression model
fit_lm <- lm(y ~ x)
summary(fit_lm)

#fit GP model
stan_model_GP <- "
data {
  int<lower=1> n;            // number of data points
  matrix[n, n] dist_matrix;  // precomputed distance matrix
  int<lower=1> p;            // number of predictors
  matrix[n, p] x;            // input observations
  vector[n] y;               // output observations
}

parameters {
  //exp model parameters
  real a;
  vector[p] b;

  //cov matrix parameters
  real<lower=0> length_scale; // length scale parameter for RBF kernel
  real<lower=0> sigma_f;      // signal variance for RBF kernel
  real<lower=0> sigma_n;      //iid noise term
}

transformed parameters {
  matrix[n, n] covmat;
  for (i in 1:n) {
    for (j in 1:n) {
      covmat[i, j] = sigma_f^2 * exp(-0.5 * dist_matrix[i, j]^2 / length_scale^2);
    }
    covmat[i, i] = covmat[i, i] + sigma_n^2;
  }

}

model {
  // Priors
  length_scale ~ std_normal();
  sigma_f ~ std_normal();
  sigma_n ~ std_normal();
  a ~ std_normal();
  b ~ std_normal();

  // Likelihood
  y ~ multi_normal_cholesky(rep_vector(a, n) + x * b, cholesky_decompose(covmat));
}
"

dat_GP <- list(n=n, 
            dist_matrix=dist_matrix,
            p=p,
            x=x, 
            y=c(y))

mod_GP <- cmdstan_model(write_stan_file(stan_model_GP))
fit_GP <- mod_GP$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                     data = dat_GP, adapt_delta = 0.9, parallel_chains = 4,
                     refresh = 100, max_treedepth = 15,
                     thin = 1)
focal_variables <- c("lp__", "a", "b", "length_scale", "sigma_f", "sigma_n")
summ_GP <- fit_GP$summary(variables = focal_variables)
summ_GP[order(summ_GP$ess_bulk),]
summ_GP[order(summ_GP$rhat, decreasing = T),]


#fit PCA covariates model
stan_model_PCA <- "
data {
  int<lower=1> n;            // number of data points
  int<lower=1> nPCs;         // number of PC axes included
  matrix[n, nPCs] PCs;       // PCs / covariates
  int<lower=1> p;            // number of predictors
  matrix[n, p] x;            // input observations
  vector[n] y;               // output observations
}

parameters {
  real a;
  vector[p] b;
  vector[nPCs] b_PCs;
  real<lower=0> sigma_e;
}

transformed parameters {
  
}

model {
  // priors
  sigma_e ~ std_normal();
  a ~ std_normal();
  b ~ std_normal();
  b_PCs ~ std_normal();

  // likelihood
  y ~ normal(rep_vector(a, n) + x * b + PCs * b_PCs, sigma_e);
}
"

dat_PCA <- list(n=n, 
               nPCs = nPCs,
               PCs = PCs,
               p=p,
               x=x, 
               y=c(y))

mod_PCA <- cmdstan_model(write_stan_file(stan_model_PCA))
fit_PCA <- mod_PCA$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                     data = dat_PCA, adapt_delta = 0.85, parallel_chains = 4,
                     refresh = 100, max_treedepth = 15, 
                     thin = 2)
summ_PCA <- fit_PCA$summary()
summ_PCA[order(summ_PCA$ess_bulk),]
summ_PCA[order(summ_PCA$rhat, decreasing = T),]

#extract parameter posteriors
samps_GP <- data.frame(as_draws_df(fit_GP$draws(variables = focal_variables)))
samps_PCA <- data.frame(as_draws_df(fit_PCA$draws()))
est_bs_GP <- subset_samps("b", as.data.table(samps_GP))
est_bs_PCA <- subset_samps("b", as.data.table(samps_PCA))

#visualize to compare to true values
par(mfrow = c(2,1))

#first draw marginal density estimates for GP model
generate_kde_plot(est_bs_GP)
abline(v = bs, lwd = 3)

#then for PCA model
generate_kde_plot(est_bs_PCA)
abline(v = bs, lwd = 3)

