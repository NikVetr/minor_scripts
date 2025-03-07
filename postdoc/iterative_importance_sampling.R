#### libraries ####
library(loo)
library(mixtools)
library(MASS)
library(cmdstanr)
library(posterior)

stan_code <- "
data {
  vector[2] blob_x;
  vector[2] blob_y;
  vector[3] probs;
}
parameters {
  real x;
  real y;
}
model {
  // Weights for the mixture
  real log_p_banana = log(probs[1]);
  real log_p_blob   = log(probs[2]);

  // 1) Banana log-density
  real lp_banana = normal_lpdf(x | 0, 1)
                 + normal_lpdf(y | x^2, 1);

  // 2) Left blob log-density
  real lp_blob1  = normal_lpdf(x | blob_x[1], 0.5)
                 + normal_lpdf(y | blob_y[1], 0.5);

  // 3) Right blob log-density
  real lp_blob2  = normal_lpdf(x |  blob_x[2], 0.5)
                 + normal_lpdf(y | blob_y[2], 0.5);

  // Mixture:
  // target += log_sum_exp( log(weight1) + lp_1, ..., log(weightN) + lp_N );
  vector[3] lp;
  lp[1] = log_p_banana + lp_banana;
  lp[2] = log_p_blob   + lp_blob1;
  lp[3] = log_p_blob   + lp_blob2;
  target += log_sum_exp(lp);
}
"

# Write code to a file "banana.stan"
model_path <- "~/scripts/minor_scripts/postdoc/evil_smiley.stan"
mod <- cmdstan_model(model_path)

#### functions #### 
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")

log_posterior <- function(theta, dat) {
  
  # Component densities
  lp_banana <- dnorm(theta[,1], dat$banana_mu[1], dat$banana_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$banana_mu[2] + theta[,1]^2, dat$banana_sd[2], log=TRUE)
  
  lp_blob1 <- dnorm(theta[,1], dat$blob_1_mu[1], dat$blob_1_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$blob_1_mu[2], dat$blob_1_sd[2], log=TRUE)
  
  lp_blob2 <- dnorm(theta[,1], dat$blob_2_mu[1], dat$blob_2_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$blob_2_mu[2], dat$blob_2_sd[2], log=TRUE)
  
  # logsumexp over these
  matrixStats::rowLogSumExps(cbind(log(dat$probs[1]) + lp_banana, 
                                   log(dat$probs[2]) + lp_blob1, 
                                   log(dat$probs[3]) + lp_blob2))
}

# 
# log_posterior <- function(theta, blob_x, blob_y, probs) {
#   x <- theta[,1]
#   y <- theta[,2]
#   
#   # Banana component
#   lp_banana <- dnorm(x, 0, 1, log=TRUE) + dnorm(y, x^2, 1, log=TRUE)
#   
#   # Left blob near (-2, -2), scale=0.3
#   lp_blob1  <- dnorm(x, blob_x[1], 0.5, log=TRUE) + dnorm(y, blob_y[1], 0.5, log=TRUE)
#   
#   # Right blob near (2, -2), scale=0.3
#   lp_blob2  <- dnorm(x,  blob_x[2], 0.5, log=TRUE) + dnorm(y, blob_y[2], 0.5, log=TRUE)
#   
#   # Mixture weights: 0.8 banana, 0.1 blob1, 0.1 blob2
#   comp1 <- log(probs[1]) + lp_banana
#   comp2 <- log(probs[2]) + lp_blob1
#   comp3 <- log(probs[3]) + lp_blob2
#   
#   #compute final density
#   matrixStats::rowLogSumExps(cbind(comp1, comp2, comp3))
# }


birth_component <- function(old_fit, data, new_weight = 0.05, up = T) {
  
  old_lambda  <- old_fit$lambda
  old_mu      <- old_fit$mu
  old_sigma   <- old_fit$sigma
  k           <- length(old_lambda)
  
  if(!up){
    best_components <- order(old_fit$lambda, decreasing = T)[1:2]
    new_lambda_vec <- old_lambda[best_components] / sum(old_lambda[best_components])
    new_mu_list    <- old_mu[best_components]
    new_sigma_list <- old_sigma[best_components]
    
    out <- list(lambda = new_lambda_vec,
                mu     = new_mu_list,
                sigma  = new_sigma_list)
    
  } else {
    
    # 1) Compute mixture PDF at each point under old_fit:
    pdf_vals <- numeric(nrow(data))
    for (j in seq_len(k)) {
      pdf_vals <- pdf_vals + old_lambda[j] *
        dmvnorm(data, mu = old_mu[[j]], sigma = old_sigma[[j]])
    }
    
    # 2) "Inverse-density" weighting => highlight underfit regions
    inv_weights <- 1 / (pdf_vals + 1e-12)
    inv_weights <- inv_weights / sum(inv_weights)
    
    # 3) Weighted average of data => new mean
    new_mu <- colSums(data * inv_weights)
    
    # 4) Simple covariance guess for the new component (identity)
    new_sigma <- diag(1, 2)
    
    # 5) Adjust old weights to sum to (1 - new_weight)
    lam_sum <- sum(old_lambda)
    old_lambda_scaled <- old_lambda * ((1 - new_weight) / lam_sum)
    new_lambda_vec <- c(old_lambda_scaled, new_weight)
    
    # Combine old means/sigmas with the new component
    new_mu_list    <- c(old_mu,    list(as.numeric(new_mu)))
    print(names(new_mu_list))
    new_sigma_list <- c(old_sigma, list(new_sigma))
    
    out <- list(lambda = new_lambda_vec,
         mu     = new_mu_list,
         sigma  = new_sigma_list)
    
  }
  
  return(out)
}


sample_from_mix <- function(n, gm_fit) {
  # gm_fit: result of mvnormalmixEM
  # randomly assign each new draw to a component j with prob=gm_fit$lambda[j]
  
  k  <- length(gm_fit$lambda)
  z  <- sample.int(k, size = n, replace = TRUE, prob = gm_fit$lambda)
  
  out <- matrix(NA_real_, n, 2)
  idx_current <- 1L
  for (j in seq_len(k)) {
    nj <- sum(z == j)
    if (nj > 0) {
      out[z == j, ] <- mvrnorm(nj, gm_fit$mu[[j]], gm_fit$sigma[[j]])
    }
  }
  colnames(out) <- c("x","y")
  out
}

plot_density_grid <- function(xy_sub, cell_dim, cols, xlim = NULL, ylim = NULL){
  plot.new()
  if(is.null(xlim)){
    xlim <- range(xy_sub$x)
  }
  if(is.null(ylim)){
    ylim <- range(xy_sub$y)
  }
  plot.window(xlim = xlim, ylim = ylim)
  axis(1); axis(2)
  mtext("x", 1, line = 2.5, cex = 2, font = 2)
  mtext("y", 2, line = 2.5, cex = 2, font = 2, )
  for(i in 1:length(cols)){
    rect(xleft = xy_sub$x[i] - cell_dim$w/2,
         ybottom = xy_sub$y[i] - cell_dim$h/2,
         xright = xy_sub$x[i] + cell_dim$w/2,
         ytop = xy_sub$y[i] + cell_dim$h/2,
         col = cols[i], border = NA
    )
  }  
}

#fit mixture model
compute_bic <- function(mixfit, data) {
  N    <- nrow(data)
  k    <- length(mixfit$lambda)   # number of components
  logL <- mixfit$loglik          # final log-likelihood from EM
  p    <- (k - 1) + 2*k + 3*k     # = 6k - 1
  bic  <- -2*logL + p*log(N)
  bic
}

log_mix_density <- function(theta, mixfit) {
  pi_vec    <- mixfit$lambda
  mu_list   <- mixfit$mu
  sigma_list<- mixfit$sigma
  k <- length(pi_vec)
  
  # We'll compute for each row i: sum_j pi_j * exp( logN(theta_i|mu_j,Sigma_j) )
  # Then take log of that sum.
  n <- nrow(theta)
  comp_vals <- matrix(NA, n, k)
  for (j in seq_len(k)) {
    comp_vals[,j] <- pi_vec[j] * dmvnorm(theta, mu=mu_list[[j]], sigma=sigma_list[[j]])
  }
  mix_dens <- rowSums(comp_vals)
  log(mix_dens)
}

density_grid <- function(n_grid = 100, xr = c(-5,5), yr = c(-10,10),
                         xy_grid = NULL, cell_dim = NULL, z_raw = NULL,
                         log_posterior = NULL, use_mass = F, theta = NULL, 
                         dat = NULL) {
  
  if(is.null(xy_grid)){
    xy_grid <- expand.grid(x = seq(xr[1],xr[2],length.out=n_grid+1), 
                           y = seq(yr[1],yr[2],length.out=n_grid+1))
  }
  if(is.null(cell_dim)){
    cell_dim <- list(w = diff(xr) / (n_grid), 
                     h = diff(yr) / (n_grid))  
  }
  if(is.null(z_raw)){
    z_raw <- z <- log_posterior(xy_grid, dat)  
  } else {
    z <- z_raw
  }
  
  #normalize z_raw to total cell count
  z_raw <- z_raw - log(sum(exp(z_raw)))
  if(!use_mass){
    #transform to an average density otherwise
    z_raw <- z_raw - log(prod(unlist(cell_dim)))
  }
  
  #for easier visual distinction
  z <- z - log(mean(exp(z)))
  z[z<0] <- -Inf
  
  #transform to colors
  z <- exp(z)
  z <- z / max(z) * 100
  xy_sub <- xy_grid[z>0,]
  z_sub <- z[z>0]
  zraw_sub <- z_raw[z>0]
  colpal <- colorRampPalette(c("white", "black"))(100)
  sqrtz <- sqrt(z_sub)
  sqrtz <- sqrtz - min(sqrtz)
  sqrtz <- sqrtz / max(sqrtz) * 99 + 1
  zcoli <- round(sqrtz)
  cols <- colpal[zcoli]
  return(list(
    xy_sub = xy_sub,
    cell_dim = cell_dim,
    cols = cols,
    colpal = colpal,
    zraw_sub = zraw_sub,
    sqrtz = sqrtz,
    xy_grid = xy_grid
  ))
}

compute_log_hist <- function(theta_final, n_grid, xr, yr, use_mass = F) {
  # 1) Some setup
  n <- nrow(theta_final)
  counts <- numeric(n_grid^2)  # 1D array of counts for each cell
  
  # 2) Grid spacing for x, y
  #    We'll treat [xr[1], xr[2]) as subdivided into n_grid uniform bins.
  #    So each bin is step_x wide in x, step_y wide in y.
  step_x <- (xr[2] - xr[1]) / n_grid
  step_y <- (yr[2] - yr[1]) / n_grid
  
  # 3) Bin each point in O(n) time
  x_data <- theta_final[,1]
  y_data <- theta_final[,2]
  
  # for i=1..n
  for (i in seq_len(n)) {
    xi <- x_data[i]
    yi <- y_data[i]
    
    # 3a) Compute bin indices along x, y
    #     i_x in [0, n_grid-1] if within range
    i_x <- floor((xi - xr[1]) / step_x)
    i_y <- floor((yi - yr[1]) / step_y)
    
    # 3b) Check boundary
    if (i_x < 0 || i_x >= n_grid) next
    if (i_y < 0 || i_y >= n_grid) next
    
    # 3c) Flatten 2D -> 1D index
    # by convention, let bin_id = i_y*n_grid + i_x + 1
    bin_id <- i_y * n_grid + i_x + 1  # +1 because R is 1-based
    counts[bin_id] <- counts[bin_id] + 1
  }
  
  # 4) Convert to log-frequency
  total_count <- sum(counts)
  if (total_count == 0) {
    # if all points were out of [xr, yr], everything is -Inf
    return(rep(-Inf, n_grid^2))
  }
  freq <- counts / total_count
  
  if(use_mass){
    log_freq <- log(freq)  # zero-count cells => -Inf
    return(log_freq)  
  } else {
    freq / (step_x * step_x)
  }
  
}

#### initial fit ####

#define log posterior
dat <- list(
  blob_1_mu = c(-1.5, 4),
  blob_1_sd = c(0.5, 0.5),
  blob_2_mu = c(1.5, 4),
  blob_2_sd = c(0.5, 0.5),
  banana_mu = c(0, 0),
  banana_sd = c(1, 1),
  probs = c(0.75, 0.125, 0.125)
)

#specify grid for plotting and plot target distribution
grid_info <- density_grid(log_posterior = log_posterior, dat = dat)
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
title(main = "Evil Smiley Face Distribution")
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, main = "target\nlogdens")


#fit normal approx with intiial importance resampling
n_draws <- 1000
num_paths <- 50
fit_pf <- mod$pathfinder(
  data = dat,
  draws = n_draws,
  num_paths = num_paths,
  psis_resample = TRUE,
  init = 8
)

#extract samples
pf_draws_all <- fit_pf$draws()  # e.g., Nx2 matrix for (x,y)
pf_draws <- as.data.frame(pf_draws_all[,c("x", "y")])  # e.g., Nx2 matrix for (x,y)

#visualize against true distribution
par(mar = c(4,4,2,6))
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
points(pf_draws, col = adjustcolor(2, 0.2))
title(main = paste0("Pathfinder Output, Number of Paths = ", num_paths))
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, main = "target\nlogdens")
legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, legend = c("pathfinder\nsamples"), pch = 19, 
       col = c(adjustcolor(2, 0.2), 4), bty = "n", xpd = NA)

#visualize the implied density
n_draws_lots <- 1E6
fit_pf_lots <- mod$pathfinder(
  data = dat,
  draws = n_draws_lots,
  num_paths = num_paths,
  psis_resample = TRUE,
  init = 8
)
pf_draws_lots <- as.data.frame(fit_pf_lots$draws(c("x", "y")))
# pf_draws_lots <- pf_draws_lots[!duplicated(pf_draws_lots),]
# pfdl_dupe <- pf_draws_lots[duplicated(pf_draws_lots),]
pf_dens <- log(compute_log_hist(theta_final = pf_draws_lots,
                                      n_grid = sqrt(nrow(grid_info$xy_grid)),
                                      xr = range(grid_info$xy_grid$x),
                                      yr = range(grid_info$xy_grid$y)))
pf_dens_grid_info <- density_grid(cell_dim = grid_info$cell_dim, xy_grid = grid_info$xy_grid,
                                   z_raw = pf_dens)
plot_density_grid(pf_dens_grid_info$xy_sub, pf_dens_grid_info$cell_dim, pf_dens_grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
title(main = paste0("Pathfinder Output, Number of Paths = ", num_paths))
add_continuous_legend(colors = rev(pf_dens_grid_info$colpal), 
                      labels = pf_dens_grid_info$zraw_sub, 
                      positions = pf_dens_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "est.\nlogdens")

#fit with MCMC
n_draws <- 1E5
fit_mcmc <- mod$sample(
  data = dat, 
  chains = 4, iter_warmup = n_draws/2, iter_sampling = n_draws/4, parallel_chains = 4, 
  refresh = 1E4
)
mcmc_draws <- as.data.frame(as_draws_df(fit_mcmc$draws()))[,c("x", "y")]
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
points(mcmc_draws[1:1E3,], col = adjustcolor(2, 0.2))
title(main = paste0("MCMC Output"))
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, main = "target\nlogdens")
legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, legend = c("MCMC\nsamples"), pch = 19, 
       col = c(adjustcolor(2, 0.2), 4), bty = "n", xpd = NA)

#plot estimated MCMC dens
mcmc_dens <- log(compute_log_hist(theta_final = mcmc_draws,
                                n_grid = sqrt(nrow(grid_info$xy_grid)),
                                xr = range(grid_info$xy_grid$x),
                                yr = range(grid_info$xy_grid$y)))
mcmc_dens_grid_info <- density_grid(cell_dim = grid_info$cell_dim, xy_grid = grid_info$xy_grid,
                                  z_raw = mcmc_dens)
plot_density_grid(mcmc_dens_grid_info$xy_sub, mcmc_dens_grid_info$cell_dim, mcmc_dens_grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
title(main = paste0("MCMC Output"))
add_continuous_legend(colors = rev(mcmc_dens_grid_info$colpal), 
                      labels = mcmc_dens_grid_info$zraw_sub, 
                      positions = mcmc_dens_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "estimated\nlogdens")

#### mixture approx ####
kmeans_for_means <- T
best_k  <- NULL
best_bic <- Inf
best_fit <- NULL
max_k   <- 6
theta_current <- pf_draws

for (k in 2:max_k) {
  
  #specify initial values
  if(kmeans_for_means){
    kmeans_out <- kmeans(theta_current, centers = k)
    kmeans_clusts <- split(1:n_draws, kmeans_out$cluster)
    init_lists <- list(lambda = kmeans_out$size / sum(kmeans_out$size),
                       mu = apply(kmeans(theta_current, centers = k)$centers, 1, identity, simplify = F),
                       sigma = lapply(kmeans_clusts, function(ci) cov(theta_current[ci,])))
  } else {
    if(k == 2){
      init_lists <- list(lambda = NULL,
                           mu = NULL,
                           sigma = NULL)
    } else {
      init_lists <- birth_component(old_fit = best_fit, data = theta_current, 
                                    new_weight = 0.05)
    }
  }
  
  #fit mixture model
  mixfit <- mvnormalmixEM(
    x      = theta_current,
    k      = k,
    lambda = init_lists$lambda,
    mu     = init_lists$mu,
    sigma  = init_lists$sigma,
    maxit   = ifelse(kmeans_for_means, 1, 200),
    epsilon = 1e-2,
    verb   = T
  )

  bic_val <- compute_bic(mixfit, theta_current)
  if (bic_val < best_bic) {
    best_bic <- bic_val
    best_fit <- mixfit
  }
}

#plot samples
resampled_draws <- sample_from_mix(n = n_draws, gm_fit = best_fit)
lp_true <- log_posterior(resampled_draws, blob_x, blob_y, probs)
lp_prop <- log_mix_density(theta = resampled_draws, 
                           mixfit = best_fit)
log_weights <- lp_true - lp_prop

#PSIS
psis_fit <- psis(log_weights)
summary(psis_fit$diagnostics$pareto_k)
w_smooth <- weights(psis_fit, normalize=TRUE)
idx  <- sample.int(n=nrow(w_smooth), size=nrow(w_smooth), 
                   replace=TRUE, prob=exp(w_smooth))
reresampled_draws <- resampled_draws[idx, ]

#plot
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
points(reresampled_draws, col = adjustcolor(2, 0.2))
points(do.call(rbind, best_fit$mu), col = 4)
title(paste0("MVN Mixture Output (iter = ", 0, ", best k = ", length(best_fit$lambda), ")"))
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, main = "target\nlogdens")
legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, legend = c("samples", "means"), pch = 19, 
       col = c(adjustcolor(2, 0.2), 4), bty = "n", xpd = NA)


#### iterate ####
#iterate over this procedure a few timed and see if it converges
max_iter <- 5
max_k   <- 8
theta_current <- reresampled_draws

for (iter in seq_len(max_iter)) {
  cat("=== Iteration", iter, "===\n")
  
  # 1) Fit mixture with a range of k, pick best by BIC
  best_bic <- Inf
  for (k in 2:max_k) {
    
    #jitter theta to prevent singular sample covariances
    sample_sds <- apply(theta_current, 2, sd)
    theta_current[,1] <- theta_current[,1] + rnorm(n_draws) * sample_sds[1]/1E3
    theta_current[,2] <- theta_current[,2] + rnorm(n_draws) * sample_sds[2]/1E3
    
    #specify initial values
    if(kmeans_for_means){
      kmeans_out <- kmeans(theta_current, centers = k)
      kmeans_clusts <- split(1:n_draws, kmeans_out$cluster)
      init_lists <- list(lambda = kmeans_out$size / sum(kmeans_out$size),
                         mu = apply(kmeans(theta_current, centers = k)$centers, 1, identity, simplify = F),
                         sigma = lapply(kmeans_clusts, function(ci) cov(theta_current[ci,])))
    } else {
      if(k == 2){
        init_lists <- list(lambda = NULL,
                           mu = NULL,
                           sigma = NULL)
      } else {
        init_lists <- birth_component(old_fit = best_fit, data = theta_current, 
                                      new_weight = 0.05)
      }
    }
    
    #fit mixture model
    mixfit <- mvnormalmixEM(
      x      = theta_current,
      k      = k,
      lambda = init_lists$lambda,
      mu     = init_lists$mu,
      sigma  = init_lists$sigma,
      maxit   = ifelse(kmeans_for_means, 1, 200),
      epsilon = 1e-2,
      verb   = T
    )
    
    #snag BIC
    bic_val <- compute_bic(mixfit, theta_current)
    if (bic_val < best_bic) {
      best_bic <- bic_val
      best_fit <- mixfit
    }
  }
  cat("  best k:", length(best_fit$lambda), ", BIC=", best_bic, "\n")
  
  #draw new samples from mixture
  new_draws <- sample_from_mix(n_draws, best_fit)
  
  #calculate importance weights
  lp_true <- log_posterior(new_draws, blob_x, blob_y, probs)
  lp_prop <- log_mix_density(new_draws, best_fit)
  log_w   <- lp_true - lp_prop
  
  #smooth weights
  psis_out <- psis(log_w)
  w_smooth <- weights(psis_out, normalize=TRUE)
  cat("  mean pareto-k:", mean(psis_out$diagnostics$pareto_k), "\n")
  
  #resample from posterior
  idx  <- sample.int(n=n_draws, size=n_draws, replace=TRUE, prob=exp(w_smooth))
  theta_current <- new_draws[idx, ]
  
  #plot the result
  plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                    xlim = c(-3,3), ylim = c(-3,7))
  points(theta_current, col = adjustcolor(2, 0.2))
  points(do.call(rbind, best_fit$mu), col = 4)
  title(paste0("MVN Mixture Output (iter = ", iter, ", best k = ", length(best_fit$lambda), ")"))
  add_continuous_legend(colors = rev(grid_info$colpal), 
                        labels = grid_info$zraw_sub, 
                        positions = grid_info$sqrtz/100, left_below = T, 
                        w = diff(par("usr")[1:2])/30, 
                        x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                        y = par("usr")[4] - diff(par("usr")[3:4])/20, main = "target\nlogdens")
  legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, legend = c("samples", "means"), pch = 19, 
         col = c(adjustcolor(2, 0.2), 4), bty = "n", xpd = NA)
  
}

# final logdens for mixture
logmixdens <- log_mix_density(theta = as.matrix(grid_info$xy_grid), mixfit = best_fit)
final_grid_info <- density_grid(cell_dim = grid_info$cell_dim, xy_grid = grid_info$xy_grid, z_raw = logmixdens)
plot_density_grid(final_grid_info$xy_sub, final_grid_info$cell_dim, final_grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
title(paste0("Final Mixture (iter = ", iter, ", k = ", length(best_fit$lambda), ")"))
add_continuous_legend(colors = rev(final_grid_info$colpal), 
                      labels = final_grid_info$zraw_sub, 
                      positions = final_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "estimated\nlogdens")
points(do.call(rbind, best_fit$mu), col = 4)
legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, 
       legend = c("means"), pch = 19, 
       col = c(4), bty = "n", xpd = NA)

#final logens for importance-weighted samples
n_draws_final <- 1E5
new_draws <- sample_from_mix(n_draws_final, best_fit)
lp_true <- log_posterior(new_draws, blob_x, blob_y, probs)
lp_prop <- log_mix_density(new_draws, best_fit)
log_w   <- lp_true - lp_prop
psis_out <- psis(log_w)
w_smooth <- weights(psis_out, normalize=TRUE)
idx  <- sample.int(n=n_draws_final, size=n_draws_final, replace=TRUE, prob=exp(w_smooth))
theta_final <- new_draws[idx, ]
iw_logmixdens <- log(compute_log_hist(theta_final = theta_final, 
                                  n_grid = sqrt(nrow(grid_info$xy_grid)), 
                                  xr = range(grid_info$xy_grid$x), 
                                  yr = range(grid_info$xy_grid$y)))
iw_final_grid_info <- density_grid(cell_dim = grid_info$cell_dim, xy_grid = grid_info$xy_grid, 
                                   z_raw = iw_logmixdens)

plot_density_grid(iw_final_grid_info$xy_sub, iw_final_grid_info$cell_dim, iw_final_grid_info$cols,
                  xlim = c(-3,3), ylim = c(-3,7))
title(paste0("Final Mixture (Importance Weighted)\n(iter = ", iter, ", k = ", 
             length(best_fit$lambda), ")"), line = -0.5)
add_continuous_legend(colors = rev(iw_final_grid_info$colpal), 
                      labels = iw_final_grid_info$zraw_sub, 
                      positions = iw_final_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/20, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "estimated\nlogdens")
points(do.call(rbind, best_fit$mu), col = 4)
legend(x = par("usr")[2], y = par("usr")[3] + diff(par("usr")[3:4])/6, 
       legend = c("means"), pch = 19, 
       col = c(4), bty = "n", xpd = NA)


