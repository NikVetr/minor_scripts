#### libraries ####
library(loo)
library(mixtools)
library(MASS)
library(cmdstanr)
library(posterior)

# Write code to a file "banana.stan"
model_path <- "/Volumes/4TB/stan_projects/evil-smiley-bug/evil_smiley.stan"
mod <- cmdstan_model(model_path)

#### functions #### 
functions_path <- "/Volumes/4TB/stan_projects/evil-smiley-bug/evil_smiley_functions.R"
source(functions_path)

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

#specify grid for plotting and plot the target distribution
par(mar = c(4,4,2,6))
xlim <- c(-3,3)
ylim <- c(-3,7)
grid_info <- density_grid(log_posterior = log_posterior, dat = dat)
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = xlim, ylim = ylim)
title(main = "Evil Smiley Face Distribution")
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/10, 
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

#visualize samples against true distribution
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = xlim, ylim = ylim)
points(pf_draws, col = adjustcolor(2, 0.2))
title(main = paste0("Pathfinder Output, Number of Paths = ", num_paths))
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/10, 
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
                  xlim = xlim, ylim = ylim)
title(main = paste0("Pathfinder Output, Number of Paths = ", num_paths))
add_continuous_legend(colors = rev(pf_dens_grid_info$colpal), 
                      labels = pf_dens_grid_info$zraw_sub, 
                      positions = pf_dens_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/10, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "estimated\nlogdens")

#fit with MCMC
n_draws <- 1E5
fit_mcmc <- mod$sample(
  data = dat, 
  chains = 4, iter_warmup = n_draws/2, iter_sampling = n_draws/4, parallel_chains = 4, 
  refresh = 1E4
)
mcmc_draws <- as.data.frame(as_draws_df(fit_mcmc$draws()))[,c("x", "y")]
plot_density_grid(grid_info$xy_sub, grid_info$cell_dim, grid_info$cols,
                  xlim = xlim, ylim = ylim)
points(mcmc_draws[1:1E3,], col = adjustcolor(2, 0.2))
title(main = paste0("MCMC Output"))
add_continuous_legend(colors = rev(grid_info$colpal), 
                      labels = grid_info$zraw_sub, 
                      positions = grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/10, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "target\nlogdens")
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
                  xlim = xlim, ylim = ylim)
title(main = paste0("MCMC Output"))
add_continuous_legend(colors = rev(mcmc_dens_grid_info$colpal), 
                      labels = mcmc_dens_grid_info$zraw_sub, 
                      positions = mcmc_dens_grid_info$sqrtz/100, left_below = T, 
                      w = diff(par("usr")[1:2])/30, 
                      x = par("usr")[2] + diff(par("usr")[1:2])/10, 
                      y = par("usr")[4] - diff(par("usr")[3:4])/20, 
                      main = "estimated\nlogdens")

#also look at marginal pf histograms
par(mfrow = c(2,1), mar = c(4,4,2,2))
hist(pf_draws_lots$x, breaks = 100)
hist(pf_draws_lots$y, breaks = 100)
