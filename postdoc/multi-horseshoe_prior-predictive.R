# minimal prior predictive with cmdstanr
library(cmdstanr)
library(posterior)

# user options
use_pf <- FALSE  # set TRUE to try pathfinder for prior draws

# build data for prior-only run
n <- nrow(X_train)
d <- ncol(X_train)

# prepare inputs
blk_names <- names(p_vec)

# Heuristic default for expected nonzeros per block (edit if you have priors):
p0_guess <- c(phenotypes = 5, molecular_analyte = 50, molecular_kPCs = 20)
p0_guess <- p0_guess[match(blk_names, names(p0_guess))]  # or write a small align helper
p0_guess[is.na(p0_guess)] <- pmax(1L, round(0.05 * p_vec[is.na(p0_guess)]))
names(p0_guess) <- blk_names

# slab scales: keep your old “phenotypes=4, others=2” by name heuristic, else 2
slab_scale <- if (length(blk_names)) {
  sapply(blk_names, function(nm) if (grepl("pheno", nm, ignore.case = TRUE)) 4 else 2)
} else {
  rep(2, G)
}
slab_df <- 5

# compute tau0 by block
tau0_block <- function(p0, d, n, sigma = 1) (p0 / (d - p0)) * (sigma / sqrt(n))
hs_tau0 <- mapply(function(p0, d) tau0_block(p0, d, n),
                  p0 = p0_guess, d = p_vec, SIMPLIFY = TRUE)
hs <- list(
  hs_tau0    = hs_tau0,
  slab_scale = slab_scale,
  slab_df    = slab_df
)


dat <- list(
  n = n,
  d = d,
  x = as.matrix(X_train),
  y = numeric(n),                  # ignored when prior_only = 1
  G = as.integer(G),
  group = as.integer(block_id),
  hs_tau0 = as.vector(hs_tau0),
  slab_scale = as.vector(slab_scale),
  slab_df = slab_df,
  prior_only = 1L
)

# compile model (the tweaked file with 'prior_only' and 'y_prior')
hs_mod <- cmdstan_model("~/scripts/minor_scripts/postdoc/multi-horseshoe_pp.stan")

# draw from priors (no likelihood)
if (use_pf) {
  message("[debug] using pathfinder for prior draws")
  fit <- hs_mod$pathfinder(data = dat, draws = 2000, num_paths = 8, refresh = 0)
} else {
  message("[debug] using NUTS for prior draws")
  fit <- hs_mod$sample(
    data = dat,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    thin = 1,
    adapt_delta = 0.95,
    max_treedepth = 12,
    parallel_chains = 4,
    refresh = 100
  )
}

# extract prior predictive for y
# observation-level prior predictive (n columns)
y_mat <- data.frame(posterior::as_draws_matrix(fit$draws("y_prior")))

# predictor-level parameters (d columns)
beta_mat  <- data.frame(posterior::as_draws_matrix(fit$draws("beta")))
kappa_mat <- data.frame(posterior::as_draws_matrix(fit$draws("kappa")))

idx_by_group <- split(seq_len(ncol(beta_mat)), block_id)
beta_by_group <- lapply(idx_by_group, function(j) unlist(beta_mat[, j, drop = FALSE]))
hist(beta_by_group[[1]], breaks = 1000, xlim = c(-5,5))
hist(beta_by_group[[2]], breaks = 1000, xlim = c(-5,5))

kappa_by_group <- lapply(idx_by_group, function(j) unlist(kappa_mat[, j, drop = FALSE]))
h <- hist(kappa_by_group[[2]], breaks = 100, plot = FALSE)
# plot with log-scaled density axis
plot(h$breaks[-1] - diff(h$breaks),
     log10(h$density),
     xlab = "kappa",
     ylab = "log10(density)",
     main = "shrinkage factor prior", type = "l", yaxt = "n")
yax_vals <- pretty(log10(h$density))
axis(2, at = yax_vals, labels = rep("", length(yax_vals)))
text(par("usr")[1], yax_vals, labels = round(10^yax_vals, 2), pos = 2, xpd = NA)
abline(h = log10(1), lty = 2, col = adjustcolor(1, 0.5))
