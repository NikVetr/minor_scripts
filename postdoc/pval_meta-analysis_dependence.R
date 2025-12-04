source("~/scripts/minor_scripts/postdoc/scatter_hist.R")

#### functions ####

binom_ci <- function(k, n, conf = 0.90) {
  alpha <- 1 - conf
  lower <- qbeta(alpha / 2, k, n - k + 1)
  upper <- qbeta(1 - alpha / 2, k + 1, n - k)
  return(c(lower = lower, upper = upper))
}

# combination methods

# fisher's method
combine_fisher <- function(p) {
  # debugging: check for invalid p-values
  if (any(is.na(p))) {
    stop("combine_fisher: NA values in p")
  }
  m <- length(p)
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  stat <- -2 * sum(log(p))
  combined_p <- pchisq(stat, df = 2 * m, lower.tail = FALSE)
  return(combined_p)
}

# harmonic mean p-value
# calibrated harmonic mean p-value (HMP, Wilson 2019)
combine_harmonic <- function(p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  p_comb <- harmonicmeanp::p.hmp(p, L = length(p), multilevel = FALSE)
  if (!is.finite(p_comb)) {
    p_comb <- eps
  }
  return(p_comb)
}

# cauchy combination test (Liu & Xie)
combine_cauchy <- function(p, weights = NULL) {
  m <- length(p)
  if (is.null(weights)) {
    weights <- rep(1 / m, m)
  } else {
    weights <- weights / sum(weights)
  }
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  t_vals <- tan((0.5 - p) * pi)
  t_stat <- sum(weights * t_vals)
  combined_p <- 0.5 - atan(t_stat) / pi
  return(combined_p)
}

combine_tippett <- function(p) {
  m <- length(p)
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  p_min <- min(p)
  # under independence: min(U_1,...,U_m) has CDF 1 - (1 - x)^m
  combined_p <- 1 - (1 - p_min)^m
  # equivalently: combined_p <- pbeta(p_min, shape1 = 1, shape2 = m, lower.tail = TRUE)
  return(combined_p)
}

combine_stouffer <- function(p, weights = NULL) {
  
  eps <- 1e-15
  m <- length(p)
  p <- pmin(pmax(p, eps), 1 - eps)
  z <- qnorm(1 - p)
  
  if (is.null(weights)) {
    z_sum <- sum(z)
    z_comb <- z_sum / sqrt(m)
  } else {
    if (length(weights) != m) {
      stop("combine_stouffer: length(weights) must equal length(p)")
    }
    # standard weighted Stouffer: Z = sum(w_i z_i) / sqrt(sum w_i^2)
    w <- as.numeric(weights)
    z_comb <- sum(w * z) / sqrt(sum(w^2))
  }
  
  # combined p is upper-tail
  p_comb <- pnorm(z_comb, lower.tail = FALSE)
  p_comb <- pmin(pmax(p_comb, eps), 1)
  
  return(p_comb)
}


combine_simes <- function(p) {
  m <- length(p)
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  
  # Simes: sort p_(1) <= ... <= p_(m), then use min_i m * p_(i) / i
  p_sorted <- sort(p, decreasing = FALSE)
  i <- seq_len(m)
  adj <- m * p_sorted / i
  p_simes <- min(adj)
  
  # keep in [eps, 1]
  p_simes <- pmin(pmax(p_simes, eps), 1)
  
  return(p_simes)
}

combine_brownCS <- function(p,
                            rho = NULL,
                            estimate_rho = TRUE) {
  # Brown/Kost-style correction to Fisher's method
  # under compound symmetry (single correlation rho).
  # If rho is NULL and estimate_rho = TRUE, we estimate
  # rho via estimate_gaussian_struct_cs(p)$coef$rho.
  
  if (any(is.na(p))) {
    stop("combine_brownCS: NA values in p")
  }
  
  m <- length(p)
  if (m == 0L) {
    stop("combine_brownCS: length(p) must be > 0")
  }
  if (m == 1L) {
    # with a single p-value, Brown reduces to identity
    return(pmin(pmax(p[1], 0), 1))
  }
  
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  
  # Fisher statistic
  X <- -2 * log(p)
  S <- sum(X)
  
  # estimate or use supplied CS correlation for latent normals
  if (is.null(rho)) {
    if (!estimate_rho) {
      stop("combine_brownCS: rho is NULL and estimate_rho = FALSE")
    }
    # your estimator expects values in (0,1).
    # feeding p or (1 - p) are equivalent for correlation;
    # both give the same rho.
    est <- estimate_gaussian_struct_cs(p)
    rho <- est$coef$rho
  }
  
  # clamp rho to a reasonable range to avoid pathological cases
  # (CS requires rho > -1/(m-1); estimator already respects this,
  # but we protect against numerical drift)
  rho_min <- -1 / (m - 1) + 1e-6
  rho_max <- 0.999
  rho <- max(min(rho, rho_max), rho_min)
  
  # Kost & McDermott approximation:
  # cov(-2 log p_i, -2 log p_j) ≈ 3.263 * rho + 0.710 * rho^2 + 0.027 * rho^3
  # for all i != j
  cov_x <- 3.263 * rho + 0.710 * rho^2 + 0.027 * rho^3
  
  # Brown's moment-matching:
  # E[S] = 2m
  # Var[S] = 4m + 2 * sum_{i<j} cov(-2 log p_i, -2 log p_j)
  #       = 4m + m(m - 1) * cov_x   under CS
  E_S   <- 2 * m
  Var_S <- 4 * m + m * (m - 1) * cov_x
  
  if (!is.finite(Var_S) || Var_S <= 0) {
    warning("combine_brownCS: Non-positive or non-finite Var_S; falling back to Fisher.")
    combined_p <- pchisq(S, df = 2 * m, lower.tail = FALSE)
    combined_p <- pmin(pmax(combined_p, eps), 1)
    return(combined_p)
  }
  
  # scaled-chi-square approximation:
  # S ~ c * chi^2_{v}, with
  #   c = Var_S / (2 * E_S)
  #   v = 2 * E_S^2 / Var_S
  c_scale <- Var_S / (2 * E_S)
  v_df    <- 2 * E_S^2 / Var_S
  
  # guard against tiny df or negative scaling
  if (!is.finite(c_scale) || c_scale <= 0 || !is.finite(v_df) || v_df <= 0) {
    warning("combine_brownCS: Invalid (c_scale, v_df); falling back to Fisher.")
    combined_p <- pchisq(S, df = 2 * m, lower.tail = FALSE)
    combined_p <- pmin(pmax(combined_p, eps), 1)
    return(combined_p)
  }
  
  # combined p-value
  combined_p <- pchisq(S / c_scale, df = v_df, lower.tail = FALSE)
  combined_p <- pmin(pmax(combined_p, eps), 1)
  
  return(combined_p)
}

combine_zcor_stouffer <- function(z, cor_mat = NULL, weights = NULL) {
  # z: vector of z-statistics (e.g. signed z from two-sided p)
  # cor_mat: estimated correlation matrix of z under the null (from estimate_null_covariance)
  #         if NULL, assumes independence (identity)
  # weights: optional weights; default = equal weights
  
  if (any(is.na(z))) {
    stop("combine_stouffer_z_cor: NA values in z")
  }
  
  m <- length(z)
  
  if (is.null(weights)) {
    w <- rep(1, m)
  } else {
    if (length(weights) != m) {
      stop("combine_stouffer_z_cor: length(weights) must equal length(z)")
    }
    w <- as.numeric(weights)
  }
  
  if (is.null(cor_mat)) {
    # independence case: var = sum w_i^2
    var_T <- sum(w^2)
  } else {
    if (!is.matrix(cor_mat) || nrow(cor_mat) != m || ncol(cor_mat) != m) {
      stop("combine_stouffer_z_cor: cor_mat must be an m x m matrix")
    }
    # ensure symmetry
    cor_mat <- (cor_mat + t(cor_mat)) / 2
    # variance of T = w' R w, assuming unit variances
    var_T <- as.numeric(t(w) %*% cor_mat %*% w)
  }
  
  if (!is.finite(var_T) || var_T <= 0) {
    warning("combine_stouffer_z_cor: non-positive or non-finite variance; returning NA")
    return(NA_real_)
  }
  
  T <- sum(w * z)
  z_comb <- T / sqrt(var_T)
  p_comb <- pnorm(z_comb, lower.tail = FALSE)
  
  p_comb
}

combine_stouffer_perm <- function(z_obs, stats_mat_perm) {
  # z_obs: length-p vector of observed z-statistics
  # stats_mat_perm: p x B matrix of z-statistics under null (each col = one permutation)
  
  if (any(is.na(z_obs))) stop("combine_stouffer_perm: NA in z_obs")
  if (!is.matrix(stats_mat_perm)) stop("combine_stouffer_perm: stats_mat_perm must be a matrix")
  
  p <- length(z_obs)
  if (nrow(stats_mat_perm) != p) stop("combine_stouffer_perm: nrow(stats_mat_perm) must equal length(z_obs)")
  
  # global Stouffer-like statistic: sum of z
  T_obs <- sum(z_obs)
  T_perm <- colSums(stats_mat_perm)
  
  B <- ncol(stats_mat_perm)
  # upper-tail p-value
  p_perm <- (1 + sum(T_perm >= T_obs)) / (B + 1)
  
  p_perm
}


estimate_gaussian_struct_cs <- function(u,
                                        estimate_mean = FALSE,
                                        estimate_variances = FALSE,
                                        verbose = FALSE) {
  # ensure u is a plain numeric vector
  if (!is.numeric(u)) {
    stop("u must be numeric")
  }
  u <- as.numeric(u)  # drop matrix/array attributes
  
  if (any(u <= 0 | u >= 1)) {
    stop("u must be strictly between 0 and 1")
  }
  
  # clip slightly away from 0 and 1 for numerical stability
  eps <- .Machine$double.eps
  u_clipped <- pmin(pmax(u, eps), 1 - eps)
  
  # transform to latent normal
  z <- qnorm(u_clipped)
  p <- length(z)
  
  if (p < 2) {
    stop("need at least 2 dimensions to estimate correlation structure")
  }
  
  if (verbose) {
    cat("dimension p:", p, "\n")
  }
  
  # choose mean of latent normal
  if (estimate_mean) {
    mu_scalar <- mean(z)
  } else {
    mu_scalar <- 0
  }
  
  # center z
  x <- z - mu_scalar
  
  # choose variance of latent normal
  if (estimate_variances) {
    sigma2 <- mean(x^2)
    if (sigma2 <= 0) {
      stop("estimated variance is non-positive; something is wrong")
    }
  } else {
    sigma2 <- 1
  }
  
  if (verbose) {
    cat("mu_scalar:", mu_scalar, "sigma2:", sigma2, "\n")
  }
  
  # precompute simple sufficient statistics
  s <- sum(x)        # sum x_i
  t <- sum(x^2)      # sum x_i^2 = ||x||^2
  
  # map unconstrained eta_rho to valid rho in (lower, upper)
  rho_from_eta <- function(eta, p) {
    lower <- -1 / (p - 1) + 1e-6
    upper <- 0.999
    t_ <- 1 / (1 + exp(-eta))  # logistic
    rho <- lower + (upper - lower) * t_
    rho
  }
  
  # derivative drho/deta for delta-method SE
  drho_deta <- function(eta, p) {
    lower <- -1 / (p - 1) + 1e-6
    upper <- 0.999
    t_ <- 1 / (1 + exp(-eta))
    dt_deta <- t_ * (1 - t_)
    (upper - lower) * dt_deta
  }
  
  # negative log-likelihood as a function of eta (scalar)
  neg_loglik <- function(eta) {
    rho <- rho_from_eta(eta, p)
    
    # compound symmetry correlation matrix:
    # R = (1 - rho) I + rho J
    a <- 1 - rho
    b <- rho
    
    # log |R| = (p-1) log(a) + log(a + b p)
    # a + b p = 1 - rho + rho p
    apb <- a + b * p
    
    if (a <= 0 || apb <= 0) {
      if (verbose) {
        cat("non positive-definite parameters for rho =", rho, "\n")
      }
      return(1e12)
    }
    
    logdetR <- (p - 1) * log(a) + log(apb)
    
    # R^{-1} = (1/a) I - (b / [a (a + b p)]) J
    # quadratic form x' R^{-1} x = (1/a) * t - (b / [a (a + b p)]) * s^2
    quad <- (1 / a) * t - (b / (a * apb)) * s^2
    
    # log-likelihood for X ~ N(mu 1, sigma2 R)
    # loglik = -0.5 * [p log(2*pi) + p log(sigma2) + log|R| + quad / sigma2]
    loglik <- -0.5 * (p * log(2 * pi) + p * log(sigma2) + logdetR + quad / sigma2)
    
    if (verbose) {
      cat("eta:", eta, "rho:", rho, "loglik:", loglik, "\n")
    }
    
    -loglik
  }
  
  # initial eta: 0 -> roughly rho ~ 0
  theta0 <- 0
  
  opt <- optim(par = theta0,
               fn = neg_loglik,
               method = "BFGS",
               hessian = TRUE,
               control = list(maxit = 1000,
                              trace = if (verbose) 1 else 0))
  
  eta_hat <- opt$par[1]
  rho_hat <- rho_from_eta(eta_hat, p)
  
  # log-likelihood at optimum
  loglik_hat <- -neg_loglik(eta_hat)
  
  # number of free parameters: only rho is optimized
  k <- 1
  AIC_val <- 2 * k - 2 * loglik_hat
  
  # approximate SE for rho via Hessian (1D)
  se_rho <- NA
  if (!is.null(opt$hessian)) {
    h11 <- opt$hessian[1, 1]
    if (is.finite(h11) && h11 > 0) {
      var_eta <- 1 / h11
      d_rho <- drho_deta(eta_hat, p)
      se_rho <- sqrt(max(var_eta, 0)) * abs(d_rho)
    } else if (verbose) {
      cat("hessian non-positive or non-finite; se_rho not available\n")
    }
  }
  
  # we are treating mu and sigma2 as plug-in, not free parameters
  se_mu <- 0
  se_sigma2 <- 0
  
  # return a minimal object (no huge matrices)
  list(
    struct = "cs",
    estimate_mean = estimate_mean,
    estimate_variances = estimate_variances,
    coef = list(
      rho = rho_hat,
      mean = mu_scalar,
      sigma2 = sigma2
    ),
    se = list(
      rho = se_rho,
      mean = se_mu,
      sigma2 = se_sigma2
    ),
    logLik = loglik_hat,
    AIC = AIC_val,
    converged = (opt$convergence == 0),
    message = opt$message,
    hessian = opt$hessian,
    # also return some handy summaries of x
    stats = list(
      n = p,
      sum_x = s,
      sum_x2 = t
    )
  )
}

# multivariate normal simulator
build_mu_vec <- function(m, mu_alt, prop_alt) {
  if (prop_alt < 0 || prop_alt > 1) {
    stop("build_mu_vec: prop_alt must be between 0 and 1")
  }
  
  mu_vec <- rep(0, m)
  num_alt <- round(m * prop_alt)
  num_alt <- max(0, min(m, num_alt))
  
  if (num_alt > 0) {
    # choose which indices are non-null; deterministic or random
    # deterministic choice:
    alt_idx <- seq_len(num_alt)
    # if you prefer random positions, comment the above and uncomment:
    # alt_idx <- sample.int(m, num_alt)
    
    mu_vec[alt_idx] <- mu_alt
  }
  
  return(mu_vec)
}

# simulate B draws of length-m normal vectors
# corr_structure: "indep", "cs", "ar1"
simulate_X <- function(B, m, mu_vec,
                       corr_structure = c("indep", "cs", "ar1"),
                       rho = 0) {
  corr_structure <- match.arg(corr_structure)
  if (length(mu_vec) != m) {
    stop("simulate_X: length(mu_vec) must equal m")
  }
  if (corr_structure == "indep") {
    # debugging
    # cat("simulate_X: using independent normals\n")
    X <- matrix(rnorm(B * m), nrow = B, ncol = m)
    X <- sweep(X, 2, mu_vec, "+")
    return(X)
  }
  
  # construct correlation matrix
  if (corr_structure == "cs") {
    # compound symmetry
    Sigma <- matrix(rho, nrow = m, ncol = m)
    diag(Sigma) <- 1
  } else if (corr_structure == "ar1") {
    # ar(1) correlation: rho^|i-j|
    idx <- 1:m
    Sigma <- outer(idx, idx, function(i, j) rho^abs(i - j))
  } else {
    stop("simulate_X: unknown corr_structure")
  }
  
  # cholesky; if not positive definite, this will error
  U <- chol(Sigma)
  
  Z <- matrix(rnorm(B * m), nrow = B, ncol = m)
  X <- Z %*% U
  X <- sweep(X, 2, mu_vec, "+")
  return(X)
}


# one simulation setting (one correlation structure + one mean vector)
run_single_setting <- function(B = 10000,
                               m = 20,
                               mu_alt = 0,
                               prop_alt = 1,
                               corr_structure = c("indep", "cs", "ar1"),
                               rho = 0,
                               alpha = 0.05,
                               comb_methods,
                               verbose = F) {
  corr_structure <- match.arg(corr_structure)
  
  mu_vec <- build_mu_vec(m = m, mu_alt = mu_alt, prop_alt = prop_alt)
  
  if(verbose){
    cat("run_single_setting:",
        "B =", B,
        "m =", m,
        "mu_alt =", mu_alt,
        "prop_alt =", prop_alt,
        "structure =", corr_structure,
        "rho =", rho, "\n")
    cat("  summary(mu_vec):\n")
    print(summary(mu_vec))
    cat("  number of nonzero means:", sum(mu_vec != 0), "\n")
  }
  
  X <- simulate_X(B, m, mu_vec,
                  corr_structure = corr_structure,
                  rho = rho)
  
  p_mat <- 1 - pnorm(X)
  
  # apply each combination method to each row of p_mat
  combined_list <- list()
  for (method_name in names(comb_methods)) {
    f <- comb_methods[[method_name]]
    combined_list[[method_name]] <- apply(p_mat, 1, f)
  }
  
  combined_p <- as.data.frame(combined_list)
  
  size_or_power <- vapply(
    names(comb_methods),
    function(method_name) {
      mean(combined_p[[method_name]] <= alpha)
    },
    numeric(1)
  )
  
  if(verbose){
    cat("  rejection rate (alpha =", alpha, "):\n")
    print(size_or_power)  
  }
  
  result <- list(
    settings = list(
      B = B,
      m = m,
      mu_alt = mu_alt,
      prop_alt = prop_alt,
      corr_structure = corr_structure,
      rho = rho,
      alpha = alpha,
      p_mat = p_mat
    ),
    rejection_rate = size_or_power,
    combined_p = combined_p
  )
  
  return(result)
}


# run the full study
run_full_study <- function(B = 10000,
                           m = 20,
                           mu_alt = 0.3,
                           prop_alt = 1,
                           rho_cs = 0.8,
                           rho_ar1 = 0.9,
                           alpha = 0.05,
                           comb_methods,
                           verbose = F) {
  structures <- c("indep", "cs", "ar1")
  results <- list()
  
  for (struct in structures) {
    if (struct == "indep") {
      rho <- 0
    } else if (struct == "cs") {
      rho <- rho_cs
    } else if (struct == "ar1") {
      rho <- rho_ar1
    } else {
      stop("run_full_study: unknown structure")
    }
    
    # null: all means zero => prop_alt is irrelevant, mu_alt unused
    if(verbose){
      cat("\n--- structure:", struct, "NULL (all mu = 0) ---\n")  
    }
    
    res_null <- run_single_setting(
      B = B,
      m = m,
      mu_alt = 0,
      prop_alt = 0,  # all coordinates truly null
      corr_structure = struct,
      rho = rho,
      alpha = alpha,
      comb_methods = comb_methods
    )
    
    # alternative: some proportion of coordinates have mu = mu_alt
    if(verbose){
      cat("\n--- structure:", struct,
          "ALT (mu_alt =", mu_alt,
          ", prop_alt =", prop_alt, ") ---\n")
    }
    res_alt <- run_single_setting(
      B = B,
      m = m,
      mu_alt = mu_alt,
      prop_alt = prop_alt,
      corr_structure = struct,
      rho = rho,
      alpha = alpha,
      comb_methods = comb_methods
    )
    
    results[[struct]] <- list(null = res_null, alt = res_alt)
  }
  
  return(results)
}


sweep_study <- function(base_vals,
                        range_vals,
                        comb_methods,
                        conf = 0.90,
                        verbose = TRUE) {
  
  sweep_param <- names(range_vals)[1]
  sweep_values <- range_vals[[1]]
  
  # unpack base values
  n_pval  <- base_vals$n_pval
  n_rep   <- base_vals$n_rep
  mu_alt  <- base_vals$mu_alt
  prop_alt <- base_vals$prop_alt
  rho_cs  <- base_vals$rho_cs
  rho_ar1 <- base_vals$rho_ar1
  alpha   <- base_vals$alpha
  
  depstructs <- c("indep", "cs", "ar1")
  methods <- names(comb_methods)
  conds <- c("null", "alt")
  
  # initialize result list: one data.frame per dependence structure
  sweep_results <- list()
  for (depstruct in depstructs) {
    df <- data.frame(value = sweep_values)
    names(df)[1] <- sweep_param
    
    for (method in methods) {
      for (cond in conds) {
        prefix <- paste(method, cond, sep = "_")
        df[[paste0(prefix, "_rejections")]] <- NA_integer_
        df[[paste0(prefix, "_rate")]]       <- NA_real_
        df[[paste0(prefix, "_low")]]        <- NA_real_
        df[[paste0(prefix, "_high")]]       <- NA_real_
      }
    }
    
    sweep_results[[depstruct]] <- df
  }
  
  # main loop over sweep_values
  for (i in seq_along(sweep_values)) {
    v <- sweep_values[i]
    
    # reset local copies of base values
    this_n_pval  <- n_pval
    this_mu_alt  <- mu_alt
    this_prop_alt <- prop_alt
    this_rho_cs  <- rho_cs
    this_rho_ar1 <- rho_ar1
    
    # apply the sweep for the chosen parameter
    if (sweep_param == "n_pval") {
      this_n_pval <- as.integer(v)
    } else if (sweep_param == "mu_alt") {
      this_mu_alt <- v
    } else if (sweep_param == "prop_alt") {
      this_prop_alt <- v
    } else if (sweep_param == "rho") {
      # use same rho for both cs and ar1, per your comment
      this_rho_cs  <- v
      this_rho_ar1 <- v
    }
    
    if (verbose) {
      cat(sprintf("[sweep_study] %s = %s (%d / %d)\n",
                  sweep_param, as.character(v), i, length(sweep_values)))
    }
    
    # run your existing full study
    res <- run_full_study(
      B = n_rep,
      m = this_n_pval,
      mu_alt = this_mu_alt,
      prop_alt = this_prop_alt,
      rho_cs = this_rho_cs,
      rho_ar1 = this_rho_ar1,
      alpha = alpha,
      comb_methods = comb_methods,
      verbose = FALSE
    )
    
    # for each dependence structure, compute rejections and CI
    for (depstruct in depstructs) {
      df <- sweep_results[[depstruct]]
      
      for (cond in conds) {
        combined_p <- res[[depstruct]][[cond]]$combined_p
        
        for (method in methods) {
          p_vec <- combined_p[[method]]
          k <- sum(p_vec <= alpha)
          rate <- k / n_rep
          ci <- binom_ci(k, n_rep, conf = conf)
          
          prefix <- paste(method, cond, sep = "_")
          df[i, paste0(prefix, "_rejections")] <- k
          df[i, paste0(prefix, "_rate")]       <- rate
          df[i, paste0(prefix, "_low")]        <- ci["lower"]
          df[i, paste0(prefix, "_high")]       <- ci["upper"]
        }
      }
      
      sweep_results[[depstruct]] <- df
    }
  }
  
  return(sweep_results)
}

run_all_sweeps <- function(base_vals,
                           range_vals,
                           comb_methods,
                           conf = 0.90,
                           verbose = TRUE) {
  all_results <- list()
  
  for (param_name in names(range_vals)) {
    if (verbose) {
      cat("\n=== Sweeping", param_name, "===\n")
    }
    this_range <- range_vals[[param_name]]
    this_range_list <- list()
    this_range_list[[param_name]] <- this_range
    
    all_results[[param_name]] <- sweep_study(
      base_vals  = base_vals,
      range_vals = this_range_list,
      comb_methods = comb_methods,
      conf       = conf,
      verbose    = verbose
    )
  }
  
  return(all_results)
}

plot_sweep_param <- function(all_sweep_results,
                             sweep_param,
                             base_vals,
                             methods = NULL,
                             depstructs = c("indep", "cs", "ar1"),
                             conf = 0.90,
                             ci_alpha = 0.5,
                             cols = NULL) {
  # all_sweep_results: output of run_all_sweeps(base_vals, range_vals, ...)
  # sweep_param: one of "n_pval", "mu_alt", "prop_alt", "rho"
  # base_vals: list(n_pval, n_rep, mu_alt, prop_alt, rho_cs, rho_ar1, alpha)
  # methods: optional subset of c("fisher","harmonic","cauchy"); if NULL, all three
  # depstructs: which dependence structures to plot
  # conf: CI level (for label only, e.g. 0.90 -> "90% CI")
  
  if (!sweep_param %in% names(all_sweep_results)) {
    stop("sweep_param not found in all_sweep_results")
  }
  
  sweep_results <- all_sweep_results[[sweep_param]]
  
  if (is.null(methods)) {
    methods <- c("fisher", "harmonic", "cauchy", "tippett", "simes")
  }
  
  # nice x-axis labels
  xlab_map <- c(
    n_pval   = "Number of p-values (m)",
    mu_alt   = "Signal mean μ*",
    prop_alt = "Proportion of non-null μ*",
    rho      = "Correlation ρ"
  )
  if (sweep_param %in% names(xlab_map)) {
    xlab <- xlab_map[[sweep_param]]
  } else {
    xlab <- sweep_param
  }
  
  # pretty label for parameter names in the right panel
  param_pretty_label <- function(name) {
    if (name == "n_pval") return("m")
    if (name == "n_rep") return("B")
    if (name == "mu_alt") return("μ*")
    if (name == "prop_alt") return("Pr(μ*)")
    if (name == "rho_cs") return("ρ_cs")
    if (name == "rho_ar1") return("ρ_ar1")
    if (name == "alpha") return("α")
    if (name == "rho") return("ρ")
    return(name)
  }
  
  # constant/variable parameter text lines for right panel
  # constant/variable parameter text lines for right panel
  build_param_lines <- function(sweep_param, base_vals, sweep_results) {
    # determine the range of the swept parameter from one depstruct
    dep_ref <- intersect(c("indep", "cs", "ar1"), names(sweep_results))[1]
    x_vals <- sweep_results[[dep_ref]][[sweep_param]]
    x_range <- range(x_vals, na.rm = TRUE)
    
    param_names <- names(base_vals)
    lines <- character(0)
    
    for (nm in param_names) {
      label <- param_pretty_label(nm)
      val <- base_vals[[nm]]
      if (is.numeric(val)) {
        val_str <- format(signif(val, 3), trim = TRUE)
      } else {
        val_str <- as.character(val)
      }
      
      # which parameters should be marked as "variable"?
      is_var <- (nm == sweep_param) ||
        (sweep_param == "rho" && nm %in% c("rho_cs", "rho_ar1"))
      
      if (is_var) {
        line <- paste0(
          label, " = variable (",
          format(signif(x_range[1], 5), trim = TRUE), "–",
          format(signif(x_range[2], 5), trim = TRUE), ")"
        )
      } else {
        line <- paste0(label, " = ", val_str)
      }
      
      lines <- c(lines, line)
    }
    
    lines
  }
  
  const_lines <- build_param_lines(sweep_param, base_vals, sweep_results)
  alpha_nom <- base_vals$alpha
  
  # figure out which methods actually exist
  dep_ref <- intersect(depstructs, names(sweep_results))[1]
  ref_df <- sweep_results[[dep_ref]]
  methods_use <- methods[
    vapply(methods, function(m) {
      paste0(m, "_null_rate") %in% names(ref_df)
    }, logical(1L))
  ]
  if (length(methods_use) == 0L) {
    stop("No methods found in sweep_results for given sweep_param.")
  }
  
  if(is.null(cols)){
    cols <- RColorBrewer::brewer.pal(max(3, length(methods_use)), "Dark2")[seq_along(methods_use)]  
  }
  ci_label <- paste0(round(conf * 100), "% CI")
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  # 2 rows (null / alt), 4 columns (3 dep panels + legend/params panel)
  par(mfrow = c(2, 4), mar = c(4, 4, 2, 1), oma = c(1,6,5,1))
  
  conds <- c("null", "alt")
  
  # precompute row-wise y-lims so each row shares the same scale
  ylim_row <- list()
  for (cond in conds) {
    all_low_row <- all_high_row <- numeric(0)
    for (depstruct in depstructs) {
      if (!depstruct %in% names(sweep_results)) next
      df_dep <- sweep_results[[depstruct]]
      for (m in methods_use) {
        low_col  <- paste0(m, "_", cond, "_low")
        high_col <- paste0(m, "_", cond, "_high")
        all_low_row  <- c(all_low_row, df_dep[[low_col]])
        all_high_row <- c(all_high_row, df_dep[[high_col]])
      }
    }
    y_min_raw <- min(all_low_row, na.rm = TRUE)
    y_max_raw <- max(all_high_row, na.rm = TRUE)
    if (!is.finite(y_min_raw) || !is.finite(y_max_raw)) {
      y_min_raw <- 0
      y_max_raw <- 1
    }
    range_len <- y_max_raw - y_min_raw
    if (range_len <= 0) range_len <- 0.02
    margin <- 0.05 * range_len
    lower <- max(0, y_min_raw - margin)
    upper <- min(1, y_max_raw + margin)
    # nudge to include 0, 1, alpha if near
    if (lower < 0.1) lower <- 0
    if (upper > 0.9) upper <- 1
    if (alpha_nom < lower && (lower - alpha_nom) < 0.05) lower <- alpha_nom
    if (alpha_nom > upper && (alpha_nom - upper) < 0.05) upper <- alpha_nom
    ylim_row[[cond]] <- c(lower, upper)
  }
  
  for (row_idx in seq_along(conds)) {
    cond <- conds[row_idx]
    ylim <- ylim_row[[cond]]
    
    # three depstruct panels
    for (depstruct in depstructs) {
      if (!depstruct %in% names(sweep_results)) {
        plot.new()
        next
      }
      
      df <- sweep_results[[depstruct]]
      x <- df[[sweep_param]]
      
      plot(x, df[[paste0(methods_use[1], "_", cond, "_rate")]],
           type = "n",
           xlab = xlab,
           ylab = "rejection rate",
           ylim = ylim)
      
      abline(h = alpha_nom, lty = 3, col = "grey60", lwd = 3)
      # pusr <- par("usr")
      # text(pusr[2], alpha_nom,
      #      # paste0("α = ", alpha_nom),
      #      paste0("α"),
      #      pos = 4, xpd = NA, col = "grey60", cex = 1.4)
      
      for (j in seq_along(methods_use)) {
        m <- methods_use[j]
        rate_col <- paste0(m, "_", cond, "_rate")
        low_col  <- paste0(m, "_", cond, "_low")
        high_col <- paste0(m, "_", cond, "_high")
        
        y <- df[[rate_col]]
        y_low <- df[[low_col]]
        y_high <- df[[high_col]]
        
        poly_x <- c(x, rev(x))
        poly_y <- c(y_low, rev(y_high))
        
        polygon(poly_x, poly_y,
                border = NA,
                col = adjustcolor(cols[j], alpha.f = ci_alpha))
        lines(x, y, col = cols[j], lwd = 2)
      }
      
      # row labels on the left (once per row)
      if (depstruct == depstructs[1]) {
        horiz_labs <- c(
          null = "Null Hypothesis\n(all μ = 0)",
          alt  = paste0("Alternative\nHypothesis")
        )
        mtext(horiz_labs[[cond]],
              side = 2, line = 5, cex = 1.25, font = 1, xpd = NA)
      }
      
      # column labels on the top (once per column, using base rho_cs/rho_ar1)
      if (cond == conds[1]) {
        if (sweep_param == "rho") {
          # when sweeping rho, show r[ij] == rho (symbolic) in the column label
          vert_labs <- list(
            indep = expression(atop("Independent", "P-Values")),
            cs    = bquote(atop("Compound Symmetry", r[ij] == rho)),
            ar1   = bquote(atop("AR(1)",              r[ij] == rho^"|i-j|"))
          )
        } else {
          # when not sweeping rho, show fixed numeric r_ij values
          vert_labs <- list(
            indep = expression(atop("Independent", "P-Values")),
            cs    = bquote(atop("Compound Symmetry", r[ij] == .(base_vals$rho_cs))),
            ar1   = bquote(atop("AR(1)",              r[ij] == .(base_vals$rho_ar1)^"|i-j|"))
          )
        }
        
        
        mtext(vert_labs[[depstruct]],
              side = 3, line = 1, cex = 1.25, font = 2, xpd = NA)
      }
      
    }
    
    # rightmost panel
    if (row_idx == 1) {
      plot.new()
      # legend + parameters
      legend_info <- legend(-0.2, 1,
             legend = c(methods_use,
                        ci_label,
                        bquote(alpha == .(base_vals$alpha))),
             col    = c(cols,
                        adjustcolor("grey40", alpha.f = ci_alpha),
                        "grey60"),
             lwd    = c(rep(3, length(methods_use)),
                        NA,
                        4),
             lty    = c(rep(1, length(methods_use)),
                        NA,
                        3),
             pch    = c(rep(NA, length(methods_use)),
                        15,
                        NA),
             pt.cex = c(rep(NA, length(methods_use)),
                        3,
                        NA),
             bty    = "n",
             cex    = 1.2, xpd = NA)
      
      
      # constants + variable info below
      dy <- max(strheight(const_lines))*2
      y_start <- legend_info$rect$top - legend_info$rect$h - dy * 2
      
      text(-0.2, y_start, "Parameters:", adj = c(0, 1), cex = 1.25, font = 2, xpd = NA)
      y <- y_start - dy
      for (line in const_lines) {
        text(-0.2, y, line, adj = c(0, 1), cex = 1.2, xpd = NA,
             font = ifelse(grepl("variable", line), 3, 1),
             col = ifelse(grepl("variable", line), 2, 1))
        y <- y - dy
      }
    } else {
      plot.new()
    }
  }
}


#### basic analysis ####
comb_methods <- list(
  fisher   = combine_fisher,
  harmonic = combine_harmonic,
  cauchy   = combine_cauchy,
  tippett     = combine_tippett,
  simes = combine_simes,
  stouffer = combine_stouffer,
  brown_cs = combine_brownCS
)

# example: actually run it
base_vals <- list(
  n_pval = 2^9,
  n_rep = 5E3,
  mu_alt = 1,
  prop_alt = 0.1,
  rho_cs = 0.75,
  rho_ar1 = 0.95,
  alpha = 0.1
)
results <- run_full_study(
  B = base_vals$n_rep,   # number of simulation replicates
  m = base_vals$n_pval,      # dimension
  mu_alt = base_vals$mu_alt,
  prop_alt = base_vals$prop_alt,
  rho_cs = base_vals$rho_cs,
  rho_ar1 = base_vals$rho_ar1,
  alpha = base_vals$alpha,
  comb_methods = comb_methods
)

#### plotting scatter-hists ####
source("~/scripts/minor_scripts/postdoc/scatter_hist.R")
plot_as_heatmap <- F
figures_dir <- "~/pval_meta-anal/"
subfigures_dir <- paste0(figures_dir, "subfigures/")
figure_paths <- c()

for(use_log10 in c(F,T)){
  
  
  depstructs <- c("indep", "cs", "ar1")
  conds <- c("null", "alt")
  
  figure_path <- paste0(subfigures_dir, "/meta-analysis_scatter-hist... ", 
                        paste0("n_pval = ", base_vals$n_pval, ", "),
                        paste0("mu_alt = ", base_vals$mu_alt, ", "),
                        paste0("prop_alt = ", base_vals$prop_alt, ", "),
                        paste0("rho_cs = ", base_vals$rho_cs, ", "),
                        paste0("rho_ar1 = ", base_vals$rho_ar1, ", "),
                        paste0("alpha = ", base_vals$alpha, ", "),
                        ifelse(use_log10, "log10transformed", "natural scale"),
                        ".pdf")
  figure_paths <- c(figure_paths, figure_path)
  cairo_pdf(filename = figure_path, width = 12, height = 8)
  par(mfcol = c(2,3), mar = c(3,3,0,0), oma = c(0,6,10,0))
  
  for(depstruct in depstructs){
    for(cond in conds){
      harmonic <- results[[depstruct]][[cond]]$combined_p$harmonic
      cauchy <-  results[[depstruct]][[cond]]$combined_p$cauchy
      alpha_thresh <- c(base_vals$alpha, base_vals$alpha)
      xbreaks <- ybreaks <- 0:10/10
      xlab <- "harmonic mean p-value"
      ylab <- "cauchy method p-value"
      if(use_log10){
        harmonic <- log10(harmonic)
        cauchy <-  log10(cauchy)
        alpha_thresh <- log10(alpha_thresh)
        xbreaks <- ybreaks <- pretty(c(harmonic, cauchy), n = 10)
        xlab <- "harmonic mean log₁₀(p-value)"
        ylab <- "cauchy method log₁₀(p-value)"
      }
      
      scatter_hist(x = harmonic, y = cauchy, col.hpt = adjustcolor(2,0.9),
                   shade_right_tail = F,
                   highlight_pt = alpha_thresh, shade_hist_tails = T, plot_pt = F,
                   target.ticks = 10, xbreaks = xbreaks, ybreaks = ybreaks,
                   xlab = xlab, ylab = ylab, asp = 1, 
                   add_one_to_one_line = T, cex.lab = 1, equal_axes = T,
                   plot_as_heatmap = plot_as_heatmap, col.pts = NULL, kde_heatmap = F)
      
      #check tl corner entries -- not meant to be run inside function
      check_tl <- F
      if(check_tl){
        tl_hits <- harmonic < 0.1 & cauchy > 0.9
        ntl_hits <- harmonic < 0.1 & cauchy < 0.1
        # ntl_hits <- !tl_hits
        tl_corner_ps <- results$indep$alt$settings$p_mat[tl_hits,]
        ntl_corner_ps <- results$indep$alt$settings$p_mat[ntl_hits,]
        
        mean(apply(tl_corner_ps, 1, var))
        mean(apply(ntl_corner_ps, 1, var))
        
        mean(apply(tl_corner_ps, 1, mean))
        mean(apply(ntl_corner_ps, 1, mean))
        
        #try minimum p value
        num_mins <- c(1, 2, 4, 8, 16)
        for(num_min in num_mins){
          stat_tlps <- apply(tl_corner_ps, 1, function(tmp) mean(sort(tmp)[1:num_min]))
          stat_ntlps <- apply(ntl_corner_ps, 1, function(tmp) mean(sort(tmp)[1:num_min]))
          
          breaks <- seq(min(c(stat_tlps, stat_ntlps)),
                        max(c(stat_tlps, stat_ntlps)),
                        length.out = 51)
          h_tlps <- hist(stat_tlps, breaks = breaks, plot = F)
          h_ntlps <- hist(stat_ntlps, breaks = breaks, plot = F)
          hist(stat_tlps, freq = F, breaks = breaks, 
               col = adjustcolor("orange", 0.5),
               xlim = range(c(stat_tlps, stat_ntlps)),
               ylim = range(c(h_tlps$density, h_ntlps$density)),
               xlab = ifelse(num_min == 1, "smallest p-value", 
                             paste0("mean of smallest ", num_min, " p-values")),
               main = "")
          hist(stat_ntlps, breaks = breaks,
               freq = F, add = T, col = adjustcolor("blue", 0.5))
          legend("topright",
                 legend = c("harmonic < 0.1 & cauchy > 0.9",
                            "harmonic < 0.1 & cauchy < 0.1"),
                 fill   = c(adjustcolor("orange", 0.5),
                            adjustcolor("blue",   0.5)),
                 border = NA,
                 bty    = "n",
                 cex    = 0.9)
          
        }
        
        
      }
      
      if(depstruct == depstructs[1]){
        horiz_labs <- c("null" = paste0("Null Hypothesis\n(all μ = 0)"),
                        "alt" = paste0("Alternative Hypothesis\n(μ* = ", base_vals$mu_alt, ", Pr(μ*) = ", base_vals$prop_alt, ")"))
        mtext(paste0(horiz_labs[cond]), outer = F, xpd = NA, side = 2, line = 4, cex = 1.5, font = 2)
      }
      if(cond == conds[1]){
        vert_labs <- c("indep" = paste0("Independent\nP-Values"),
                       "cs" = paste0("Compound Symmetry\n(rᵢⱼ = ", base_vals$rho_cs, ")"),
                       "ar1" = paste0("AR(1)\n(rᵢⱼ = ", base_vals$rho_ar1, ")"))
        mtext(paste0(vert_labs[depstruct]), outer = F, xpd = NA, side = 3, line = 0, cex = 1.5, font = 2)
      }
    }
  }
  mtext(paste0("Comparing Aggregation Methods for ", base_vals$n_pval," P-Values Generated Using Φ⁻¹(MVN) at α = ", base_vals$alpha), 
        outer = T, xpd = NA, side = 3, line = 6, cex = 1.65)
  
  
  dev.off()
  
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
  
}

all_figures_path <- paste0(figures_dir, "meta-analysis_scatter-hist... ", 
                           paste0("n_pval = ", base_vals$n_pval, ", "),
                           paste0("mu_alt = ", base_vals$mu_alt, ", "),
                           paste0("prop_alt = ", base_vals$prop_alt, ", "),
                           paste0("rho_cs = ", base_vals$rho_cs, ", "),
                           paste0("rho_ar1 = ", base_vals$rho_ar1, ", "),
                           paste0("alpha = ", base_vals$alpha),
                           ".pdf")
pdftools::pdf_combine(input = figure_paths, output = all_figures_path)

#### sweep analysis ####

#base values
base_vals <- list(
  n_pval = 2^9,
  n_rep = 1E4,
  mu_alt = 1,
  prop_alt = 0.1,
  rho_cs = 0.75,
  rho_ar1 = 0.95,
  alpha = 0.1
)

#ranges of values to check
range_vals <- list(
  n_pval = 2^(1:10),
  mu_alt = (0:10)/4,
  prop_alt = c((0.5)^(0:10), 0),
  rho = c(0, 1 - 0.5^(1:8))
)

sweep_results <- run_all_sweeps(
  base_vals  = base_vals,
  range_vals = range_vals,
  conf = 0.90,
  comb_methods = comb_methods,
  verbose = TRUE
)

#### plot sweep results ####
names(comb_methods)
methods_inds <- c(2,3)
all_methods <- names(comb_methods)
all_cols <- RColorBrewer::brewer.pal(8, "Dark2")
methods <- all_methods[methods_inds]
cols <- all_cols[methods_inds]
figure_paths <- c()
for(sp in names(range_vals)){
  
  figure_path <- paste0(subfigures_dir, "/meta-analysis_coverage-at-alpha-=",
                        base_vals$alpha, "_", paste0(methods, collapse = "-"), "_", sp, ".pdf")
  figure_paths <- c(figure_paths, figure_path)
  cairo_pdf(filename = figure_path, width = 10, height = 5)

  plot_sweep_param(
    all_sweep_results = sweep_results,
    methods = methods,
    sweep_param       = sp,
    base_vals         = base_vals,
    cols = cols
  )
  
  dev.off()
  
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
  
}

all_figures_path <- paste0(figures_dir, "/meta-analysis_coverage-at-alpha-=",
                           base_vals$alpha, ".pdf")
pdftools::pdf_combine(input = figure_paths, output = all_figures_path)


#### plots for presentation ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")
# identity, cs, and ar(1) corrmats
p <- 20
rho_cs <- 0.5
rho_ar1 <- 0.95
R_titles <- list(
  iden = bquote("Identity:" ~ r[ij] == 0 ~~ "(" * i != j * ")"),
  cs = bquote(
    "Compound symmetry:" ~ 
      r[ij] == rho ~~ "(" * i != j * ")" ~~ " =" ~~ .(rho_cs)
  ),
  ar1 = bquote(
    "AR(1):" ~ 
      r[ij] == rho^{ group("|", i - j, "|") } ~~ 
      " =" ~~ .(rho_ar1)^{ group("|", i - j, "|") }
  )
)

for(Rn in names(Rs)){
  figure_path <- paste0(subfigures_dir, "/presentations_matrix-", Rn, ".pdf")
  cairo_pdf(filename = figure_path, width = 5, height = 5)
  par(mar = c(1,1,6,1))
  my_heatmap(Rs[[Rn]], diag_matters = T, asp = 1, plot_diagonal_labels = F, 
             legend_title = expression(r[ij]), plot_labels = F, 
             plot_numbers = T, plot_guiding_lines = F, mat_size_rule = "abs", 
             number_col = "white", outer_border = T, main = R_titles[[Rn]])
  dev.off()
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
}


#some extra tiny plots
pv <- rbeta(1E4, 1, 1.05)
hist(pv, xlab = "Uniform(0,1) p-values", main = "", freq = F, col = "white", border = "white")
rect(0,0,1,1, col = "grey")


qmax <- 0.999
maxexp <- qexp(qmax, 0.5)
curve(dexp(x, rate = 0.5), xlim = c(0,maxexp), xlab = "x", ylab = "density", lwd = 2, 
      main = "exponential(rate = 0.5)", cex.main = 1)
x <- seq(0, maxexp, length.out = 1E3)
d <- dexp(x, rate = 0.5)
polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")

for(np in 1:3){
  maxchi2 <- qchisq(qmax, df = 2 * np)
  curve(dchisq(x, df = 2 * np), xlim = c(0,maxchi2), xlab = "x", ylab = "density", lwd = 2, 
        main = paste0("χ²(df = ", 2 * np, ")"), cex.main = 1)
  x <- seq(0, maxchi2, length.out = 1E3)
  d <- dchisq(x, df = 2 * np)
  polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")
}


curve(dnorm(x, mean = 1), xlim = c(-4,4), xlab = "x", ylab = "density", lwd = 2, 
      main = "normal(1,1)", cex.main = 1)
x <- seq(-4, 4, length.out = 1E3)
d <- dnorm(x, 1, 1)
polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")



hmean <- function(p) {
  1/(mean(1/p))
}
np <- 5000
hist(replicate(1E4, log10(hmean(runif(np)))), freq = F, breaks = 120,
     main = paste0("num. p-vals = ", np))
abline(h = 1, lty = 2, lwd = 2)
text(par("usr")[2], 1, labels = "1", pos = 4, xpd = NA)

#### extract dependence ####

#can we get z matrix from single dataset?
nrep <- 2E2
n <- 100
p <- 5E2
r <- 0.9
R <- outer(1:p, 1:p, function(i, j) r^abs(i - j))

#constant params
b <- rnorm(p) * 0
x <- rnorm(n)

#variable params
sim_dat <- function(){
  e <- simulate_X(n, p, rep(0, p), "ar1", r)
  y <- outer(x, b, "*") + e
  return(y)
}
sim_res <- function(){
  y <- sim_dat()
  out <- as.data.frame(do.call(rbind, lapply(1:p, function(i) summary(lm(y[,i] ~ x))$coefficients[2,c(1,4)])))
  return(out)
}

outs <- parallel::mclapply(1:nrep, function(nri) sim_res(), mc.cores = 8)
ests <- do.call(cbind, lapply(outs, function(res) res$Estimate))
pvals <- do.call(cbind, lapply(outs, function(res) res$`Pr(>|t|)`))
cor_ests <- cor(t(ests))
cor_pvals <- cor(t(pvals))
cor_log10pvals <- cor(t(log10(pvals)))
cor_zs <- cor(t(qnorm(pvals)))
cor_szs <- cor(t(sign(ests) * qnorm(pvals)))

plot(R[upper.tri(R)], cor_ests[upper.tri(cor_ests)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_pvals[upper.tri(cor_pvals)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_log10pvals[upper.tri(cor_log10pvals)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_zs[upper.tri(cor_zs)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

cor(R[upper.tri(R)], cor_ests[upper.tri(cor_ests)])
cor(R[upper.tri(R)], cor_pvals[upper.tri(cor_pvals)])
cor(R[upper.tri(R)], cor_log10pvals[upper.tri(cor_log10pvals)])
cor(R[upper.tri(R)], cor_zs[upper.tri(cor_zs)])

#estimate covariance matrix by permutation

compute_tstats <- function(y, x) {
  # y: n x p matrix (one dataset)
  # x: length-n vector (predictor)
  # returns: length-p vector of t-statistics for slope in lm(y[,j] ~ x)
  
  n <- nrow(y)
  p <- ncol(y)
  tstats <- numeric(p)
  
  for (j in seq_len(p)) {
    fit <- lm(y[, j] ~ x)
    coefs <- summary(fit)$coefficients
    tstats[j] <- coefs[2, "t value"]
  }
  
  tstats
}

estimate_null_covariance <- function(y, x,
                                     B = 500,
                                     transform = c("t", "z", "p", "z_from_p"),
                                     mc.cores = 8) {
  # y: n x p matrix (one observed dataset)
  # x: length-n vector (predictor)
  # B: number of null resamples (permutations of x)
  # transform:
  #   "t"         = use raw t-statistics
  #   "z"         = treat t as z (large n)
  #   "p"         = work with two-sided p-values
  #   "z_from_p"  = signed z from two-sided p
  # mc.cores: number of cores for mclapply
  
  transform <- match.arg(transform)
  
  n <- nrow(y)
  p <- ncol(y)
  df <- n - 2
  
  # one null resample → one length-p vector of transformed stats
  one_null_draw <- function(b) {
    x_perm <- sample(x)  # permutation of labels
    tstats <- compute_tstats(y, x_perm)
    
    if (transform == "t" || transform == "z") {
      return(tstats)
    }
    
    # two-sided p-values from t
    pvals <- 2 * pt(-abs(tstats), df = df)
    
    if (transform == "p") {
      return(pvals)
    }
    
    # signed z from two-sided p
    z_signed <- sign(tstats) * qnorm(pvals / 2, lower.tail = FALSE)
    z_signed
  }
  
  # parallel null draws: list of length B, each element length p
  stats_list <- parallel::mclapply(seq_len(B), one_null_draw, mc.cores = mc.cores)
  
  # p x B matrix
  stats_mat <- do.call(cbind, stats_list)
  
  # covariance and correlation across tests (rows) under null
  cov_mat <- t(cov(stats_mat))
  cor_mat <- t(cor(stats_mat))
  
  list(
    stats_mat = stats_mat,  # p x B, optional to inspect
    cov = cov_mat,          # p x p covariance of transformed stats under null
    cor = cor_mat           # p x p correlation of transformed stats under null
  )
}


# estimate covariance of null z-statistics via permutation of x
n <- 100
p <- 5E2
r <- 0.9
R <- outer(1:p, 1:p, function(i, j) r^abs(i - j))

#constant params
b <- rnorm(p) * 0
x <- rnorm(n)
y_obs <- sim_dat()

#estimate covariance matrix
null_est <- estimate_null_covariance(
  y         = y_obs,
  x         = x,
  B         = 500,
  transform = "z_from_p"
)

t_obs <- compute_tstats(y_obs, x)
df    <- nrow(y_obs) - 2
p_two <- 2 * pt(-abs(t_obs), df = df)
z_obs <- sign(t_obs) * qnorm(p_two / 2, lower.tail = FALSE)

combine_stouffer_perm(z_obs, null_est$stats_mat)

#quick check of local correlation thing
n <- 1E7
rho <- 0.95
X <- simulate_X(B = n, m = 2, mu_vec = c(0,0), 
                corr_structure = "cs", rho = rho)
cor(X)[1,2]
nbreaks <- 10
eps <- 1E-6
breaks <- seq(-5, 5, length.out = nbreaks + 1)
midp <- breaks[-1] - diff(breaks)/2
cutX1 <- as.numeric(cut(X[,1], breaks))
cutX2 <- as.numeric(cut(X[,2], breaks))
outX1 <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX1==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
outX2 <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX2==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
outXb <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX1==i) & (cutX2==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
plot(midp, outX1$cor, type = "l", ylim = range(c(outX1$cor,
                                                 outX2$cor,
                                                 outXb$cor,
                                                 rho)),
     col = "red")
lines(midp, outXb$cor, col = "purple")
lines(midp, outX2$cor, col = "blue")
abline(h=rho, lty=2, col=adjustcolor(1,0.5), lwd = 2)

#### quick tests ####
# m <- 1E3
# rho <- 0.9
# mu_vec <- rep(1, m)
# sigma2 <- 1
# corr_structure <- c("cs", "ar1")[1]
# x <- simulate_X(1, 
#            m = m, 
#            mu_vec = mu_vec, 
#            rho = rho, 
#            corr_structure = corr_structure) * sqrt(sigma2)
# p <- pnorm(x)
# hist(p, breaks = 0:20/20)
# est <- estimate_gaussian_struct_cs(p)
# est$coef

#### end ####
