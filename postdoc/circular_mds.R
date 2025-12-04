xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

circ_mds <- function(dists, maxit = 2000, n_starts = 4, init_sigma = NULL,
                     trace = F) {
  
  as_full <- function(d) {
    if (inherits(d, "dist")) return(as.matrix(d))
    if (is.matrix(d)) return(d)
    stop("dists must be a matrix or a 'dist'")
  }
  D <- as_full(dists)
  n <- nrow(D); diag(D) <- 0
  if (n < 3) stop("need at least 3 points")
  if (max(abs(D - t(D))) > 1e-10) stop("distance matrix must be symmetric")
  if (any(D < 0)) stop("distances must be nonnegative")
  ut <- which(upper.tri(D), arr.ind = TRUE)
  
  shortest_arc <- function(theta) {
    raw <- outer(theta, theta, "-")
    m <- abs(raw) %% (2*pi)
    pmin(m, 2*pi - m)
  }
  alpha_hat <- function(theta) {
    A <- shortest_arc(theta)
    a <- A[upper.tri(A)]; d <- D[upper.tri(D)]
    den <- sum(a*a); if (den <= 0) return(0)
    max(sum(a*d)/den, 0)
  }
  stress_val <- function(theta) {
    A <- shortest_arc(theta)
    a <- A[upper.tri(A)]; d <- D[upper.tri(D)]
    ah <- alpha_hat(theta); sum((ah*a - d)^2)
  }
  
  # --- better initializations ---
  init_spectral <- function() {
    # choose sigma by median of nonzero distances if not provided
    if (is.null(init_sigma)) {
      nz <- D[D > 0]
      s <- if (length(nz)) median(nz) else 1
    } else s <- init_sigma
    S <- exp(-(D / s)^2)
    diag(S) <- 0
    # normalized laplacian eigenmaps (top 2 nontrivial eigenvectors of D^-1/2 S D^-1/2)
    ddeg <- rowSums(S)
    ddeg[ddeg == 0] <- 1
    Dm12 <- 1 / sqrt(ddeg)
    K <- (Dm12 * S) * rep(Dm12, each = n)
    es <- eigen((K + t(K))/2, symmetric = TRUE, only.values = FALSE)
    # skip the largest eigenvector (~constant); take next two
    if (ncol(es$vectors) < 3) {
      if (trace) cat("# spectral init: fallback to random (insufficient eigens)\n")
      return(runif(n, 0, 2*pi))
    }
    v1 <- es$vectors[,2]; v2 <- es$vectors[,3]
    ang <- atan2(v2, v1); (ang + 2*pi) %% (2*pi)
  }
  init_equal_spaced_from_order <- function(order_idx) {
    ang <- numeric(n)
    ang[order_idx] <- seq(0, 2*pi*(1 - 1/n), length.out = n)
    ang
  }
  init_random <- function() runif(n, 0, 2*pi)
  
  # build a small pool of candidates
  cand <- list()
  # spectral
  th_s <- init_spectral()
  cand[[length(cand)+1]] <- th_s
  # equally spaced using spectral order (circular seriation bootstrap)
  ord <- order(th_s)
  cand[[length(cand)+1]] <- init_equal_spaced_from_order(ord)
  # a few random restarts
  for (k in seq_len(max(0, n_starts - length(cand)))) cand[[length(cand)+1]] <- init_random()
  
  # objective for optim with anchored angle to remove rotation
  anchor <- 1L
  obj <- function(par, n, anchor) {
    theta <- numeric(n); theta[anchor] <- 0; theta[-anchor] <- par %% (2*pi)
    st <- stress_val(theta)
    if (trace && runif(1) < 0.02) cat("# debug: stress=", format(st, digits = 7), "\n")
    st
  }
  
  best <- NULL
  for (i in seq_along(cand)) {
    th0 <- cand[[i]]
    th0 <- (th0 - th0[anchor]) %% (2*pi)  # align anchor to 0
    par0 <- th0[-anchor]
    opt <- optim(par0, obj, method = "L-BFGS-B", lower = rep(0, n-1), upper = rep(2*pi, n-1),
                 control = list(maxit = maxit, factr = 1e7), n = n, anchor = anchor)
    theta_fit <- numeric(n); theta_fit[anchor] <- 0; theta_fit[-anchor] <- opt$par %% (2*pi)
    sse <- stress_val(theta_fit)
    if (trace) cat("# start", i, ": final_sse=", format(sse, digits = 7), "\n")
    if (is.null(best) || sse < best$sse) best <- list(theta = theta_fit, sse = sse)
  }
  
  theta_fit <- best$theta
  ah <- alpha_hat(theta_fit)
  Afit <- shortest_arc(theta_fit)
  a <- Afit[upper.tri(Afit)]; d <- D[upper.tri(D)]
  pred <- ah * a
  rss <- sum((pred - d)^2); tss <- sum((d - mean(d))^2)
  r2 <- 1 - rss/tss; cor_ud <- cor(pred, d)
  coords <- data.frame(x = cos(theta_fit), y = sin(theta_fit))
  
  if (trace) {
    cat("# finished: stress=", format(rss, digits = 7),
        " alpha=", format(ah, digits = 6),
        " R2=", format(r2, digits = 5),
        " cor(upper)=", format(cor_ud, digits = 5), "\n", sep = "")
  }
  
  list(theta = theta_fit,
       coords = coords,
       alpha = ah,
       stress = rss,
       r2 = r2,
       correlation_upper = cor_ud,
       fitted_upper = pred,
       target_upper = d)
}

# sunburst labelling with leaders to the text's nearest radial edge (not the center)
label_sunburst <- function(coords, labels,
                           r_anchor = 1.08, r_text = 1.16,
                           min_gap_deg = 8, cex = 0.8,
                           rotate = c("none","radial","tangent"),
                           upright_left = TRUE,
                           leader_extend = 0.04, leader_min = 0.06,
                           link_to_text = TRUE,
                           line_lwd = 0.8, line_col = "gray40",
                           point_pch = 16, point_cex = 0.7, xyt = c(0,0)) {
  rotate <- match.arg(rotate)
  if (is.null(colnames(coords))) colnames(coords) <- c("x","y")
  x <- coords[, "x"]; y <- coords[, "y"]; n <- length(x)
  
  # base angles in [0, 2pi)
  theta <- atan2(y, x); theta[theta < 0] <- theta[theta < 0] + 2*pi
  
  # minimal-perturbation spacing
  ord <- order(theta); th <- theta[ord]
  min_gap <- min_gap_deg * pi/180
  for (i in 2:n) if (th[i] - th[i-1] < min_gap) th[i] <- th[i-1] + min_gap
  span <- th[n] - th[1]
  if (span > 2*pi) th <- th[1] + (th - th[1]) * (2*pi/span)
  theta_lab <- th[order(ord)]
  
  # geometry helpers
  urx <- cos(theta_lab); ury <- sin(theta_lab)                   # radial unit vectors
  r_leader_end <- pmax(r_anchor + leader_extend, 1 + leader_min) # enforce min length
  x_lead <- r_leader_end * urx; y_lead <- r_leader_end * ury
  x_text0 <- r_text * urx; y_text0 <- r_text * ury
  
  # draw points
  points(x + xyt[1], y + xyt[2], pch = point_pch, cex = point_cex)
  
  # base rotation (before upright adjustment)
  srt_base <- switch(rotate,
                     none    = rep(0, n),
                     radial  = theta_lab * 180/pi,
                     tangent = (theta_lab * 180/pi + 90) %% 360)
  
  # outward padding (keep labels outside the circle)
  pad <- strwidth(labels, cex = cex) * 0.55
  x_lab <- x_text0 + pad * urx
  y_lab <- y_text0 + pad * ury
  
  # upright flip for rotated text
  srt_use <- srt_base
  if (upright_left && rotate != "none") {
    ang <- ((srt_base + 180) %% 360) - 180          # (-180, 180]
    flip <- abs(ang) > 90
    ang[flip] <- ang[flip] - 180 * sign(ang[flip])
    srt_use <- (ang + 360) %% 360
  }
  
  # compute the **nearest radial edge point** of each text label
  # measure text size (in user coords) and project onto radial direction
  w <- strwidth(labels, cex = cex)
  h <- strheight(labels, cex = cex)
  
  # distance from label center to its inner radial edge:
  # - radial: half width along radial axis
  # - tangent: half height along radial axis
  # - none: projection of half-rect onto the radial direction
  edge_inset <- switch(rotate,
                       radial  = 0.5 * w,
                       tangent = 0.5 * h,
                       none    = 0.5 * (w * abs(urx) + h * abs(ury))
  )
  
  # nearest edge (toward origin) = center - inset * radial_unit
  x_edge <- x_lab - edge_inset * urx
  y_edge <- y_lab - edge_inset * ury
  
  # leaders: point -> leader end (beyond anchor)
  segments(x + xyt[1], y + xyt[2], x_lead + xyt[1], y_lead + xyt[2], 
           lwd = line_lwd, col = line_col)
  
  # leaders: leader end -> text nearest radial edge (clear association)
  if (link_to_text) {
    segments(x_lead + xyt[1], y_lead + xyt[2], x_edge + xyt[1], y_edge + xyt[2], 
             lwd = line_lwd, col = line_col)
  }
  
  # horizontal justification for rotate="none"
  adj_h <- ifelse(urx >= 0, 0, 1)
  
  # draw labels
  for (i in seq_len(n)) {
    if (rotate == "none") {
      text(x_lab[i] + xyt[1], y_lab[i] + xyt[2], labels[i], cex = cex, adj = c(adj_h[i], 0.5), xpd = NA)
    } else {
      text(x_lab[i] + xyt[1], y_lab[i] + xyt[2], labels[i], cex = cex, adj = c(0.5, 0.5), srt = srt_use[i], xpd = NA)
    }
  }
}

flip_signs <- function(R, objective = c("sum", "count"), max_iter = 100000, verbose = TRUE, seed = NULL) {
  objective <- match.arg(objective)
  if (!isSymmetric(R)) stop("R must be symmetric")
  p <- nrow(R)
  if (!is.null(seed)) set.seed(seed)
  
  # precompute the matrix used in the quadratic form d' M d (zero diagonal)
  if (objective == "sum") {
    M <- R
    diag(M) <- 0
  } else {
    M <- sign(R)
    diag(M) <- 0
  }
  
  # spectral init: sign of leading eigenvector (ties -> +1)
  ev <- eigen(M, symmetric = TRUE, only.values = FALSE)
  d <- ifelse(ev$vectors[, 1] >= 0, 1L, -1L)
  
  # helper: objective value (upper triangle)
  obj_value <- function(d) {
    # d' M d / 2 equals sum_{i<j} M_ij d_i d_j
    as.numeric(t(d) %*% M %*% d) / 2
  }
  
  # greedy coordinate ascent: flip when it increases the objective
  # flipping i changes objective by: Δ = -2 * d_i * sum_{j != i} M_ij d_j
  best <- obj_value(d)
  if (verbose) cat("# init objective:", best, "\n")
  
  it <- 0L
  improved <- TRUE
  while (improved && it < max_iter) {
    improved <- FALSE
    # compute "field" h_i = sum_j M_ij d_j
    h <- as.numeric(M %*% d)
    # candidates where flipping helps: d_i * h_i < 0
    to_flip <- which(d * h < 0)
    if (length(to_flip) > 0) {
      # flip the one with largest improvement first (steepest ascent)
      gains <- -2 * d[to_flip] * h[to_flip]
      i <- to_flip[which.max(gains)]
      d[i] <- -d[i]
      best <- best + max(gains)
      improved <- TRUE
      it <- it + 1L
      if (verbose && it %% 10L == 0L) cat("# iter", it, "objective:", best, "\n")
    }
  }
  if (verbose) cat("# finished after", it, "iterations. objective:", best, "\n")
  
  # construct flipped matrix and summaries
  D <- diag(d)
  R_flip <- D %*% R %*% D
  ut <- upper.tri(R_flip)
  avg_corr <- mean(R_flip[ut])
  frac_pos  <- mean(R_flip[ut] > 0)
  
  list(
    signs = d,
    R_flipped = R_flip,
    objective_value = best,
    average_correlation = avg_corr,
    fraction_positive = frac_pos
  )
}

make_cov_rhs <- function(covars, interact_with = "sex") {
  if (!interact_with %in% covars) return(paste(covars, collapse = " + "))
  others <- setdiff(covars, interact_with)
  if (length(others) == 0) return(interact_with)  # only 'sex' in covars
  paste0(interact_with, " * (", paste(others, collapse = " + "), ")")
}

#### viewpoint helpers ####

bezier_curve <- function(x0, y0, x1, y1,
                         p = 0.5, n = 128, k = 1,
                         ends = c("flat", "steep")[1],
                         col = 1, return_points = FALSE, debug = FALSE, ...) {
  # choose behavior: "flat" = horizontal tangents at ends (original);
  #                  "steep" = vertical tangents at ends (the “other direction”)
  if(is.numeric(ends)){
    ends <- c("flat", "steep")[ends]
  }
  ends <- match.arg(ends)  
  
  
  # handle vector inputs (recurse over elementwise args)
  arg_vals <- as.list(environment())
  nargs <- length(arg_vals)
  arg_lens <- sapply(arg_vals, length)
  if (any(arg_lens > 1)) {
    for (i in 1:max(arg_lens)) {
      arg_inds <- (i - 1) %% arg_lens + 1
      arg_vals_i <- lapply(setNames(1:nargs, names(arg_vals)),
                           function(j) arg_vals[[j]][arg_inds[j]])
      do.call(bezier_curve, arg_vals_i)
    }
    return(invisible(NULL))
  }
  
  # control points
  if (ends == "flat") {
    # horizontal tangents at endpoints (original behavior)
    control_x1 <- x0 + p * (x1 - x0) * k
    control_y1 <- y0
    control_x2 <- x1 - (1 - p) * (x1 - x0) * k
    control_y2 <- y1
  } else {
    # vertical tangents at endpoints (steep at start/end)
    control_x1 <- x0
    control_y1 <- y0 + p * (y1 - y0) * k
    control_x2 <- x1
    control_y2 <- y1 - (1 - p) * (y1 - y0) * k
  }
  
  if (debug) {
    cat("# control points:\n")
    cat(sprintf("# P0=(%.4f, %.4f)  P1=(%.4f, %.4f)\n", x0, y0, control_x1, control_y1))
    cat(sprintf("# P2=(%.4f, %.4f)  P3=(%.4f, %.4f)\n", control_x2, control_y2, x1, y1))
    # endpoint derivatives of a cubic Bézier: 3*(P1-P0) at t=0; 3*(P3-P2) at t=1
    d0x <- 3 * (control_x1 - x0); d0y <- 3 * (control_y1 - y0)
    d1x <- 3 * (x1 - control_x2); d1y <- 3 * (y1 - control_y2)
    cat(sprintf("# dB/dt at t=0  ≈ (%.4f, %.4f)\n", d0x, d0y))
    cat(sprintf("# dB/dt at t=1  ≈ (%.4f, %.4f)\n", d1x, d1y))
  }
  
  # sample curve
  t <- seq(0, 1, length.out = n)
  curve_x <- (1 - t)^3 * x0 +
    3 * (1 - t)^2 * t * control_x1 +
    3 * (1 - t) * t^2 * control_x2 +
    t^3 * x1
  curve_y <- (1 - t)^3 * y0 +
    3 * (1 - t)^2 * t * control_y1 +
    3 * (1 - t) * t^2 * control_y2 +
    t^3 * y1
  
  # draw
  lines(curve_x, curve_y, col = col, ...)
  
  if (return_points) {
    return(invisible(list(x = curve_x, y = curve_y)))
  } else {
    return(invisible(NULL))
  }
}

pretty_large_number <- function(num, nsf = NA) {
  if(length(num) > 1){
    return(sapply(num, pretty_large_number, nsf = nsf))
  }
  
  if (is.na(num)) return(NA)
  
  if (num == 0) return("0")
  
  if(num < 1) return(as.character(round(num, 2)))
  
  suffixes <- c("", "K", "M", "B", "T")
  index <- floor(log10(num) / 3)
  
  divisor <- 10^(index * 3)
  rounded_num <- round(num / divisor, 1)
  if(!is.na(nsf) & !grepl("\\.", rounded_num)){
    rounded_num <- paste0(rounded_num, ".0")
  }
  
  return(paste0(rounded_num, suffixes[index + 1]))
}

fit_one <- function(df, outcome, covars, newdata, min_n = 5) {
  ok <- !is.na(df[[outcome]])
  for (v in covars) ok <- ok & !is.na(df[[v]])
  n_ok <- sum(ok)
  if (n_ok < min_n) {
    return(list(mu = NA_real_, res = numeric(0), sd = NA_real_, n = 0L))
  }
  fmla <- as.formula(paste(outcome, "~", paste(covars, collapse = " + ")))
  fit <- lm(fmla, data = df[ok, , drop = FALSE])
  mu  <- as.numeric(predict(fit, newdata = newdata))
  res <- residuals(fit)
  list(mu = mu, res = res, sd = sd(res), n = length(res))
}

# incorporating heteroskedasticity
# fit_one(): location–scale mean/dispersion, same return shape as before
# backend = "auto" tries dglm -> gamlss -> twostage -> homoskedastic
fit_one <- function(df, outcome, covars, newdata, min_n = 5,
                    backend = c("auto","dglm","gamlss","twostage","homoskedastic"),
                    verbose = FALSE,
                    interact_with = "",          # new: who to interact with
                    interact_in_scale = TRUE) {     # new: also add interactions in scale?
  backend <- match.arg(backend)
  
  ok <- !is.na(df[[outcome]])
  for (v in covars) ok <- ok & !is.na(df[[v]])
  if (sum(ok) < min_n) return(list(mu = NA_real_, res = numeric(0), sd = NA_real_, n = 0L))
  
  dsub <- df[ok, , drop = FALSE]
  
  # build RHS once; includes sex + others + all sex:other terms
  cov_text <- make_cov_rhs(covars, interact_with = interact_with)
  
  fmla_mu  <- as.formula(paste(outcome, "~", cov_text))
  fmla_sc  <- if (interact_in_scale) {
    as.formula(paste("y_sc ~", cov_text))
  } else {
    as.formula(paste("y_sc ~", paste(covars, collapse = " + ")))
  }
  
  .homosked <- function() {
    fit <- lm(fmla_mu, data = dsub)
    mu  <- as.numeric(predict(fit, newdata = newdata))
    res <- residuals(fit)
    list(mu = mu, res = as.numeric(res), sd = sd(res), n = length(res))
  }
  
  .twostage <- function() {
    fit_mu <- lm(fmla_mu, data = dsub)
    mu     <- as.numeric(predict(fit_mu, newdata = newdata))
    r      <- as.numeric(residuals(fit_mu))
    
    eps <- 1e-12 + 1e-8 * var(r, na.rm = TRUE)
    y_sc <- log(pmax(r^2, eps))
    
    fit_sc <- try(lm(fmla_sc, data = transform(dsub, y_sc = y_sc)), silent = TRUE)
    if (inherits(fit_sc, "try-error")) return(.homosked())
    
    logvar_i <- as.numeric(fitted(fit_sc))
    logvar_0 <- as.numeric(predict(fit_sc, newdata = transform(newdata, y_sc = NA)))
    sigma_i  <- sqrt(pmax(exp(logvar_i), eps))
    sigma_0  <- sqrt(pmax(exp(logvar_0), eps))
    if (!is.finite(sigma_0) || any(!is.finite(sigma_i))) return(.homosked())
    
    r_adj <- r * (sigma_0 / sigma_i)
    if (verbose) {
      rsq <- summary(fit_sc)$r.squared
      cat("# twostage scale R2 =", round(rsq, 4), " sigma0 =", round(sigma_0, 6), "\n")
    }
    list(mu = mu, res = as.numeric(r_adj), sd = sigma_0, n = length(r_adj))
  }
  
  .dglm <- function() {
    if (!requireNamespace("dglm", quietly = TRUE)) stop("no dglm")
    fit <- dglm::dglm(fmla_mu,
                      dformula = as.formula(paste("~", cov_text)),
                      data = dsub, family = gaussian())
    mu     <- as.numeric(predict(fit, newdata = newdata, type = "response"))
    phi_i  <- as.numeric(predict(fit$dispersion.fit, type = "response"))
    phi_0  <- as.numeric(predict(fit$dispersion.fit, newdata = newdata, type = "response"))
    sigma_i <- sqrt(pmax(phi_i, 1e-12))
    sigma_0 <- sqrt(pmax(phi_0, 1e-12))
    res     <- as.numeric(residuals(fit, type = "response"))
    r_adj   <- res * (sigma_0 / sigma_i)
    if (verbose) {
      cat("# dglm sigma range =", paste(round(range(sigma_i), 6), collapse = " .. "),
          " sigma0 =", round(sigma_0, 6), "\n")
    }
    list(mu = mu, res = r_adj, sd = sigma_0, n = length(r_adj))
  }
  
  .gamlss <- function() {
    if (!requireNamespace("gamlss", quietly = TRUE) ||
        !requireNamespace("gamlss.dist", quietly = TRUE)) stop("no gamlss")
    fam <- gamlss.dist::NO()
    fit <- gamlss::gamlss(fmla_mu,
                          sigma.fo = if (interact_in_scale)
                            as.formula(paste("~", cov_text))
                          else
                            as.formula(paste("~", paste(covars, collapse = " + "))),
                          family = fam, data = dsub,
                          control = gamlss::gamlss.control(trace = FALSE))
    mu      <- as.numeric(gamlss::predict(fit, newdata = newdata, what = "mu",    type = "response"))
    sigma_0 <- as.numeric(gamlss::predict(fit, newdata = newdata, what = "sigma", type = "response"))
    mu_i    <- as.numeric(gamlss::fitted(fit, what = "mu"))
    sigma_i <- as.numeric(gamlss::fitted(fit, what = "sigma"))
    res     <- as.numeric(dsub[[outcome]] - mu_i)
    sigma_0 <- sqrt(pmax(sigma_0^2, 1e-12))
    sigma_i <- sqrt(pmax(sigma_i^2, 1e-12))
    r_adj   <- res * (sigma_0 / sigma_i)
    if (verbose) {
      cat("# gamlss sigma range =", paste(round(range(sigma_i), 6), collapse = " .. "),
          " sigma0 =", round(sigma_0, 6), "\n")
    }
    list(mu = mu, res = r_adj, sd = sigma_0, n = length(r_adj))
  }
  
  backends <- switch(backend,
                     auto = c("dglm","gamlss","twostage","homoskedastic"),
                     dglm = "dglm",
                     gamlss = "gamlss",
                     twostage = "twostage",
                     homoskedastic = "homoskedastic"
  )
  
  for (b in backends) {
    out <- try(switch(b,
                      dglm = .dglm(),
                      gamlss = .gamlss(),
                      twostage = .twostage(),
                      homoskedastic = .homosked()),
               silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    if (verbose) cat("# backend", b, "failed; trying next\n")
  }
  
  .homosked()
}


# empirical CDF with sensible fallback if too few residuals
cdf_empirical_or_norm <- function(residuals, sd_fallback) {
  rr <- residuals[is.finite(residuals)]
  if (length(rr) >= 10 && sd(rr, na.rm = TRUE) > 0) {
    ecdf(rr)
  } else {
    # fallback to normal with mean 0 and sd = sd_fallback (or 1 if 0)
    s <- if (is.finite(sd_fallback) && sd_fallback > 0) sd_fallback else 1
    function(z) pnorm(z, mean = 0, sd = s)
  }
}

# small helpers for safety
.safe_q <- function(muA, muB, Fref) {
  if (!is.null(Fref) && is.finite(muA) && is.finite(muB)) Fref(muA - muB) else NA_real_
}
.safe_z <- function(muA, muB, sdref) {
  if (is.finite(muA) && is.finite(muB) && is.finite(sdref) && sdref > 0) (muA - muB)/sdref else NA_real_
}

make_norm_cdf <- function(s) {
  if (is.finite(s) && s > 0) function(z) pnorm(z, mean = 0, sd = s) else NULL
}

z_to_unit <- function(z) {
  # map z in [-zr, zr] to [0,1], clipping outside
  pmin(1, pmax(0, 0.5 + 0.5 * (z / zr)))
}

invert_phantom <- function(expr) {
  if (inherits(expr, "expression")) {
    expr[[1]] <- invert_phantom(expr[[1]])
    expr
  } else if (is.call(expr) && expr[[1]] == as.name("phantom")) {
    expr[[2]]
  } else if (is.atomic(expr)) {
    substitute(phantom(TEXT), list(TEXT = expr))
  } else if (is.name(expr)) {
    expr
  } else if (is.call(expr)) {
    new_args <- lapply(expr[-1], invert_phantom)
    as.call(c(expr[[1]], new_args))
  } else {
    invert_phantom(as.list(expr))
  }
}

interp_x_for_y <- function(y_target, x1, y1, x2, y2) {
  if (y1 == y2) {
    return((x1 + x2) / 2)
  }
  t <- (y_target - y1) / (y2 - y1)
  x1 + t * (x2 - x1)
}
