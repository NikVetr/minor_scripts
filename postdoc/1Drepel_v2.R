#minsep_fit is the clear winner here
test_methods <- F

if(test_methods){
  
  adjust_spacing <- function(x, w, bounds = NULL, method = "greedy",
                             w_sep = 1e6, w_sim = 1) {
    
    # Validate inputs
    if (!is.numeric(x) || !is.numeric(w)) {
      stop("x and w must be numeric")
    }
    if (length(w) != 1) {
      stop("w must be a scalar")
    }
    if (!is.null(bounds)) {
      if (length(bounds) != 2 || !is.numeric(bounds)) {
        stop("bounds must be NULL or a numeric vector of length 2")
      }
      if (bounds[2] < bounds[1]) {
        stop("Upper bound must be >= lower bound")
      }
    }
    
    n <- length(x)
    if (n == 0) return(numeric(0))
    if (n == 1) {
      if (!is.null(bounds)) {
        return(pmax(bounds[1], pmin(bounds[2], x)))
      }
      return(x)
    }
    
    min_spacing <- w
    
    # Check if problem is feasible with bounds
    if (!is.null(bounds)) {
      range_available <- bounds[2] - bounds[1]
      range_needed <- (n - 1) * min_spacing
      if (range_needed > range_available) {
        warning("Cannot fit all points with required spacing within bounds. ",
                "Spacing will be compressed.")
      }
    }
    
    if (method == "greedy") {
      # Greedy two-pass algorithm (fast and effective)
      
      # Sort x and track original order
      ord <- order(x)
      x_sorted <- x[ord]
      
      # Initialize y with sorted x
      y_sorted <- x_sorted
      
      # Apply bounds to initial positions if provided
      if (!is.null(bounds)) {
        y_sorted <- pmax(bounds[1], pmin(bounds[2], y_sorted))
      }
      
      # Forward pass: ensure each point is at least min_spacing from previous
      for (i in 2:n) {
        min_pos <- y_sorted[i - 1] + min_spacing
        if (y_sorted[i] < min_pos) {
          y_sorted[i] <- min_pos
        }
        # Respect upper bound
        if (!is.null(bounds) && y_sorted[i] > bounds[2]) {
          y_sorted[i] <- bounds[2]
        }
      }
      
      # Backward pass: pull points back toward originals while maintaining spacing
      for (i in (n - 1):1) {
        max_pos <- y_sorted[i + 1] - min_spacing
        if (y_sorted[i] > max_pos) {
          y_sorted[i] <- max_pos
        }
        # Try to get closer to original if possible
        if (y_sorted[i] > x_sorted[i]) {
          y_sorted[i] <- max(x_sorted[i], max_pos)
        }
        # Respect lower bound
        if (!is.null(bounds) && y_sorted[i] < bounds[1]) {
          y_sorted[i] <- bounds[1]
        }
      }
      
      # Restore original order
      y <- numeric(n)
      y[ord] <- y_sorted
      
    } else if (method == "optim") {
      # Optimization-based approach
      
      ord <- order(x)
      x_sorted <- x[ord]
      
      # Parameterize as displacements from minimum possible positions
      # Start position
      start_pos <- if (!is.null(bounds)) bounds[1] else min(x_sorted)
      
      # Objective: displacements[i] represents position[i] - (start_pos + i * min_spacing)
      # This ensures monotonicity and makes optimization easier
      
      # Loss function
      loss_fn <- function(displacements) {
        positions <- start_pos + (0:(n-1)) * min_spacing + displacements
        
        # Apply bounds penalty
        bound_penalty <- 0
        if (!is.null(bounds)) {
          bound_penalty <- sum(pmax(0, bounds[1] - positions)^2) + 
            sum(pmax(0, positions - bounds[2])^2)
        }
        
        # Similarity penalty: want positions close to x_sorted
        sim_penalty <- sum((positions - x_sorted)^2)
        
        # Separation penalty: penalize negative displacements relative to required
        # (displacements should be >= 0 to maintain spacing)
        sep_penalty <- sum(pmax(0, -displacements)^2)
        
        return(w_sim * sim_penalty + w_sep * sep_penalty + 1e8 * bound_penalty)
      }
      
      # Initial guess: try to match x_sorted
      init_disp <- pmax(0, x_sorted - (start_pos + (0:(n-1)) * min_spacing))
      
      # Optimize
      result <- optim(par = init_disp, fn = loss_fn, method = "L-BFGS-B",
                      lower = rep(0, n))
      
      y_sorted <- start_pos + (0:(n-1)) * min_spacing + result$par
      
      # Apply bounds
      if (!is.null(bounds)) {
        y_sorted <- pmax(bounds[1], pmin(bounds[2], y_sorted))
      }
      
      # Restore original order
      y <- numeric(n)
      y[ord] <- y_sorted
      
    } else {
      stop("method must be 'greedy' or 'optim'")
    }
    
    return(y)
  }
  
  
  space_out_vector_optim <- function(x, w, bounds = NULL, w_sim = 1.0, 
                                     w_sep = 1.0, initial_y = NULL) {
    if (!is.numeric(x) || !is.numeric(w)) {
      stop("Both 'x' and 'w' must be numeric.")
    }
    n <- length(x)
    if (n <= 1) {
      return(x)
    }
    
    min_sep <- w
    
    # We will optimize the displacements from the original sorted positions.
    original_order <- order(x)
    x_sorted <- x[original_order]
    
    # Loss function to be minimized by optim()
    loss_function <- function(y) {
      # 1. Similarity loss: squared distance from original positions
      loss_sim <- sum((y - x_sorted)^2)
      
      # 2. Separation loss: penalize points that are too close
      displacements <- diff(y)
      # We only care about displacements that are less than the minimum separation
      violations <- min_sep - displacements
      # Penalize only the violations (where violations > 0)
      # Using a smooth penalty function (e.g., squared)
      loss_sep <- sum(ifelse(violations > 0, violations^2, 0))
      
      # Total weighted loss
      return(w_sim * loss_sim + 10 * w_sep * loss_sep)
    }
    
    # Use the sorted x as the initial guess for y
    if(is.null(initial_y)){
      initial_y <- x_sorted  
    }
    
    # Set up optimization parameters
    # L-BFGS-B method allows for box constraints (bounds)
    lower_bounds <- if (!is.null(bounds)) rep(bounds[1], n) else -Inf
    upper_bounds <- if (!is.null(bounds)) rep(bounds[2], n) else Inf
    
    # Run the optimization
    # We pass initial_y as the starting point for the 'par'ameters (our y vector)
    opt_result <- optim(
      par = initial_y,
      fn = loss_function,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds
    )
    
    # Extract the optimized y vector
    y_optimized <- opt_result$par
    
    # Restore original order
    y_final <- y_optimized[order(original_order)]
    
    return(y_final)
  }
  
}

# project x onto vectors with min pairwise spacing w/2, as close (L2) as possible
# x: numeric vector
# w: numeric scalar (>= 0)
# bounds: NULL or numeric length-2 c(lower, upper)
# verbose: logical; print debugging info
minsep_fit <- function(x, w, bounds = NULL, verbose = FALSE) {
  # input checks
  if (!is.numeric(x) || !is.vector(x)) stop("x must be a numeric vector")
  if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w < 0) stop("w must be a finite, nonnegative scalar")
  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) stop("bounds must be c(lower, upper) with finite numbers")
    if (bounds[1] > bounds[2]) stop("bounds[1] must be <= bounds[2]")
  }
  
  n <- length(x)
  if (n == 0L) return(x)
  s <- w
  
  # feasibility check for bounds + spacing
  if (!is.null(bounds)) {
    L <- bounds[1]; U <- bounds[2]
    if ((U - L) < s * (n - 1)) {
      stop("infeasible: upper - lower < (w/2) * (n - 1); cannot fit all points with required spacing within bounds")
    }
  }
  
  # sort x (stable tie-break by index), remember permutation
  idx <- order(x, seq_along(x))
  xs <- x[idx]
  
  # shift to isotonic form: z_i = y_i - i*s should be nondecreasing
  i_seq <- seq_len(n)
  x_shift <- xs - s * i_seq
  
  # per-index bounds for the shifted space
  if (is.null(bounds)) {
    lower <- rep(-Inf, n)
    upper <- rep( Inf, n)
  } else {
    lower <- rep(bounds[1], n) - s * i_seq
    upper <- rep(bounds[2], n) - s * i_seq
  }
  
  # bounded isotonic regression via pooled blocks
  # returns z (nondecreasing) minimizing sum (z - x_shift)^2 subject to lower <= z <= upper
  bounded_isoreg <- function(xv, lo, hi) {
    n <- length(xv)
    if (length(lo) != n || length(hi) != n) stop("bounds must match length of x")
    if (any(lo > hi)) stop("elementwise lower > upper encountered")
    # initialize one-element blocks
    b_start <- seq_len(n)
    b_end   <- seq_len(n)
    b_sum   <- xv
    b_w     <- rep(1.0, n)
    b_lo    <- lo
    b_hi    <- hi
    # block value is the block mean, clipped to the intersection of bounds
    b_val   <- pmin(pmax(b_sum / b_w, b_lo), b_hi)
    
    m <- n
    i <- 1L
    # pav with bounds: merge while violating nondecreasing order
    while (i < m) {
      if (b_val[i] > b_val[i + 1L]) {
        # merge block i with i+1
        b_end[i] <- b_end[i + 1L]
        b_sum[i] <- b_sum[i] + b_sum[i + 1L]
        b_w[i]   <- b_w[i]   + b_w[i + 1L]
        b_lo[i]  <- max(b_lo[i], b_lo[i + 1L])
        b_hi[i]  <- min(b_hi[i], b_hi[i + 1L])
        if (b_lo[i] > b_hi[i]) stop("infeasible constraints while merging blocks (check global feasibility)")
        # recompute value, then remove block i+1
        mu <- b_sum[i] / b_w[i]
        if (mu < b_lo[i]) mu <- b_lo[i]
        if (mu > b_hi[i]) mu <- b_hi[i]
        b_val[i] <- mu
        
        # drop block i+1
        if (i + 1L < m) {
          b_start <- b_start[-(i + 1L)]
          b_end   <- b_end[-(i + 1L)]
          b_sum   <- b_sum[-(i + 1L)]
          b_w     <- b_w[-(i + 1L)]
          b_lo    <- b_lo[-(i + 1L)]
          b_hi    <- b_hi[-(i + 1L)]
          b_val   <- b_val[-(i + 1L)]
        } else {
          b_start <- b_start[-(i + 1L)]
          b_end   <- b_end[-(i + 1L)]
          b_sum   <- b_sum[-(i + 1L)]
          b_w     <- b_w[-(i + 1L)]
          b_lo    <- b_lo[-(i + 1L)]
          b_hi    <- b_hi[-(i + 1L)]
          b_val   <- b_val[-(i + 1L)]
        }
        m <- m - 1L
        # after merging, step back if needed
        if (i > 1L) {
          i <- i - 1L
        }
      } else {
        i <- i + 1L
      }
    }
    
    # expand block values back to length n
    z <- numeric(n)
    for (k in seq_len(m)) {
      z[b_start[k]:b_end[k]] <- b_val[k]
    }
    z
  }
  
  if (verbose) {
    cat("# n:", n, "\n")
    cat("# s (min gap):", s, "\n")
    if (!is.null(bounds)) cat("# bounds on y:", bounds[1], bounds[2], "\n")
  }
  
  z <- bounded_isoreg(x_shift, lower, upper)
  y_sorted <- z + s * i_seq
  
  # restore original order
  y <- numeric(n)
  y[idx] <- y_sorted
  
  if (verbose) {
    diffs_sorted <- diff(y_sorted)
    min_gap <- if (length(diffs_sorted)) min(diffs_sorted) else Inf
    cat("# min adjacent gap in sorted y:", min_gap, "\n")
    if (!is.null(bounds)) {
      cat("# min(y), max(y):", min(y), max(y), "\n")
    }
    cat("# mean squared deviation:", mean((y - x)^2), "\n")
  }
  
  y
}

# project x onto vectors with min pairwise spacing w (adjacent gap), as close (L2) as possible
# x: numeric vector
# w: numeric scalar (>= 0)
# bounds: NULL or numeric length-2 c(lower, upper)
# verbose: logical; print debugging info
# resize: logical; if TRUE and infeasible within bounds, shrink w just enough to fit exactly
minsep_fit <- function(x, w, bounds = NULL, verbose = FALSE, resize = FALSE) {
  # input checks
  if (!is.numeric(x) || !is.vector(x)) stop("x must be a numeric vector")
  if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w < 0) stop("w must be a finite, nonnegative scalar")
  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) stop("bounds must be c(lower, upper) with finite numbers")
    if (bounds[1] > bounds[2]) stop("bounds[1] must be <= bounds[2]")
  }
  
  n <- length(x)
  if (n == 0L) return(x)
  s <- w
  scale_factor <- NULL
  
  # numeric tolerance used for feasibility and bound intersections
  rtol <- 64 * .Machine$double.eps
  
  # feasibility check for bounds + spacing; optionally shrink spacing to fit exactly
  if (!is.null(bounds)) {
    L <- bounds[1]; U <- bounds[2]
    have <- U - L
    need <- s * (n - 1)
    
    if (resize) {
      # cap s to the largest feasible gap; if already <= cap (within tol), leave s alone
      gap_cap <- if (n > 1L) have / (n - 1L) else Inf
      # shrink only if s exceeds cap by more than tolerance
      if (is.finite(gap_cap) && s > gap_cap * (1 + rtol)) {
        scale_factor <- if (s > 0) gap_cap / s else 1
        s <- gap_cap
        if (verbose) cat("# resized to fit: scale_factor:", scale_factor, " new s:", s, "\n")
      } else if (have + rtol * max(1, abs(have), abs(need)) < need) {
        stop("infeasible: upper - lower < w * (n - 1); cannot fit all points with required spacing within bounds")
      }
    } else {
      if (have + rtol * max(1, abs(have), abs(need)) < need) {
        stop("infeasible: upper - lower < w * (n - 1); cannot fit all points with required spacing within bounds")
      }
    }
  }
  
  # sort x (stable tie-break by index), remember permutation
  idx <- order(x, seq_along(x))
  xs <- x[idx]
  
  # shift to isotonic form: z_i = y_i - i*s should be nondecreasing
  i_seq <- seq_len(n)
  x_shift <- xs - s * i_seq
  
  # per-index bounds for the shifted space
  if (is.null(bounds)) {
    lower <- rep(-Inf, n)
    upper <- rep( Inf, n)
  } else {
    lower <- rep(bounds[1], n) - s * i_seq
    upper <- rep(bounds[2], n) - s * i_seq
  }
  
  # bounded isotonic regression via pooled blocks
  bounded_isoreg <- function(xv, lo, hi) {
    n <- length(xv)
    if (length(lo) != n || length(hi) != n) stop("bounds must match length of x")
    if (any(lo > hi & (lo - hi) > rtol * pmax(1, abs(lo), abs(hi)))) stop("elementwise lower > upper encountered")
    # initialize one-element blocks
    b_start <- seq_len(n)
    b_end   <- seq_len(n)
    b_sum   <- xv
    b_w     <- rep(1.0, n)
    b_lo    <- lo
    b_hi    <- hi
    # block value is the block mean, clipped to the intersection of bounds
    b_val   <- pmin(pmax(b_sum / b_w, b_lo), b_hi)
    
    m <- n
    i <- 1L
    # pav with bounds: merge while violating nondecreasing order
    while (i < m) {
      if (b_val[i] > b_val[i + 1L]) {
        # merge block i with i+1
        b_end[i] <- b_end[i + 1L]
        b_sum[i] <- b_sum[i] + b_sum[i + 1L]
        b_w[i]   <- b_w[i]   + b_w[i + 1L]
        # intersect bounds and tolerate tiny inversions
        new_lo <- max(b_lo[i], b_lo[i + 1L])
        new_hi <- min(b_hi[i], b_hi[i + 1L])
        if (new_lo > new_hi) {
          if (new_lo - new_hi <= rtol * max(1, abs(new_lo), abs(new_hi))) {
            mid <- (new_lo + new_hi) / 2
            new_lo <- mid
            new_hi <- mid
          } else {
            stop("infeasible constraints while merging blocks (check global feasibility)")
          }
        }
        b_lo[i]  <- new_lo
        b_hi[i]  <- new_hi
        
        # recompute value, then remove block i+1
        mu <- b_sum[i] / b_w[i]
        if (mu < b_lo[i]) mu <- b_lo[i]
        if (mu > b_hi[i]) mu <- b_hi[i]
        b_val[i] <- mu
        
        # drop block i+1
        if (i + 1L < m) {
          b_start <- b_start[-(i + 1L)]
          b_end   <- b_end[-(i + 1L)]
          b_sum   <- b_sum[-(i + 1L)]
          b_w     <- b_w[-(i + 1L)]
          b_lo    <- b_lo[-(i + 1L)]
          b_hi    <- b_hi[-(i + 1L)]
          b_val   <- b_val[-(i + 1L)]
        } else {
          b_start <- b_start[-(i + 1L)]
          b_end   <- b_end[-(i + 1L)]
          b_sum   <- b_sum[-(i + 1L)]
          b_w     <- b_w[-(i + 1L)]
          b_lo    <- b_lo[-(i + 1L)]
          b_hi    <- b_hi[-(i + 1L)]
          b_val   <- b_val[-(i + 1L)]
        }
        m <- m - 1L
        if (i > 1L) i <- i - 1L
      } else {
        i <- i + 1L
      }
    }
    
    # expand block values back to length n
    z <- numeric(n)
    for (k in seq_len(m)) {
      z[b_start[k]:b_end[k]] <- b_val[k]
    }
    z
  }
  
  if (verbose) {
    cat("# n:", n, "\n")
    cat("# s (min gap):", s, "\n")
    if (!is.null(bounds)) cat("# bounds on y:", bounds[1], bounds[2], "\n")
  }
  
  z <- bounded_isoreg(x_shift, lower, upper)
  y_sorted <- z + s * i_seq
  
  # restore original order
  y <- numeric(n)
  y[idx] <- y_sorted
  
  if (verbose) {
    diffs_sorted <- diff(y_sorted)
    min_gap <- if (length(diffs_sorted)) min(diffs_sorted) else Inf
    cat("# min adjacent gap in sorted y:", min_gap, "\n")
    if (!is.null(bounds)) {
      cat("# min(y), max(y):", min(y), max(y), "\n")
    }
    cat("# mean squared deviation:", mean((y - x)^2), "\n")
  }
  
  if (!is.null(scale_factor)) {
    return(list(y = y, scale_factor = scale_factor))
  }
  y
}


if(test_methods){
  n <- 50
  x <- rnorm(n)
  w <- 0.5
  
  #optim approach
  y1 <- space_out_vector_optim(x, w)
  y1 <- space_out_vector_optim(x, w, w_sim = 1, w_sep = 10,
                               initial_y = y1)
  
  #isotonic regression approach
  y2 <- minsep_fit(x, w)
  
  #greedy algorithm approach
  y3 <- adjust_spacing(x, w)
  
  plot(NULL, ylim = c(-0.1,1.1), xlim = range(c(x,y1, y2, y3)), 
       xlab = "", ylab = "")
  
  points(x, rep(0,n))
  points(y1, rep(0.15,n))
  segments(x,0,y1,0.15)
  
  points(x, rep(0.25,n))
  points(y2, rep(0.4,n))
  segments(x,0.25,y2,0.4)
  
  
  points(x, rep(0.5,n))
  points(y3, rep(0.65,n))
  segments(x,0.5,y3,0.65)
  
  #how well did we hit the spacing?
  summary(diff(sort(y1)) / (w/2))
  summary(diff(sort(y2)) / (w/2))
  summary(diff(sort(y3)) / (w/2))
  
  #how close are we to the original marks?
  mean(abs(y1-x) / (w/2))
  mean(abs(y2-x) / (w/2))
  mean(abs(y3-x) / (w/2))
}

bezier_curve <- function(x0, y0, x1, y1,
                         p = 0.5, n = 128, k = 1,
                         ends = c("flat", "steep"),
                         col = 1, return_points = FALSE, debug = FALSE, ...) {
  # choose behavior: "flat" = horizontal tangents at ends (original);
  #                  "steep" = vertical tangents at ends (the “other direction”)
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

detect_blocks_igraph <- function(mat_cols, cor_thresh = 0, seed = 1L) {
  # mat_cols: symmetric similarity/correlation matrix with row/colnames
  # cor_thresh: drop edges with |cor| < cor_thresh (0 = no threshold)
  
  if(is.null(rownames(mat_cols)) || is.null(colnames(mat_cols))) {
    stop("mat_cols must have row and column names")
  }
  
  if(!all(rownames(mat_cols) == colnames(mat_cols))) {
    stop("row and column names of mat_cols must match and be in the same order")
  }
  
  w <- mat_cols
  
  # remove self-edges
  diag(w) <- 0
  
  # optional thresholding to sharpen community structure
  if(!is.null(cor_thresh) && !is.na(cor_thresh)) {
    w[abs(w) < cor_thresh] <- 0
  }
  
  # build graph
  g <- igraph::graph_from_adjacency_matrix(
    w,
    mode   = "undirected",
    weighted = TRUE,
    diag   = FALSE
  )
  
  set.seed(seed)
  comm <- igraph::cluster_louvain(g)
  
  memb <- igraph::membership(comm)
  
  # ensure the membership is ordered by vertex name
  memb <- memb[igraph::V(g)$name]
  
  # sanity: names should match rownames(mat_cols)
  if(!all(names(memb) == rownames(mat_cols))) {
    stop("membership names do not match matrix rownames")
  }
  
  return(memb)   # named integer vector: variable -> block id
}
