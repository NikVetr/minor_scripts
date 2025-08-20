# Load necessary libraries
library(expm)
library(optimx)
library(pracma) # For cross product

# --- Helper Functions ---
# Helper function to vectorize the upper triangle of a skew-symmetric matrix
vectorize_skew_symmetric <- function(A) {
  k <- nrow(A)
  if (k < 2) return(numeric(0))
  p <- A[upper.tri(A)]
  return(p)
}

# Helper function to reconstruct a k x k skew-symmetric matrix from its vectorized form
reconstruct_skew_symmetric <- function(p, k) {
  expected_len <- k * (k - 1) / 2
  if (length(p) != expected_len) {
    stop(paste("Incorrect number of parameters. Expected", expected_len, "but got", length(p), "for dimension k =", k))
  }
  A <- matrix(0, k, k)
  if (k < 2) return(A)
  A[upper.tri(A)] <- p
  A <- A - t(A)
  return(A)
}

# Generate a random skew-symmetric matrix in k dimensions
random_skew_symmetric <- function(k) {
  M <- matrix(rnorm(k * k), nrow = k)
  A <- M - t(M)
  return(A)
}

skew_to_rot <- function(A){
  k <- dim(A)[1]
  R <- expm(A)
  svd_res <- svd(R)
  R <- svd_res$u %*% diag(c(rep(1, k-1), det(svd_res$u %*% t(svd_res$v)))) %*% t(svd_res$v)
  return(R)
}

# Generate a random rotation matrix in SO(k)
random_rotation_matrix <- function(k) {
  A <- random_skew_symmetric(k)
  return(skew_to_rot(A))
}

rot_to_skew <- function(R){
  logR <- expm::logm(R)
  skew <- 0.5 * (logR - t(logR))
  return(skew)
}

# Calculate magnitude of angular velocity for a rotation R over time delta_t
angular_velocity <- function(R, delta_t = 1) {
  k <- nrow(R)
  svd_res <- svd(R)
  det_uv <- det(svd_res$u %*% t(svd_res$v))
  diag_mat <- diag(c(rep(1, k - 1), det_uv))
  R_proj <- svd_res$u %*% diag_mat %*% t(svd_res$v)
  
  if (norm(R_proj - diag(k), type = "F") < sqrt(.Machine$double.eps) * k) {
    return(0)
  }
  logm_R <- try(logm(R_proj), silent = TRUE)
  if(inherits(logm_R, "try-error")){
    warning("logm failed in angular_velocity")
    # What to return? 0 might bias variance down, NA needs handling. Let's try 0.
    return(0) 
  }
  A_skew <- 0.5 * (logm_R - t(logm_R))
  omega_mag <- norm(A_skew, type = "F") / sqrt(2)
  return(max(0, omega_mag / delta_t))
}

# Log-Exp interpolation between R1 and R2 with weight w
blend_skew_symmetric <- function(A1, A2, w, method = c("lie_group", "linear")) {
  method <- match.arg(method)
  
  if (!all(dim(A1) == c(3, 3)) || !all(dim(A2) == c(3, 3))) {
    stop("Both A1 and A2 must be 3x3 matrices.")
  }
  if (abs(w) > 1) stop("w should be between 0 and 1.")
  
  if (method == "linear") {
    return((1 - w) * A1 + w * A2)
  }
  
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("Please install the 'expm' package: install.packages('expm')")
  }
  
  R1 <- expm::expm(A1)
  R2 <- expm::expm(A2)
  
  # Interpolated rotation matrix
  R_blend <- R1 %*% expm::expm(w * expm::logm(t(R1) %*% R2))
  
  # Return the skew-symmetric matrix that generates this
  A_blend <- 0.5 * (expm::logm(R_blend) - t(expm::logm(R_blend)))
  return(A_blend)
}


# Blending weights function (logistic or linear)
blending_weights <- function(n, shape = 0) {
  stopifnot(n >= 1) # Can have n=1
  if (n == 1) return(0) 
  t <- seq(0, 1, length.out = n)
  if (is.na(shape) || shape == 0) {
    w <- t
  } else {
    shape <- max(-5, min(5, shape))
    a <- 10^shape
    raw_w <- 1 / (1 + exp(-a * (t - 0.5)))
    min_w <- raw_w[1]
    max_w <- raw_w[n]
    if (abs(max_w - min_w) < .Machine$double.eps) {
      w <- t
    } else {
      w <- (raw_w - min_w) / (max_w - min_w)
    }
  }
  w[1] <- 0
  w[n] <- 1 # Ensure ends are exact
  return(w)
}

blending_simplex_path <- function(k, n, shape = 2, pow_ends = 2) {
  
  # Steepness of the logistic transitions
  a <- 10^shape
  
  # Transition points between dimensions
  t_grid <- seq(0, 1, length.out = n)
  t_breaks <- seq(0, 1, length.out = k + 1)
  
  # Logistic sigmoid
  sigmoid <- function(x) 1 / (1 + exp(-x))
  
  # Initialize output matrix
  weights <- matrix(0, nrow = n, ncol = k)
  
  for (i in seq_len(k)) {
    # For each dimension, compute its logistic "bump"
    lower <- sigmoid(a * (t_grid - t_breaks[i]))
    upper <- sigmoid(a * (t_grid - t_breaks[i + 1]))
    weights[, i] <- lower - upper
  }
  
  #raise weights to a power
  weights_pow <- abs(t_grid * 2 - 1)^pow_ends * (max(0, pow_ends-1)) + 1
  weights <- weights^weights_pow
  
  # Normalize rows to sum to 1
  row_sums <- rowSums(weights)
  weights <- weights / row_sums
  
  return(weights)
}


#test out simplex function
test_simplex_func <- F
if(test_simplex_func){
  n_simplex <- 2
  n_pts <- 100
  shape_simplex <- 0.8
  pow_ends <- 2
  simplex_weights <- blending_simplex_path(n_simplex, n_pts, shape_simplex, pow_ends)
  simplex_cumweights <- t(apply(simplex_weights, 1, cumsum))
  
  par(mfrow = c(2,1), mar = c(4,4,2,2))
  plot(x = 1:n_pts, simplex_weights[,1], type = "l", ylim = c(0,1))
  for(i in 2:n_simplex){
    lines(x = 1:n_pts, simplex_weights[,i], col = i)  
  }
  
  simplex_cumweights <- cbind(0, simplex_cumweights)
  plot(x = 1:n_pts, simplex_cumweights[,1], type = "l", ylim = c(0,1), col = 0)
  for(i in 2:(n_simplex+1)){
    polygon(x = c(1:n_pts, n_pts:1), c(simplex_cumweights[,i], rev(simplex_cumweights[,i-1])), col = i-1)  
  }
}


# Geodesic distance on SO(k)
rotation_distance <- function(R1, R2) {
  k <- nrow(R1)
  R_rel <- try(t(R1) %*% R2, silent=TRUE)
  if(inherits(R_rel, "try-error")){
    warning("Matrix multiplication failed in rotation_distance (t(R1) %*% R2)")
    return(pi) # Return a large distance like pi?
  }
  
  svd_res <- svd(R_rel)
  det_uv <- det(svd_res$u %*% t(svd_res$v))
  diag_mat <- diag(c(rep(1, k - 1), det_uv))
  R_rel_proj <- svd_res$u %*% diag_mat %*% t(svd_res$v)
  
  if (norm(R_rel_proj - diag(k), type = "F") < sqrt(.Machine$double.eps) * k) {
    return(0)
  }
  
  A <- try(logm(R_rel_proj), silent=TRUE)
  if(inherits(A, "try-error")){
    warning("logm failed in rotation_distance")
    return(pi) # Return a large distance
  }
  A_skew <- 0.5 * (A - t(A))
  d <- norm(A_skew, type = "F") / sqrt(2)
  return(max(0, d)) # Ensure non-negative
}

# dkrot <- function(k, thetas) {
#   planes <- t(combn(k, 2))  # same ordering as in `dkrot()`
#   planes <- planes[order(apply(planes, 1, max)),]
#   angles <- data.frame(cbind(planes, theta = thetas))
#   rotmats <- lapply(angles$theta, d2rot)
#   big_rotmats <- lapply(1:choose(k, 2), function(i){
#     bigrm <- diag(k)
#     bigrm[angles$V1[i], angles$V1[i]] <- rotmats[[i]][1,1]
#     bigrm[angles$V1[i], angles$V2[i]] <- rotmats[[i]][1,2]
#     bigrm[angles$V2[i], angles$V1[i]] <- rotmats[[i]][2,1]
#     bigrm[angles$V2[i], angles$V2[i]] <- rotmats[[i]][2,2]
#     bigrm
#   })
#   rotmat <- Reduce("%*%", big_rotmats)
#   rotmat
# }

compute_rotmat <- function(angles, rotmat, k){
  bigrm <- diag(k)
  bigrm[angles$V1, angles$V1] <- rotmat[1,1]
  bigrm[angles$V1, angles$V2] <- rotmat[1,2]
  bigrm[angles$V2, angles$V1] <- rotmat[2,1]
  bigrm[angles$V2, angles$V2] <- rotmat[2,2]
  bigrm
}

compute_rotmat_theta <- function(angles, theta, k){
  bigrm <- diag(k)
  bigrm[angles$X1, angles$X1] <- cos(theta)
  bigrm[angles$X1, angles$X2] <- -sin(theta)
  bigrm[angles$X2, angles$X1] <- sin(theta)
  bigrm[angles$X2, angles$X2] <- cos(theta)
  bigrm
}

dkrot <- function(k, thetas) {
  planes <- t(combn(k, 2))  # same ordering as in `dkrot()`
  planes <- data.frame(planes[order(apply(planes, 1, max)),])
  big_rotmats <- lapply(1:choose(k, 2), function(i){
    compute_rotmat_theta(planes[i,], thetas[i], k)
  })
  rotmat <- Reduce("%*%", big_rotmats)
  rotmat
}

undokrot <- function(rotmat) {
  k <- nrow(rotmat)
  # list of plane pairs in the same order as dkrot
  planes <- t(combn(k, 2))  # same ordering as in `dkrot()`
  planes <- planes[order(apply(planes, 1, max)),]
  N <- nrow(planes)
  
  # We will have: rotmat = G_1 * G_2 * ... * G_N
  # We want to peel off G_1, then G_2, etc., from the LEFT side in that order.
  
  # Start with the matrix we want to factor
  R <- rotmat
  thetas <- numeric(N)
  
  for (i in seq_len(N)) {
    i1 <- planes[i,1]
    i2 <- planes[i,2]
    
    # Solve for G_i from the left factor:
    #   R = G_i * M   ==>   M = G_i^{-1} * R = G_i^T * R
    # But we first need to figure out the angle used by G_i in the original (i1,i2) plane.
    
    # G_i^T is also a rotation in the same (i1,i2) plane, but we have to measure it correctly
    # from R. The sub-block R[i1, c(i1,i2); i2, c(i1,i2)] is already "contaminated" by G_2..G_N.
    #
    # There's no trivial one-liner to do this. We must find the angle that un-rotates R's columns
    # i1 and i2 back to something orthonormal in the original plane. The safest approach is to
    # re-embed the basis as we go.
    
    # One approach: look at columns i1 and i2 of R to see how the standard basis e_i1, e_i2
    # got rotated. Then deduce the rotation angle from those 2 vectors.
    #
    # Let v1 = R[,i1], v2 = R[,i2]. In an "ideal" scenario (kD rotation), these should be orthonormal
    # if R is orthonormal. The angle in plane (i1,i2) is the angle between the projections
    # of v1, v2 onto that plane and the original axes. This is complicated to do in a short snippet.
    
    # Pseudocode approach (not fully fleshed out):
    #
    # 1) let e1 = unit vector with 1 in position i1, 0 in others
    #    let e2 = unit vector with 1 in position i2, 0 in others
    # 2) find out how those got mapped by R => r1 = R %*% e1, r2 = R %*% e2
    # 3) in principle, G_i was the rotation that took e1->r1', e2->r2' before
    #    the next factors changed them further. So we have to solve or guess from geometry.
    #
    # Because so many frames have changed, you'd have to do a partial backward approach
    # or store partial products.
    
    # => In practice, doing the exact left-to-right peeling to get the same angles is quite complex.
    #    Usually you store the partial products if you need the same angles.
    
    # (for illustration, let's pretend we do a "standard" Givens extraction, but
    #  it won't produce the same angles as the original.)
    
    # Here is the typical "standard" approach that zeroes out R[i2,i1]:
    r  <- sqrt(R[i1,i1]^2 + R[i2,i1]^2)
    c  <- R[i1,i1]/r
    s  <- R[i2,i1]/r
    theta <- atan2(s, c)
    thetas[i] <- theta
    
    # Left-multiply R by G_i^T to peel it off
    GiT <- diag(k)
    GiT[i1,i1] <- c
    GiT[i1,i2] <- s
    GiT[i2,i1] <- -s
    GiT[i2,i2] <- c
    
    R <- GiT %*% R  # now we've peeled off G_i, conceptually
  }
  
  # 'thetas' are the angles from a standard row-by-row Givens factorization approach
  # with the same plane order, but it still won't match the original angles from 'dkrot()'.
  
  return(thetas)
}

blendRotations <- function(thetas1, thetas2, w = 0.5,
                           method = c("linear", "geodesic")) {
  method <- match.arg(method)
  
  # Sanity checks
  if (length(thetas1) != length(thetas2)) {
    stop("theta1 and theta2 must have the same length.")
  }
  if (w < 0 || w > 1) {
    stop("w must be in [0,1].")
  }
  
  # Infer k from length(thetas1) = choose(k,2)
  # Solve k*(k-1)/2 = length(thetas1)
  # E.g. find integer k that satisfies the equation
  nTheta <- length(thetas1)
  # Not the most robust but works for typical small k
  k <- which( (seq_len(100)*(seq_len(100)-1))/2 == nTheta )
  if (length(k) != 1) {
    stop("Could not infer a unique k from length(thetas1). Please specify or check your inputs.")
  }
  
  if (method == "linear") {
    # Naive angle-wise interpolation
    theta_blend <- (1 - w)*thetas1 + w*thetas2
    Rblend <- dkrot(k, theta_blend)
    return(Rblend)
  } else if (method == "geodesic") {
    # Shortest path on SO(k)
    # R(t) = R1 * exp(t * log(R1^T * R2)), for t in [0,1]
    
    # 1) Construct the two rotation matrices
    R1 <- dkrot(k, thetas1)
    R2 <- dkrot(k, thetas2)
    
    # 2) Use matrix log/expm for the geodesic
    #    Rblend = R1 * exp( w * log( R1^T * R2 ) )
    library(expm)         # for expm(), logm()
    M <- logm(t(R1) %*% R2)
    Rblend <- R1 %*% expm(w * M)
    
    return(Rblend)
  }
}

k <- 5
thetas <- runif(choose(k, 2), 0, 2*pi)
rotmat <- dkrot(k, thetas)
thetas_recovered <- undokrot(rotmat)
dkrot(k, thetas_recovered) - rotmat

blend_givens <- function(thetas1, thetas2, w = 0.5,
                         method = c("linear", "geodesic")) {
  method <- match.arg(method)
  
  # Sanity checks
  if (length(thetas1) != length(thetas2)) {
    stop("theta1 and theta2 must have the same length.")
  }
  if (w < 0 || w > 1) {
    stop("w must be in [0,1].")
  }
  
  # Infer k from length(thetas1) = choose(k,2)
  # Solve k*(k-1)/2 = length(thetas1)
  # E.g. find integer k that satisfies the equation
  nTheta <- length(thetas1)
  # Not the most robust but works for typical small k
  k <- which( (seq_len(100)*(seq_len(100)-1))/2 == nTheta )
  if (length(k) != 1) {
    stop("Could not infer a unique k from length(thetas1). Please specify or check your inputs.")
  }
  
  if (method == "linear") {
    # Naive angle-wise interpolation
    theta_blend <- (1 - w)*thetas1 + w*thetas2
    Rblend <- dkrot(k, theta_blend)
    return(Rblend)
  } else if (method == "geodesic") {
    # Shortest path on SO(k)
    # R(t) = R1 * exp(t * log(R1^T * R2)), for t in [0,1]
    
    # 1) Construct the two rotation matrices
    R1 <- dkrot(k, thetas1)
    R2 <- dkrot(k, thetas2)
    
    # 2) Use matrix log/expm for the geodesic
    #    Rblend = R1 * exp( w * log( R1^T * R2 ) )
    library(expm)         # for expm(), logm()
    M <- logm(t(R1) %*% R2)
    Rblend <- R1 %*% expm(w * M)
    
    return(Rblend)
  }
}

dkrot_info <- function(k, thetas) {
  planes <- t(combn(k, 2))  # same ordering as in `dkrot()`
  planes <- data.frame(planes[order(apply(planes, 1, max)),])
  big_rotmats <- lapply(1:choose(k, 2), function(i){
    compute_rotmat_theta(planes[i,], thetas[i], k)
  })
  return(list(planes = planes,
              big_rotmats = big_rotmats))
}

dkrot_premult <- function(k, 
                          thetas_remaining, 
                          planes_remaining, 
                          precomputed_bg) {
  big_rotmats_remaining <- lapply(1:length(thetas_remaining), function(i){
    compute_rotmat_theta(planes_remaining[i,], thetas_remaining[i], k)
  })
  rotmat <- Reduce("%*%", c(list(precomputed_bg), big_rotmats_remaining))
  rotmat
}

mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}

#### try collapsing from the end of the forward rotations ####
# hcube_init <- tail(full_hcubes, 1)[[1]]
# R_init <- tail(rot_mats, 1)[[1]]
# k <- dim(R_init)[1]
# A_init <- rot_to_skew(R_init)
# n_steps <- 100
# shape_blend <- 1
# w <- blending_weights(n_steps, shape = shape_blend)
# A_random <- random_skew_symmetric(k) / 20
# A_mats <- lapply(w, function(wi){
#   blend_skew_symmetric(A_init, A_random, wi, method = "lie_group")
# })
# R_mats <- lapply(A_mats, skew_to_rot)
# cumulative_R <- Reduce("%*%", rev(R_mats))
# hcube_end <- t(cumulative_R %*% t(hcube_init))
# 
# #compute loss function
# angular_velocities <- sapply(R_mats, angular_velocity)
# av_loss <- sd(angular_velocities)
# alignment_loss <- sum((abs(hcube_end[,ncol(hcube_end)]) - 0.5)^2)
# total_loss <- av_loss + alignment_loss
# 
# #now try optimization:
# p_init <- vectorize_skew_symmetric(A_random)
# 
# # Now define the objective in terms of p
# 
# collapse_objective <- function(
    #     p,                  # Current guess for skew-symmetric parameters
#     R_init,             # Your "starting" rotation
#     hcube_init,         # The hypercube vertices from the end of forward rotations
#     n_steps,            # Number of frames (100, say)
#     angular_vel_fun,    # Your function: angular_velocity(R)
#     k
# ) {
#   # 1) Reconstruct the candidate skew-symmetric matrix
#   A_candidate <- reconstruct_skew_symmetric(p, k)
#   
#   # 2) Exponentiate once for "step matrix":
#   #    S = exp( A_candidate / n_steps )
#   S <- expm(A_candidate / n_steps)
#   
#   # 3) Generate each R_i by left-multiplying R_init with S^i
#   #    We'll measure the angular velocity for each R_i (or between frames).
#   #    Then, the final frame is R_n = R_init * S^n.
#   #    We'll do an iterative approach to avoid computing S^i by repeated multiplication.
#   R_frames <- vector("list", n_steps)
#   Ri_power <- diag(k)  # This holds S^i; starts at S^0 = I.
#   
#   for(i in seq_len(n_steps)) {
#     Ri_power <- Ri_power %*% S         # S^i
#     R_frames[[i]] <- R_init %*% Ri_power
#   }
#   
#   # 4) alignment_loss from the final frame
#   R_final <- R_frames[[n_steps]]
#   
#   # rotate hcube_init
#   hcube_end <- t(R_final %*% t(hcube_init))  # 2^k x k
#   # Example alignment: last dimension near ±0.5
#   alignment_loss <- sum((abs(hcube_end[, k]) - 0.5)^2)
#   
#   # 5) angular velocities
#   # If your "angular_velocity()" function is per-frame,
#   # you might consider consecutive pairs (R_frames[i], R_frames[i-1]) 
#   # or just measure each R_frames[i] vs identity. 
#   # For consistency with your earlier code, let's do one velocity per R_frames[i].
#   angular_velocities <- sapply(R_frames, angular_vel_fun)
#   av_loss <- sd(angular_velocities)
#   
#   # 6) total_loss
#   total_loss <- av_loss + alignment_loss
#   
#   # 7) Print for debugging
#   print(total_loss)
#   
#   return(total_loss)
# }
# 
# 
# # Initialize p (the skew-symmetric parameters for A_candidate)
# p_init <- vectorize_skew_symmetric( random_skew_symmetric(k) / 20 )
# res <- optimx(
#   par = p_init,
#   fn = collapse_objective,
#   R_init = R_init,
#   hcube_init = hcube_init,
#   n_steps = n_steps,
#   angular_vel_fun = angular_velocity,
#   k = k,
#   method = "L-BFGS-B",
#   control = list(
#     maxit = 200,
#     trace = 1
#   )
# )
# 
# # Extract best parameter vector
# p_opt <- as.numeric(res[1, paste0("p", seq_along(p_init))])
# A_opt <- reconstruct_skew_symmetric(p_opt, k)
# 
# # You can now build the final rotation from R_init:
# R_final <- R_init %*% expm(A_opt)
# 
# R_mats_opt <- lapply(A_mats_opt, skew_to_rot)
# cumulative_R <- Reduce("%*%", rev(R_mats_opt))
# hcube_end <- t(R_final %*% t(hcube_init))
# hcube_end

#...fail, let's try it with givens rotations instead?

#### givens approach ####


#now for the optimization

# Assume you have:
#   - undokrot(R) -> returns a vector of angles (length choose(k,2))
#   - dkrot(k, theta) -> returns a kxk rotation from Givens angles in each plane
#   - blend_givens(theta_init, theta_cand, w) -> do (1-w)*theta_init + w*theta_cand
#   - angular_velocity(R) as before
#   - blending_weights(n_steps, shape) for your w

invlogit <- function(x) exp(x) / (1+exp(x))

collapse_objective_givens <- function(
    p,                   # current guess for Givens angles (theta_cand)
    n_p, # number of p-vectors to blend together
    theta_init,          # Givens angles from R_init
    p_inds, # indices for the p angles to optimize over
    p_fixed, # vector of length choose(k,2) with all the fixed + 0 angles 
    hcube_init,          # hypercube from the forward rotation
    angular_vel_fun,     # your function "angular_velocity(R)"
    dkrot,           # function that builds a rotation from Givens angles
    k, #dimensionality of the hypercube
    precomputed_bg = precomputed_bg, #precomputed big matrix for givens rotations
    planes_remaining = planes_remaining, #planes we are still optimizing over
    print_loss = F, #should we print the loss at the current optimization step?
    loss_weights, #weights for the different components of our loss function
    n_steps #the number of discrete rotations over which the blending will occur
) {
  
  #construct the full matrix of 'p's
  length_full_p <- choose(k, 2)
  # all_ps <- matrix(theta_init, n_steps, length_full_p, byrow = T)
  n_angles_p <- prod(n_p-1, length(p_inds))
  sub_ps <- rbind(theta_init[p_inds], 
                  matrix(p[1:n_angles_p], n_p-1, length(p_inds), byrow = T))
  w <- blending_simplex_path(n_p, n_steps, 
                             shape =  0.8 + invlogit(p["shape_blend"]) * 0.4, pow_ends = 2)
  sub_ps_weighted <- w %*% sub_ps
  # all_ps[,p_inds] <- sub_ps_weighted
  
  #compose rotation matrices (make efficient later)
  R_mats <- vector("list", n_steps)
  for(i in seq_len(n_steps)) {
    # R_mats[[i]] <- dkrot(k, all_ps[i,])
    R_mats[[i]] <- dkrot_premult(k = k, 
                                 thetas_remaining = sub_ps_weighted[i,], 
                                 planes_remaining, 
                                 precomputed_bg)
    
  }
  
  # multiply them together
  cumulative_R <- Reduce("%*%", rev(R_mats))
  
  # 4) Rotate the hypercube
  hcube_end <- t(cumulative_R %*% t(hcube_init))
  
  # compute alignment_loss
  alignment_loss <- sum((abs(hcube_end[, k]) - 0.5)^2) * loss_weights[1]
  
  # regularize shape term of blending function
  shape_blend_loss <- abs(p["shape_blend"] + 1)
  
  # what if we let it align any dimension?
  # no, this is cool, but it just fixes it in a different axis :/
  # for it to properly collapse, it needs to be the last coordinate I think
  # alignment_loss <- min(apply(hcube_end, 2, function(x) sum((abs(x) - 0.5)^2)))
  
  # 6) angular_velocities
  # angular_vels <- sapply(R_mats, angular_vel_fun)
  # av_loss <- sd(angular_vels)
  av_loss <- 0
  
  # or just compare to existing theta_init[p_inds]?
  dev_loss <- sum((t(t(sub_ps[-1,]) - theta_init[p_inds]))^2) * loss_weights[2]
  
  total_loss <- alignment_loss + av_loss + dev_loss + shape_blend_loss
  
  #print to screen
  if(print_loss){
    print(total_loss)  
  }
  
  #return value
  return(total_loss)
}


#### more preprocessing ####
#which rotation angles are pinned? which can we optimize?
run_test <- F
if(run_test){
  
#set the stage
hcube_init <- tail(full_hcubes, 1)[[1]]
R_init <- tail(rot_mats, 1)[[1]]
k <- dim(R_init)[1]
givens_init <- undokrot(R_init)

#try out a random matrix
# n_steps <- 100
# shape_blend <- 1
# w <- blending_weights(n_steps, shape = shape_blend)
# givens_random <- runif(choose(k,2), 0, 2*pi)
# R_mats <- lapply(w, function(wi){
#   blend_givens(givens_init, givens_random, wi, method = "linear")
# })
# cumulative_R <- Reduce("%*%", rev(R_mats))
# hcube_end <- t(cumulative_R %*% t(hcube_init))
# hcube_end
# 
# #compute loss function
# angular_velocities <- sapply(R_mats, angular_velocity)
# av_loss <- sd(angular_velocities)
# alignment_loss <- sum((abs(hcube_end[,ncol(hcube_end)]) - 0.5)^2)
# total_loss <- av_loss + alignment_loss


n_angles <- choose(k, 2)
rm_fracs_givens <- lapply(rm_fracs, undokrot)
rmfg_0_angles <- lapply(rm_fracs_givens, function(x) which(abs(x) < 1E-6))[-1]
rmfg_n0_angles <- lapply(rm_fracs_givens, function(x) which(abs(x) > 1E-6))
rmfg_free_angles <- lapply(2:k, function(i) setdiff(rmfg_n0_angles[[i]], rmfg_n0_angles[[i-1]]))
rmfg_fixed_angles <- lapply(2:k, function(i) setdiff(1:n_angles, 
                                                     c(rmfg_free_angles[[i-1]], 
                                                       rmfg_0_angles[[i-1]])))

# specify inputs
R_init <- tail(rot_mats, 1)[[1]]
hcube_init <- tail(full_hcubes, 1)[[1]]
k <- nrow(R_init)
givens_matrices <- dkrot_info(k, tail(rm_fracs_givens, 1)[[1]])
precomputed_bg <- Reduce("%*%", givens_matrices$big_rotmats[rmfg_fixed_angles[[k-1]]])
planes_remaining <- givens_matrices$planes[rmfg_free_angles[[k-1]],]

# dkrot(k, tail(rm_fracs_givens, 1)[[1]]) -
# dkrot_premult(k, 
#               tail(rm_fracs_givens, 1)[[1]][p_inds], 
#               planes_remaining, 
#               precomputed_bg)

# initialize angles to last rotation
theta_init <- undokrot(R_init)      # length choose(k,2)

# specify parameters to optimize, as well as their starting value
n_p <- 2
p_inds <- tail(rmfg_free_angles, 1)[[1]]
p_init <- c(rep(theta_init[p_inds], n_p), 
            shape_blend = 0)

# specify iterative rotation metadata
n_steps <- 80
# shape_blend <- 1.2
# w <- blending_simplex_path(3, n_steps, shape = shape_blend, pow_ends = 2)

# optimize (first pass)
res <- optim(
  par = p_init,
  fn = collapse_objective_givens,
  theta_init = theta_init,
  n_p = n_p,
  p_inds = p_inds,
  hcube_init = hcube_init,
  angular_vel_fun = angular_velocity,
  dkrot = dkrot,     # your function
  k = k,
  precomputed_bg = precomputed_bg,
  planes_remaining = planes_remaining,
  n_steps = n_steps,
  loss_weights = c(coords = 1, theta = 20),
  method = "L-BFGS-B",
  control = list(
    maxit = 50,
    trace = 1
  )
)

# best_par <- as.numeric(res[1, paste0("p", seq_along(p_init))])
best_par <- res$par

# 6) Optimize with L-BFGS-B, focusing on hitting the target
res <- optim(
  par = best_par,
  fn = collapse_objective_givens,
  theta_init = theta_init,
  n_p = n_p,
  p_inds = p_inds,
  hcube_init = hcube_init,
  angular_vel_fun = angular_velocity,
  dkrot = dkrot,     # your function
  k = k,
  precomputed_bg = precomputed_bg,
  planes_remaining = planes_remaining,
  n_steps = n_steps,
  loss_weights = c(coords = 10, theta = 5),
  method = "L-BFGS-B",
  control = list(
    maxit = 50,
    trace = 1
  )
)

best_par <- res$par

# Build the final path:
#construct the full matrix of 'p's
length_full_p <- choose(k, 2)
n_angles_p <- prod(n_p-1, length(p_inds))
all_ps <- matrix(theta_init, n_steps, length_full_p, byrow = T)
sub_ps <- rbind(theta_init[p_inds], 
                matrix(best_par[1:n_angles_p], n_p-1, length(p_inds), byrow = T))
w <- blending_simplex_path(n_p, n_steps, 
                           shape =  0.8 + invlogit(best_par["shape_blend"]) * 0.4, pow_ends = 2)
sub_ps_weighted <- w %*% sub_ps
all_ps[,p_inds] <- sub_ps_weighted

#compose rotation matrices (make efficient later)
R_mats <- vector("list", n_steps)
for(i in seq_len(n_steps)) {
  R_mats[[i]] <- dkrot(k, all_ps[i,])
}

# multiply them together
cumulative_R <- Reduce("%*%", rev(R_mats))

# 4) Rotate the hypercube
hcube_end <- t(cumulative_R %*% t(hcube_init))


hcube_end[,ncol(hcube_end)]
final_rotation <- align_hypercube_axis(hcube_end)
t(final_rotation$rotation_matrix %*% t(hcube_end))[,ncol(hcube_end)]
det(final_rotation$rotation_matrix)

#woot woot! works well.

#compose rotation matrices (make efficient later)
R_mats[[n_steps]] <- final_rotation$rotation_matrix %*% R_mats[[n_steps]]

# multiply them together
cumulative_R <- Reduce("%*%", rev(R_mats))
hcube_end <- t(cumulative_R %*% t(hcube_init))
hcube_end

}


#### functions ####
#perfection! Now let's put it into a function.
# apply final correction
align_hypercube_axis <- function(hcube_coords, tol = 1e-9) {
  # --- Input Validation ---
  if (!is.matrix(hcube_coords) || !is.numeric(hcube_coords)) {
    stop("Input 'hcube_coords' must be a numeric matrix.")
  }
  n <- nrow(hcube_coords)
  k <- ncol(hcube_coords)
  if (k < 1) {
    stop("Input matrix must have at least one column (k >= 1).")
  }
  # Check if it looks like a hypercube (optional)
  # is_power_of_2 <- n > 0 && k >= 0 && ceiling(log2(n)) == floor(log2(n))
  # if (k > 0 && n != 2^k && is_power_of_2 && log2(n) != k) {
  #   warning(paste("Number of rows (", n, ") suggests dimension", log2(n), "but number of columns is k =", k))
  # } else if (k > 0 && !is_power_of_2) {
  #   warning(paste("Number of rows (", n, ") is not 2^k for k =", k, ". Input might not be a hypercube."))
  # }
  
  if (k == 1) {
    warning("k=1, rotation is trivial (identity matrix).")
    # Calculate side length even for k=1
    min_val <- min(hcube_coords[,1])
    max_val <- max(hcube_coords[,1])
    side_length = max_val - min_val
    return(list(rotation_matrix = matrix(1, 1, 1),
                rotated_coords = hcube_coords,
                side_length = side_length))
  }
  
  # --- Identify Hypercube Axis corresponding to last coordinate separation ---
  # Use the median of the last coordinate to split the points
  last_col_median <- median(hcube_coords[, k])
  group1_indices <- which(hcube_coords[, k] < last_col_median)
  group2_indices <- which(hcube_coords[, k] >= last_col_median)
  
  # Check if split is reasonable
  if (length(group1_indices) == 0 || length(group2_indices) == 0) {
    stop("Could not separate points based on the last coordinate. Is the input degenerate or already aligned perfectly?")
  }
  # For a hypercube, expect n/2 points in each group
  if (n == 2^k && (length(group1_indices) != n/2 || length(group2_indices) != n/2)) {
    warning(paste("Split based on median of last column resulted in uneven groups (",
                  length(group1_indices), "vs", length(group2_indices),
                  "). This might indicate precision issues or non-standard hypercube input.",
                  "Proceeding with calculated centroids."))
  }
  
  center1 <- colMeans(hcube_coords[group1_indices, , drop = FALSE])
  center2 <- colMeans(hcube_coords[group2_indices, , drop = FALSE])
  
  # This vector represents the direction (and length) of the hypercube's axis
  # that currently separates the points along the k-th dimension.
  u_direction <- center2 - center1
  
  # Calculate the side length along this axis
  side_length <- sqrt(sum(u_direction^2))
  if (side_length < tol) {
    stop("Hypercube seems collapsed along the k-th dimension axis (side length is near zero).")
  }
  
  # Normalize the direction vector
  u_norm <- u_direction / side_length
  
  # --- Define Target Vector ---
  # We want to align u_norm with the k-th standard basis vector e_k
  e_k <- c(rep(0, k - 1), 1)
  
  # --- Calculate Rotation Matrix R ---
  # R maps u_norm to e_k. We use the method based on reflections.
  # R = H_ek %*% H_prime, where H_prime maps u_norm -> -e_k, and H_ek maps -e_k -> e_k
  
  # Handle edge case: u_norm is already aligned with e_k
  if (sqrt(sum((u_norm - e_k)^2)) < tol) {
    R <- diag(k) # Identity matrix
    # Handle edge case: u_norm is anti-aligned with e_k
  } else if (sqrt(sum((u_norm + e_k)^2)) < tol) {
    # Need 180-degree rotation mapping -e_k to e_k. det(R) must be 1.
    # Rotate 180 degrees in the (x1, xk) plane.
    R <- diag(k)
    R[1, 1] <- -1
    R[k, k] <- -1
    # General case: use the reflection method
  } else {
    # H_prime = I - 2*outer(v,v)/dot(v,v), where v = u_norm + e_k
    v_prime <- u_norm + e_k
    v_prime_outer <- outer(v_prime, v_prime)
    v_prime_dot <- sum(v_prime^2) # Guaranteed > 0 because u_norm != -e_k
    
    H_prime <- diag(k) - 2 * v_prime_outer / v_prime_dot
    
    # H_ek reflects across the hyperplane orthogonal to e_k (flips k-th coord)
    H_ek <- diag(k)
    H_ek[k,k] <- -1
    # H_ek = diag(c(rep(1, k - 1), -1)) # Equivalent
    
    # The rotation R = H_ek %*% H_prime maps u_norm to e_k
    R <- H_ek %*% H_prime
  }
  
  rotated_coords <- hcube_coords %*% t(R)
  
  # --- Return Results ---
  return(list(rotation_matrix = R,
              rotated_coords = rotated_coords,
              side_length = side_length))
}

collapse_hcube <- function(hcube_init, R_init, 
                           n_steps = 100, n_p = 2, nrep = 10, 
                           n_steps_range = 0.1, trace = 0){
  
  #get useful variables 
  k <- dim(R_init)[1]
  n_angles <- choose(k, 2)
  rmfg_fixed_angles <- 1:choose(k-1, 2)
  rmfg_free_angles <- (choose(k-1, 2) + 1):n_angles
  
  # specify inputs for optimization
  theta_init <- undokrot(R_init)
  givens_matrices <- dkrot_info(k, theta_init)
  precomputed_bg <- Reduce("%*%", givens_matrices$big_rotmats[rmfg_fixed_angles])
  planes_remaining <- givens_matrices$planes[rmfg_free_angles,]
  
  # specify parameters to optimize, as well as their starting value
  p_inds <- rmfg_free_angles
  p_init <- c(rep(theta_init[p_inds], n_p), 
              shape_blend = -1)
  n_steps_vec <- round(n_steps * seq(1-n_steps_range, 1+n_steps_range, length.out = nrep))
  
  # compute blending weights (can also put inside optimization)
  # w <- blending_simplex_path(n_p, n_steps, shape = shape_blend, pow_ends = 2)
  
  # optimize (first pass), focusing on getting minimally deviation to rotation
  print("Performing Optimization:")
  res <- parallel::mclapply(1:nrep, function(nri){
    
    mcprint(paste0(nri, ".1 "))
    
    if(nri > 1){
      p_init <- (p_init + runif(length(p_init), -1, 1)) %% (2*pi)
    }
    
    res_init <- optim(
      par = p_init,
      fn = collapse_objective_givens,
      theta_init = theta_init,
      n_p = n_p,
      p_inds = p_inds,
      hcube_init = hcube_init,
      angular_vel_fun = angular_velocity,
      dkrot = dkrot,
      k = k,
      precomputed_bg = precomputed_bg,
      planes_remaining = planes_remaining,
      n_steps = n_steps_vec[nri],
      loss_weights = c(1, 5),
      method = "L-BFGS-B",
      print_loss = F,
      control = list(
        maxit = 40,
        trace = trace
      )
    )
    
    mcprint(paste0(nri, ".2 "))
    
    # optimize again, focusing on hitting the target
    res_final <- optim(
      par = res_init$par,
      fn = collapse_objective_givens,
      theta_init = theta_init,
      n_p = n_p,
      p_inds = p_inds,
      hcube_init = hcube_init,
      angular_vel_fun = angular_velocity,
      dkrot = dkrot,
      k = k,
      precomputed_bg = precomputed_bg,
      planes_remaining = planes_remaining,
      n_steps = n_steps_vec[nri],
      loss_weights = c(5, 1),
      method = "L-BFGS-B",
      print_loss = F,
      control = list(
        maxit = 20,
        trace = trace
      )
    )
    
    #modify final value to assess goodness of rotation
    value <- collapse_objective_givens(p = res_final$par,
                                       theta_init = theta_init,
                                       n_p = n_p,
                                       p_inds = p_inds,
                                       hcube_init = hcube_init,
                                       angular_vel_fun = angular_velocity,
                                       dkrot = dkrot,
                                       k = k,
                                       precomputed_bg = precomputed_bg,
                                       planes_remaining = planes_remaining,
                                       n_steps = n_steps,
                                       loss_weights = c(1, 5))
    
    return(list(res_final = res_final,
                value = value))
    
  }, mc.cores = 10)
  
  #extract parameter vector
  best_val <- which.min(sapply(res, function(x) x$value))
  n_steps <- n_steps_vec[best_val]
  best_res <- res[[best_val]]$res_final
  best_par <- best_res$par
  
  #build the final result's worth of matrices
  length_full_p <- choose(k, 2)
  n_angles_p <- prod(n_p-1, length(p_inds))
  all_ps <- matrix(theta_init, n_steps, length_full_p, byrow = T)
  sub_ps <- rbind(theta_init[p_inds], 
                  matrix(best_par[1:n_angles_p], n_p-1, length(p_inds), byrow = T))
  w <- blending_simplex_path(n_p, n_steps, 
                             shape =  0.8 + invlogit(best_par["shape_blend"]) * 0.4, pow_ends = 2)
  sub_ps_weighted <- w %*% sub_ps
  
  #compose rotation matrices (make efficient later)
  R_mats <- vector("list", n_steps)
  for(i in seq_len(n_steps)) {
    R_mats[[i]] <- dkrot_premult(k = k, 
                                 thetas_remaining = sub_ps_weighted[i,], 
                                 planes_remaining, 
                                 precomputed_bg)
    
  }
  
  
  # multiply them together
  cumulative_R <- Reduce("%*%", rev(R_mats))
  
  # compute provisional rotated hcube
  hcube_end <- t(cumulative_R %*% t(hcube_init))
  
  #make a final adjustment for exactness
  final_rotation <- align_hypercube_axis(hcube_end)
  
  #compose rotation matrices
  R_mats[[n_steps]] <- final_rotation$rotation_matrix %*% R_mats[[n_steps]]
  
  # perform the final multiplication
  cumulative_R <- Reduce("%*%", rev(R_mats))
  hcube_end <- t(cumulative_R %*% t(hcube_init))
  
  #get some "performance" indicators
  mean_abs_deviation_from_R_init <- sapply(R_mats, function(x) mean(abs(x - R_init)))
  mean_abs_deviation_from_givens <- apply(sub_ps_weighted, 1, function(x) mean(abs(x - theta_init[p_inds])))
  
  #compile object to return and exit
  out <- list(
    hcube_end = hcube_end,
    R_mats = R_mats,
    cumulative_R = cumulative_R,
    final_rotation = final_rotation,
    mean_abs_deviation_from_R_init = mean_abs_deviation_from_R_init,
    mean_abs_deviation_from_givens = mean_abs_deviation_from_givens,
    n_steps = n_steps,
    res = res,
    best_res = best_res
  )
  
  return(out)
}

test_collapse_function <- F
if(test_collapse_function){
  hcube_init <- tail(full_hcubes, 1)[[1]]
  R_init <- tail(rot_mats, 1)[[1]]
  collapsed_hcube <- collapse_hcube(hcube_init, R_init,
                                    n_steps = 100,
                                    n_p = 2)
}


#other functions
find_min_directed_rotation <- function(curr_line, starting_line, curr_rot_angle, tolerance = 1e-9) {
  
  # --- Pre-calculations ---
  # Check for zero rotation case: only possible if already aligned
  if (abs(curr_rot_angle) < tolerance) {
    delta_coords <- starting_line - curr_line # Check direct match P1->S1, P2->S2
    delta_coords_swap <- starting_line[2:1, ] - curr_line # Check swapped match P1->S2, P2->S1
    
    aligned <- (all(abs(delta_coords) < tolerance) || all(abs(delta_coords_swap) < tolerance))
    return(if (aligned) 0 else NA_real_)
  }
  
  # Calculate angles and radii using vectorized operations
  a_curr <- atan2(curr_line[, 2], curr_line[, 1])
  a_start <- atan2(starting_line[, 2], starting_line[, 1])
  r_curr <- sqrt(rowSums(curr_line^2))
  r_start <- sqrt(rowSums(starting_line^2))
  
  possible_angles <- c()
  direction_sign <- sign(curr_rot_angle) # +1 for CCW, -1 for CW
  
  # --- Loop through pairings (curr index i to start index j) ---
  for (i in 1:2) {
    for (j in 1:2) {
      # Check radius match first
      if (abs(r_curr[i] - r_start[j]) < tolerance) {
        
        # Calculate raw angle difference (start - curr)
        delta_angle_raw <- a_start[j] - a_curr[i]
        
        # Normalize to [0, 2*pi) for CCW angle
        delta_ccw <- (delta_angle_raw %% (2 * pi) + 2 * pi) %% (2 * pi)
        
        # Calculate equivalent CW angle (-2*pi, 0]
        # Handle case where delta_ccw is ~0 (raw angle was 0 or multiple of 2*pi)
        delta_cw <- if (abs(delta_ccw) < tolerance) 0 else delta_ccw - 2 * pi
        
        # Add the angle corresponding to the desired rotation direction
        angle_to_consider <- if (direction_sign > 0) delta_ccw else delta_cw
        
        # Add to list if it's non-zero OR if raw difference was zero
        # Also filter out angles in the "wrong" direction (e.g., negative if positive required)
        # Allow angles very close to zero.
        if ( (abs(angle_to_consider) > tolerance || abs(delta_angle_raw) < tolerance) && 
             (sign(angle_to_consider) == direction_sign || abs(angle_to_consider) < tolerance) ) {
          possible_angles <- c(possible_angles, angle_to_consider)
        }
      }
    }
  }
  
  # --- Select Minimum Magnitude Angle in Correct Direction ---
  if (length(possible_angles) == 0) {
    warning("No possible alignment found for the specified direction and tolerance.")
    return(NA_real_)
  }
  
  # Remove duplicates caused by floating point noise near 0/2pi/-2pi etc.
  possible_angles <- unique(round(possible_angles / tolerance)) * tolerance
  
  # Filter again just to be safe (shouldn't be needed with check inside loop)
  possible_angles <- possible_angles[sign(possible_angles) == direction_sign | abs(possible_angles) < tolerance]
  
  if (length(possible_angles) == 0) {
    warning("No possible alignment found after final filtering.")
    return(NA_real_)
  }
  
  # Find the angle with the minimum absolute value *in the desired direction*
  result_angle <- if (direction_sign > 0) {
    min(possible_angles) # Smallest positive angle
  } else {
    max(possible_angles) # Largest negative angle (closest to 0 from below)
  }
  
  # Ensure exactly zero if result is within tolerance
  return(if (abs(result_angle) < tolerance) 0 else result_angle)
}

generate_smooth_rotation_steps <- function(curr_rot_angle,
                                                       starting_rot_angle,
                                                       rotation_needed,
                                                       tolerance = 1e-9) {
  
  # --- Input Validation and Edge Cases (mostly same as fixed_ends) ---
  # Check for Zero Rotation Needed
  if (abs(rotation_needed) < tolerance) {
    if (abs(curr_rot_angle) < tolerance && abs(starting_rot_angle) < tolerance) {
      return(numeric(0))
    } else if (abs(curr_rot_angle + starting_rot_angle) < tolerance) {
      # Achievable with N=2 if start = -end
      return(c(curr_rot_angle, starting_rot_angle))
    } else {
      warning("Target rotation is zero, but start/end angles don't allow for simple cancellation.")
      return(NA_real_)
    }
  }
  
  avg_angle_initial = (curr_rot_angle + starting_rot_angle) / 2
  
  # Check for Direction Conflict
  if (abs(avg_angle_initial) > tolerance && sign(rotation_needed) != sign(avg_angle_initial)) {
    warning("Required rotation direction conflicts with the average step direction.")
    return(NA_real_)
  }
  
  # --- Estimate Number of Steps (N) ---
  n_steps <- NA
  
  # Case 1: Start and End angles are the same
  if (abs(curr_rot_angle - starting_rot_angle) < tolerance) {
    if (abs(curr_rot_angle) < tolerance) {
      warning("Cannot reach non-zero rotation if start/end angles are zero.")
      return(NA_real_)
    }
    n_steps_est <- rotation_needed / curr_rot_angle
    if (n_steps_est < 1 - tolerance) {
      warning("Rotation needed cannot be achieved with constant positive angle steps.")
      return(NA_real_)
    }
    n_steps <- max(1, round(n_steps_est))
    if (n_steps == 1 && abs(rotation_needed - curr_rot_angle) > tolerance) {
      warning("Single step required but target rotation doesn't match start/end angle.")
      return(NA_real_)
    }
    # If start=end, the only smooth profile is constant.
    # Return constant steps if possible.
    return(rep(rotation_needed / n_steps, n_steps))
    
  } else {
    # Case 2: Start and End angles differ
    if (abs(avg_angle_initial) < tolerance) {
      if (abs(curr_rot_angle + starting_rot_angle - rotation_needed) < tolerance) {
        n_steps <- 2 # N=2 is exact and the only possibility
      } else {
        warning("Average step angle is zero, cannot achieve this total rotation.")
        return(NA_real_)
      }
    } else {
      n_steps_est <- rotation_needed / avg_angle_initial
      n_steps <- max(2, round(n_steps_est))
      if(n_steps < 2) n_steps <- 2 # Force N>=2 for transition
    }
    # Check if N=2 works, otherwise force N=3 minimum for adjustment
    if (n_steps == 2 && abs(curr_rot_angle + starting_rot_angle - rotation_needed) > tolerance) {
      n_steps <- 3
    }
  }
  
  # --- Generate Final Step Angles ---
  
  # Handle N=1 (Should only be reached if start=end=needed)
  if (n_steps == 1) {
    # This case should have been fully handled by the start=end logic above
    return(curr_rot_angle)
  }
  
  # Handle N=2
  if (n_steps == 2) {
    # N=2 is only valid if start + end == needed (checked above)
    return(c(curr_rot_angle, starting_rot_angle))
  }
  
  # Handle N >= 3: Use smooth base + smooth adjustment
  if (n_steps >= 3) {
    # 1. Generate ideal smooth steps using cosine blend
    t <- (0:(n_steps - 1)) / (n_steps - 1)
    blend_factor <- (1 - cos(pi * t)) / 2
    ideal_step_angles <- curr_rot_angle + (starting_rot_angle - curr_rot_angle) * blend_factor
    
    # Verify ends (should be correct by construction)
    # ideal_step_angles[1] should be curr_rot_angle
    # ideal_step_angles[n_steps] should be starting_rot_angle
    
    # 2. Calculate sum error
    sum_ideal <- sum(ideal_step_angles)
    delta_sum <- rotation_needed - sum_ideal
    
    # If error is negligible, we are done
    if (abs(delta_sum) < tolerance) {
      return(ideal_step_angles)
    }
    
    # 3. Generate smooth adjustment profile (Sine based, zero at ends)
    # Profile shape: sin(pi * t)
    adjustment_profile <- sin(pi * t)
    # Ensure ends are exactly zero due to potential float issues
    adjustment_profile[1] <- 0
    adjustment_profile[n_steps] <- 0
    
    sum_profile <- sum(adjustment_profile)
    
    # 4. Calculate scaling factor for adjustment
    # Avoid division by zero if profile sum is somehow zero (e.g., N=2?)
    if (abs(sum_profile) < tolerance) {
      # This shouldn't happen for N>=3 unless delta_sum is also zero (handled above)
      # If delta_sum is non-zero but profile sum is zero, we cannot adjust smoothly.
      warning("Cannot apply smooth adjustment: Adjustment profile sum is zero.")
      # Fallback: return the unadjusted ideal steps? Or NA?
      # Returning ideal steps fails the sum constraint. NA is safer.
      return(NA_real_)
    }
    k <- delta_sum / sum_profile
    
    # 5. Apply adjustment
    final_step_angles <- ideal_step_angles + k * adjustment_profile
    
    # Ensure start/end are *exactly* the target values after float arithmetic
    final_step_angles[1] <- curr_rot_angle
    final_step_angles[n_steps] <- starting_rot_angle
    
    # Optional Verification: Check final sum
    # final_sum <- sum(final_step_angles)
    # print(paste("Final Sum Verification:", final_sum, "vs", rotation_needed))
    
    return(final_step_angles)
  }
  
  # Fallback
  warning("Failed to determine valid step count or configuration.")
  return(NA_real_)
}
