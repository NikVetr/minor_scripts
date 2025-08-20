# Keep necessary libraries and helper functions
library(expm)
library(optimx)

# --- Helper Functions (vectorize/reconstruct_skew_symmetric, random*, angular_velocity, blending_weights, rotation_distance) ---
# (Include the robust versions from the previous response)

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

skew2rot <- function(A) {
  R <- expm(A)
  svd_res <- svd(R)
  R <- svd_res$u %*% diag(c(rep(1, k-1), det(svd_res$u %*% t(svd_res$v)))) %*% t(svd_res$v)
  return(R)
}

# Generate a random rotation matrix in SO(k)
random_rotation_matrix <- function(k) {
  skew2rot(random_skew_symmetric(k))
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
blend_rotations_lie <- function(R1, R2, w) {
  k <- nrow(R1)
  # Ensure inputs are SO(k) - optional but safer
  # R1 = project_so_k(R1)
  # R2 = project_so_k(R2)
  
  R_rel <- try(t(R1) %*% R2, silent=TRUE)
  if(inherits(R_rel, "try-error")) {
    warning("Matrix multiplication failed in blend_rotations_lie (t(R1) %*% R2)")
    return(diag(k)) # Return identity on failure? Or R1? Needs thought. Let's return R1.
    return(R1)
  }
  
  # Project R_rel onto SO(k)
  svd_res <- svd(R_rel)
  det_uv <- det(svd_res$u %*% t(svd_res$v))
  diag_mat <- diag(c(rep(1, k - 1), det_uv))
  R_rel_proj <- svd_res$u %*% diag_mat %*% t(svd_res$v)
  
  if (norm(R_rel_proj - diag(k), type = "F") < sqrt(.Machine$double.eps) * k) {
    return(R1) # No relative rotation, return R1
  }
  
  A_rel <- try(logm(R_rel_proj), silent = TRUE)
  if(inherits(A_rel, "try-error")){
    warning("logm failed in blend_rotations_lie")
    return(R1) # Return R1 if log fails
  }
  A_rel_skew <- 0.5 * (A_rel - t(A_rel))
  A_scaled <- w * A_rel_skew
  
  R_exp_scaled <- try(expm(A_scaled), silent=TRUE)
  if(inherits(R_exp_scaled, "try-error")){
    warning("expm failed in blend_rotations_lie")
    return(R1) # Return R1 if exp fails
  }
  
  R_blend <- try(R1 %*% R_exp_scaled, silent = TRUE)
  if(inherits(R_blend, "try-error")){
    warning("Matrix multiplication failed in blend_rotations_lie (R1 %*% expm(...))")
    return(R1) # Return R1 if final multiply fails
  }
  
  # Optional: Project final result onto SO(k)
  svd_blend <- svd(R_blend)
  det_uv_blend <- det(svd_blend$u %*% t(svd_blend$v))
  diag_mat_blend <- diag(c(rep(1, k - 1), det_uv_blend))
  R_blend_proj <- svd_blend$u %*% diag_mat_blend %*% t(svd_blend$v)
  
  return(R_blend_proj)
}

# Blending weights function (logistic or linear)
blending_weights <- function(n, shape = 0) {
  stopifnot(n >= 1) # Can have n=1
  if (n == 1) return(0) # Or maybe 0.5? Let's assume 0 for M1=R_start. Check objective logic later.
  # If w[1]=0, M1=blend(R_start, M_target, 0) = R_start. Yes, w[1]=0 is correct.
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

# --- New Optimization Function ---

optimize_step_sequence_product <- function(R_start_step, R_target_product, n, lambda,
                                           initial_A_step_target = NULL, initial_shape = 0,
                                           optim_method = "L-BFGS-B",
                                           optim_control = list(trace = 1),
                                           maxit = 100) {
  
  k <- nrow(R_start_step)
  num_A_params <- k * (k - 1) / 2
  
  # Ensure R_start_step and R_target_product are valid rotation matrices
  # (Add projection here if needed, or assume valid inputs)
  
  # Define the objective function
  objective_fn <- function(params, R_start_step, R_target_product, n, lambda, k) {
    
    # --- 1. Decode Parameters ---
    if(length(params) != num_A_params + 1) {
      warning("Incorrect number of parameters in objective function.")
      return(1e10)
    }
    A_step_target_params <- params[1:num_A_params]
    shape <- params[num_A_params + 1]
    
    A_step_target <- try(reconstruct_skew_symmetric(A_step_target_params, k), silent = TRUE)
    if (inherits(A_step_target, "try-error")) {
      warning("Failed to reconstruct A_step_target.")
      return(1e10)
    }
    
    # --- 2. Generate Sequence of *Step* Matrices ---
    M_step_target <- try(expm(A_step_target), silent = TRUE)
    if (inherits(M_step_target, "try-error")) {
      warning("expm(A_step_target) failed.")
      return(1e10)
    }
    # Project M_step_target onto SO(k) just in case
    svd_mtarget <- svd(M_step_target)
    det_uv_mtarget <- det(svd_mtarget$u %*% t(svd_mtarget$v))
    diag_mat_mtarget <- diag(c(rep(1, k - 1), det_uv_mtarget))
    M_step_target_proj <- svd_mtarget$u %*% diag_mat_mtarget %*% t(svd_mtarget$v)
    
    w_vec <- blending_weights(n, shape)
    
    M_sequence <- vector("list", n)
    valid_steps = TRUE
    for (i in 1:n) {
      M_i <- blend_rotations_lie(R_start_step, M_step_target_proj, w_vec[i])
      # Check if blend returned a valid matrix (though blend should handle internal errors)
      if (!is.matrix(M_i) || nrow(M_i) != k || ncol(M_i) != k || any(!is.finite(M_i))) {
        warning(paste("Invalid matrix generated for step", i))
        valid_steps = FALSE
        break # Stop generating if one step fails badly
      }
      M_sequence[[i]] <- M_i
    }
    
    if(!valid_steps) { return(1e10) } # Return large loss if sequence generation failed
    
    # --- 3. Calculate Loss Components ---
    
    # a) Distance of Product to Target
    # Check if n=0 edge case (shouldn't happen if n>=1 enforced)
    if (n == 0) {
      R_product <- diag(k) # Product of no matrices is identity
    } else {
      R_product <- try(Reduce("%*%", M_sequence), silent = TRUE)
      if (inherits(R_product, "try-error")) {
        warning("Reduce('%*%', M_sequence) failed.")
        return(1e10) # Product failed
      }
      # Project final product onto SO(k) before distance calc
      svd_prod <- svd(R_product)
      det_uv_prod <- det(svd_prod$u %*% t(svd_prod$v))
      diag_mat_prod <- diag(c(rep(1, k - 1), det_uv_prod))
      R_product_proj <- svd_prod$u %*% diag_mat_prod %*% t(svd_prod$v)
      
      R_product <- R_product_proj # Use projected version
    }
    
    dist_target_product <- rotation_distance(R_product, R_target_product)
    
    # b) Variance in Step Angular Velocity
    if (n >= 2) { # Need at least 2 steps for variance
      step_angular_velocities <- sapply(M_sequence, angular_velocity, delta_t = 1)
      
      # Check for non-finite velocities returned by angular_velocity()
      valid_vels <- is.finite(step_angular_velocities)
      if(sum(valid_vels) < 2){
        var_ang_vel <- 0 # Not enough valid points for variance
      } else {
        var_ang_vel <- var(step_angular_velocities[valid_vels])
      }
      if(is.na(var_ang_vel)) var_ang_vel <- 0 # Handle case where all valid vels might be identical
    } else {
      var_ang_vel <- 0 # Not enough steps
    }
    
    # --- 4. Combine Loss ---
    loss <- dist_target_product + lambda * var_ang_vel
    
    if (!is.finite(loss)) {
      warning("Non-finite loss encountered.")
      loss <- 1e10
    }
    
    # Optional intermediate print:
    # cat("Shape:", signif(shape, 3), " A_target norm:", signif(norm(A_step_target,"F"),3),
    #     " -> Loss:", signif(loss, 5), " (Dist:", signif(dist_target_product, 5),
    #     ", VarVel:", signif(var_ang_vel, 5), ")\n")
    
    return(loss)
  }
  
  # --- Set up Initial Parameters ---
  if (is.null(initial_A_step_target)) {
    # Guess: If all steps were identical, say M, then M^n = R_target_product.
    # So M approx expm(logm(R_target_product)/n).
    # Let initial A_step_target be logm(R_target_product)/n.
    
    # Project R_target_product first
    svd_rtarget <- svd(R_target_product)
    det_uv_rtarget <- det(svd_rtarget$u %*% t(svd_rtarget$v))
    diag_mat_rtarget <- diag(c(rep(1, k - 1), det_uv_rtarget))
    R_target_proj <- svd_rtarget$u %*% diag_mat_rtarget %*% t(svd_rtarget$v)
    
    if (norm(R_target_proj - diag(k), type = "F") < sqrt(.Machine$double.eps) * k) {
      initial_A_step_target <- matrix(0, k, k) # Target is identity
    } else {
      A_target_total <- try(logm(R_target_proj), silent=TRUE)
      if(inherits(A_target_total, "try-error")) {
        warning("logm(R_target_product) failed for initial guess. Using zero matrix.")
        initial_A_step_target <- matrix(0, k, k)
      } else {
        initial_A_step_target <- (0.5 * (A_target_total - t(A_target_total))) / n # Ensure skew-symm and scale
      }
    }
  }
  
  initial_A_params <- vectorize_skew_symmetric(initial_A_step_target)
  initial_params <- c(initial_A_params, initial_shape)
  
  # Check initial parameters
  initial_loss <- objective_fn(initial_params, R_start_step, R_target_product, n, lambda, k)
  if (!is.finite(initial_loss)) {
    warning("Initial parameters produce non-finite loss (", initial_loss, "). Optimization might struggle. Trying zero A_step_target.")
    initial_A_step_target <- matrix(0,k,k)
    initial_A_params <- vectorize_skew_symmetric(initial_A_step_target)
    initial_params <- c(initial_A_params, initial_shape)
    initial_loss <- objective_fn(initial_params, R_start_step, R_target_product, n, lambda, k)
    if (!is.finite(initial_loss)) {
      stop("Even zero initial A_step_target produces non-finite loss. Check inputs (R_start_step?).")
    }
  }
  cat("Initial Loss:", initial_loss, "\n")
  
  # --- Run Optimization ---
  cat("Starting optimization...\n")
  opt_results <- optimx(par = initial_params,
                        fn = objective_fn,
                        method = optim_method,
                        control = c(optim_control, list(maxit=maxit, kkt=FALSE)),
                        R_start_step = R_start_step,
                        R_target_product = R_target_product,
                        n = n,
                        lambda = lambda,
                        k = k)
  cat("Optimization finished.\n")
  
  # --- Process Results ---
  print(summary(opt_results))
  valid_results <- opt_results[is.finite(opt_results$value), ]
  if (nrow(valid_results) == 0) {
    warning("Optimization failed to find any valid solution. Returning initial guess.")
    optimal_A_step_target <- initial_A_step_target
    optimal_shape <- initial_shape
    best_result_value <- initial_loss
  } else {
    best_result_idx <- which.min(valid_results$value)
    best_result <- valid_results[best_result_idx, ]
    best_result_value <- best_result$value
    
    optimal_params <- as.numeric(best_result[1:(num_A_params + 1)])
    optimal_A_params <- optimal_params[1:num_A_params]
    optimal_shape <- optimal_params[num_A_params + 1]
    optimal_A_step_target <- reconstruct_skew_symmetric(optimal_A_params, k)
  }
  
  # --- Post-Optimization Analysis ---
  optimal_M_step_target <- expm(optimal_A_step_target)
  # Project just in case
  svd_mtarget_opt <- svd(optimal_M_step_target)
  det_uv_mtarget_opt <- det(svd_mtarget_opt$u %*% t(svd_mtarget_opt$v))
  diag_mat_mtarget_opt <- diag(c(rep(1, k - 1), det_uv_mtarget_opt))
  optimal_M_step_target_proj <- svd_mtarget_opt$u %*% diag_mat_mtarget_opt %*% t(svd_mtarget_opt$v)
  
  optimal_w_vec <- blending_weights(n, optimal_shape)
  optimal_M_sequence <- vector("list", n)
  for (i in 1:n) {
    optimal_M_sequence[[i]] <- blend_rotations_lie(R_start_step, optimal_M_step_target_proj, optimal_w_vec[i])
  }
  
  final_R_product <- diag(k)
  if (n > 0) {
    final_R_product_calc <- try(Reduce("%*%", optimal_M_sequence), silent=TRUE)
    if(!inherits(final_R_product_calc, "try-error")) {
      svd_finalprod <- svd(final_R_product_calc)
      det_uv_finalprod <- det(svd_finalprod$u %*% t(svd_finalprod$v))
      diag_mat_finalprod <- diag(c(rep(1, k - 1), det_uv_finalprod))
      final_R_product <- svd_finalprod$u %*% diag_mat_finalprod %*% t(svd_finalprod$v)
    } else {
      warning("Failed to calculate final product for reporting.")
      # final_R_product remains identity
    }
  }
  
  final_dist_target_product <- rotation_distance(final_R_product, R_target_product)
  
  if (n >= 2) {
    final_step_angular_velocities <- sapply(optimal_M_sequence, angular_velocity, delta_t = 1)
    valid_final_vels <- is.finite(final_step_angular_velocities)
    if(sum(valid_final_vels) < 2){
      final_var_ang_vel <- 0
    } else {
      final_var_ang_vel <- var(final_step_angular_velocities[valid_final_vels])
    }
    if(is.na(final_var_ang_vel)) final_var_ang_vel <- 0
  } else {
    final_var_ang_vel <- 0
  }
  
  recalculated_final_loss <- final_dist_target_product + lambda * final_var_ang_vel
  
  # Return useful information
  return(list(
    optimal_A_step_target = optimal_A_step_target,
    optimal_M_step_target = optimal_M_step_target_proj, # The target step matrix
    optimal_shape = optimal_shape,
    optimal_weights = optimal_w_vec,
    optimal_step_sequence = optimal_M_sequence, # The sequence M_1, ..., M_n
    final_product = final_R_product,             # The product M_1 %*% ... %*% M_n
    final_loss_reported = best_result_value,
    final_loss_recalculated = recalculated_final_loss,
    final_distance_product_to_target = final_dist_target_product,
    final_variance_step_angular_velocity = final_var_ang_vel,
    initial_A_step_target = initial_A_step_target,
    initial_shape = initial_shape,
    optim_results_summary = summary(opt_results)
    # optim_results_raw = opt_results
  ))
}


#### test out ####
set.seed(456)
k <- 3 # Dimension
n <- 20 # Number of steps in the sequence

# Define the *target* product matrix
R_target_product <- random_rotation_matrix(k)

# Define the *initial step* matrix (intended to be small)
# Method 1: Fraction of the target product's rotation
# logm_target <- logm(R_target_product) # Assume R_target_product is SO(k)
# A_start_step <- (0.5 * (logm_target - t(logm_target))) / n
# R_start_step <- expm(A_start_step)

# Method 2: Fraction of some other base rotation (as per user)
R_base_for_step <- random_rotation_matrix(k)
logm_base <- logm(R_base_for_step) # Assume SO(k)
A_start_step <- (0.5 * (logm_base - t(logm_base))) / n
R_start_step <- expm(A_start_step)
# Project R_start_step just to be sure
svd_rstart <- svd(R_start_step)
R_start_step <- svd_rstart$u %*% diag(c(rep(1, k-1), det(svd_rstart$u %*% t(svd_rstart$v)))) %*% t(svd_rstart$v)


print("Norm of R_start_step generator (target is small rotation):")
print(norm(A_start_step, "F"))


# Set the trade-off parameter
lambda <- 0.05 # Adjust this: higher lambda -> smoother steps, potentially less accurate product

# Perform the optimization
optim_output_steps <- optimize_step_sequence_product(
  R_start_step = R_start_step,
  R_target_product = R_target_product,
  n = n,
  lambda = lambda,
  optim_method = "L-BFGS-B", # Or try others like "nlminb", "BFGS"
  optim_control = list(trace = 0, maxit=200), # Set trace > 0 for progress
  maxit = 200
)

# --- Verification ---
cat("\n--- Verifying the Optimized Step Sequence ---\n")

final_M_sequence <- optim_output_steps$optimal_step_sequence
final_R_product_actual <- diag(k)
if (length(final_M_sequence) > 0) {
  final_R_product_actual_calc <- try(Reduce("%*%", final_M_sequence), silent=TRUE)
  if(!inherits(final_R_product_actual_calc, "try-error")){
    svd_finalprod_act <- svd(final_R_product_actual_calc)
    det_uv_finalprod_act <- det(svd_finalprod_act$u %*% t(svd_finalprod_act$v))
    diag_mat_finalprod_act <- diag(c(rep(1, k - 1), det_uv_finalprod_act))
    final_R_product_actual <- svd_finalprod_act$u %*% diag_mat_finalprod_act %*% t(svd_finalprod_act$v)
  }
}


# Check product distance
dist_prod_target <- rotation_distance(final_R_product_actual, R_target_product)
cat(sprintf("Distance between actual product (M1 * ... * Mn) and R_target_product: %.4e\n", dist_prod_target))
cat(sprintf("(Should match reported distance: %.4e)\n", optim_output_steps$final_distance_product_to_target))
if (abs(dist_prod_target - optim_output_steps$final_distance_product_to_target) < 1e-7) {
  cat("Verification PASSED: Product distance matches reported value.\n")
} else {
  cat("Verification WARNING: Product distance differs from reported value.\n")
}
cat("\n")

# Check first step
if (n > 0) {
  first_step_M1 <- final_M_sequence[[1]]
  dist_m1_rstart <- rotation_distance(first_step_M1, R_start_step)
  cat(sprintf("Distance between first generated step (M1) and input R_start_step: %.4e\n", dist_m1_rstart))
  # Since w[1]=0, M1 = blend(R_start, M_target, 0) = R_start. This should be near zero.
  if (abs(dist_m1_rstart) < 1e-9) {
    cat("Verification PASSED: First step M1 is effectively R_start_step.\n")
  } else {
    cat("Verification WARNING: First step M1 differs significantly from R_start_step.\n")
  }
}
cat("\n")

# Check variance
if (n >= 2) {
  step_vels <- sapply(final_M_sequence, angular_velocity, delta_t=1)
  valid_step_vels <- step_vels[is.finite(step_vels)]
  actual_var = 0
  if(length(valid_step_vels) >=2){
    actual_var <- var(valid_step_vels)
  }
  cat(sprintf("Variance of calculated step angular velocities: %.4e\n", actual_var))
  cat(sprintf("(Should match reported variance: %.4e)\n", optim_output_steps$final_variance_step_angular_velocity))
  if (abs(actual_var - optim_output_steps$final_variance_step_angular_velocity) < 1e-7) {
    cat("Verification PASSED: Step variance matches reported value.\n")
  } else {
    cat("Verification WARNING: Step variance differs from reported value.\n")
  }
}

cat("\n--- Verification Complete ---\n")

# --- Plotting ---
# Plot blending weights
plot(1:n, optim_output_steps$optimal_weights, type = 'b', pch=19,
     xlab = "Step Index (i)", ylab = "Blending Weight (w_i)",
     main = paste("Optimized Blending Weights (Shape =", round(optim_output_steps$optimal_shape, 3), ")"))

# Plot step angular velocities
if (n >= 1) {
  step_vels <- sapply(optim_output_steps$optimal_step_sequence, angular_velocity, delta_t=1)
  valid_vel_indices <- which(is.finite(step_vels))
  if(length(valid_vel_indices) > 0) {
    plot(valid_vel_indices, step_vels[valid_vel_indices], type = 'b', pch = 19,
         xlab = "Step Index (i)", ylab = "Step Angular Velocity Magnitude (rad/step)",
         main = paste("Step Angular Velocities (Variance =", signif(optim_output_steps$final_variance_step_angular_velocity, 4), ")"),
         xlim = c(1, n))
  } else {
    print("Could not plot step angular velocities (no valid velocities).")
  }
  
  # Compare to angular velocity of R_start_step and M_step_target
  abline(h = angular_velocity(R_start_step), col = "blue", lty = 2)
  abline(h = angular_velocity(optim_output_steps$optimal_M_step_target), col = "red", lty = 2)
  legend("topright", legend = c("Step Velocities", "Vel(R_start_step)", "Vel(M_target_step)"),
         col = c("black", "blue", "red"), lty = c(1, 2, 2), pch = c(19, NA, NA), bg="white")
}


# --- (Previous code: All function definitions and the optimization call) ---

# Perform the optimization for the step sequence product
# optim_output_steps <- optimize_step_sequence_product(...)

# --- Verification Section (for optimize_step_sequence_product output) ---

cat("\n--- Verifying the Optimized Step Sequence ---\n")

# Recover the optimized sequence of STEP matrices (M_1, ..., M_n)
optimal_M_sequence <- optim_output_steps$optimal_step_sequence
n_seq <- length(optimal_M_sequence)

# Recover the inputs and key outputs for comparison
R_start_step_input <- optim_output_steps$initial_A_step_target # Wait, no, the input R_start_step wasn't returned. Need to use the one passed to the function.
# Let's assume R_start_step is available from the example scope.
R_target_product_input <- R_target_product # Assume available from example scope
optimal_M_step_target <- optim_output_steps$optimal_M_step_target # The step target found by optimization

if (n_seq == 0) {
  cat("Error: The optimized sequence of steps is empty.\n")
} else {
  
  # 1. Compare the FIRST STEP matrix (M_1) to the input R_start_step
  cat("--- Check 1: First Step vs R_start_step ---\n")
  first_step_M1 <- optimal_M_sequence[[1]]
  dist_m1_to_rstart <- rotation_distance(first_step_M1, R_start_step) # Use R_start_step from the calling scope
  cat(sprintf("Distance between the first generated step (M1) and the input R_start_step: %.4e\n", dist_m1_to_rstart))
  # Because the blending weight w_1 = 0, M1 = blend(R_start_step, M_target, 0) which should equal R_start_step.
  if (abs(dist_m1_to_rstart) < 1e-9) {
    cat("Verification PASSED: The first step M1 is effectively R_start_step (as expected since w_1=0).\n")
  } else {
    cat("Verification WARNING: The first step M1 differs significantly from R_start_step.\n")
  }
  cat("\n")
  
  # 2. Observe the Progression of Steps
  cat("--- Check 2: Progression of Step Matrices ---\n")
  cat("Distances of each step M_i from the initial step R_start_step:\n")
  for (i in 1:n_seq) {
    dist_mi_to_rstart <- rotation_distance(optimal_M_sequence[[i]], R_start_step)
    cat(sprintf("  Distance(M_%d, R_start_step): %.4e\n", i, dist_mi_to_rstart))
  }
  cat("Distances of each step M_i from the optimized target step M_step_target:\n")
  for (i in 1:n_seq) {
    dist_mi_to_mtarget <- rotation_distance(optimal_M_sequence[[i]], optimal_M_step_target)
    cat(sprintf("  Distance(M_%d, M_step_target): %.4e\n", i, dist_mi_to_mtarget))
  }
  cat("(Expect Distance(M_i, R_start_step) to generally increase as i increases)\n")
  cat("(Expect Distance(M_i, M_step_target) to generally decrease as i increases)\n")
  # Note: The exact pattern depends heavily on the optimal shape parameter.
  cat("\n")
  
  
  # 3. Compare the CUMULATIVE PRODUCT of steps to R_target_product
  cat("--- Check 3: Cumulative Product vs R_target_product ---\n")
  # Recalculate the product robustly
  final_product_recalculated <- diag(k)
  if (n_seq > 0) {
    product_calc <- try(Reduce("%*%", optimal_M_sequence), silent = TRUE)
    if(!inherits(product_calc, "try-error")) {
      # Project final product onto SO(k)
      svd_prod <- svd(product_calc)
      det_uv_prod <- det(svd_prod$u %*% t(svd_prod$v))
      diag_mat_prod <- diag(c(rep(1, k - 1), det_uv_prod))
      final_product_recalculated <- svd_prod$u %*% diag_mat_prod %*% t(svd_prod$v)
    } else {
      cat("Warning: Failed to recalculate product using Reduce for verification.\n")
    }
  }
  
  # Compare recalculated product to the target input
  dist_product_to_target <- rotation_distance(final_product_recalculated, R_target_product_input)
  # cat(sprintf("Distance between the cumulative product (M1 %*% ... %*% Mn) and R_target_product: %.4e\n", dist_product_to_target))
  cat(sprintf("(This should ideally be close to the reported 'final_distance_product_to_target': %.4e)\n", optim_output_steps$final_distance_product_to_target))
  
  if (abs(dist_product_to_target - optim_output_steps$final_distance_product_to_target) < 1e-7) {
    cat("Verification CONSISTENT: The calculated product distance matches the reported final distance.\n")
  } else {
    cat("Verification WARNING: The calculated product distance differs significantly from the reported final distance.\n")
  }
  # Assess how well the optimization met the primary goal (getting the product right)
  if (abs(dist_product_to_target) < 1e-6) { # Threshold for being "close"
    cat("Verification PASSED: The cumulative product is close to R_target_product.\n")
  } else {
    cat("Verification NOTE: The cumulative product is NOT very close to R_target_product.\n")
    cat("(This might be acceptable depending on lambda and optimization success.)\n")
  }
  
}

cat("\n--- Verification Complete ---\n")

optimal_M_sequence[[2]]
R_start_step
