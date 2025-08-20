# --- Load Libraries & Ensure Helper Functions are Defined ---
library(expm)
library(pracma) # If using cross product in find_rotation_aligning for k=3

# Need these functions defined (robust versions):
# find_rotation_aligning (using Householder is good)
# project_so_k
# blending_weights
# rotation_distance (for verification)
# logm (base or from expm)

# --- Inputs ---
# Ensure these variables exist from your previous setup:
# rot_mats: List of step matrices from unfolding phase
# cumrot_mats: List of cumulative matrices from unfolding
# k: Dimension
# num_collapse_steps_per_dim: How many animation frames/steps per dim collapse
# comp_shape: Shape parameter for blending_weights (e.g., 1 or 2 for ease-in/out)

cat("Setting up Composition collapse phase...\n")

# --- Initialization ---
all_collapse_steps_comp <- list()
if (length(cumrot_mats) == 0) { stop("cumrot_mats list is empty.") }
R_current_comp <- project_so_k(cumrot_mats[[length(cumrot_mats)]]) # Start orientation SO(k)

if (length(rot_mats) == 0) { R_trend_step <- diag(k) } else { R_trend_step <- project_so_k(rot_mats[[length(rot_mats)]]) } # Initial trend SO(k)

# --- Composition Collapse Loop ---
start_dim_collapse <- k
end_dim_collapse <- 3 # Or 2

if (start_dim_collapse < end_dim_collapse) {
  cat("No dimensions to collapse.\n")
} else {
  for (dim_to_collapse in start_dim_collapse:end_dim_collapse) {
    cat(sprintf("\n--- Starting Composition Phase for Dimension %d ---\n", dim_to_collapse))
    
    # 1. Define Start Orientation and Axes
    R_start_phase <- R_current_comp
    e_dim_local <- numeric(k); e_dim_local[dim_to_collapse] <- 1 # Hypercube's axis
    e_dim_world <- numeric(k); e_dim_world[dim_to_collapse] <- 1 # Target world axis
    
    current_axis_orientation_world <- R_start_phase %*% e_dim_local
    
    # 2. Calculate Total Relative Alignment Rotation needed for the phase
    cat("Calculating total alignment rotation...\n")
    R_align_target <- find_rotation_aligning(current_axis_orientation_world, e_dim_world)
    # find_rotation_aligning projects to SO(k)
    
    dist_to_target <- rotation_distance(R_start_phase %*% e_dim_local, e_dim_world) # Check initial angle
    cat(sprintf("Total alignment needed: %.4f radians\n", dist_to_target))
    
    phase_steps_i <- list()
    
    if (dist_to_target < 1e-7) {
      cat(sprintf("Alignment already met for dim %d. Adding trend steps.\n", dim_to_collapse))
      phase_steps_i <- rep(list(R_trend_step), num_collapse_steps_per_dim)
      R_product_actual_phase <- Reduce('%*%', phase_steps_i) # Calculate product of trend steps
    } else {
      cat("Generating composed steps...\n")
      # 3. Get Alignment Generator and Weights
      logm_attempt <- try(logm(R_align_target), silent=TRUE)
      if (inherits(logm_attempt, "try-error")) {
        warning("logm(R_align_target) failed for dim ", dim_to_collapse, ". Adding trend steps. Error: ", logm_attempt)
        phase_steps_i <- rep(list(R_trend_step), num_collapse_steps_per_dim)
      } else {
        A_align_total <- 0.5 * (logm_attempt - t(logm_attempt)) # Ensure skew-symmetric
        w_vec <- blending_weights(num_collapse_steps_per_dim, shape = comp_shape)
        
        # 4. Generate Steps via Composition
        for (i in 1:num_collapse_steps_per_dim) {
          w_prev <- if (i == 1) 0 else w_vec[i-1]
          w_curr <- w_vec[i]
          delta_w <- w_curr - w_prev
          # Safety check for non-monotonic weights (shouldn't happen with blend func)
          if(delta_w < -1e-10) warning("Negative delta_w encountered at step ", i)
          delta_w <- max(0, delta_w) # Ensure non-negative weight increment
          
          A_align_step_i <- delta_w * A_align_total
          expm_attempt <- try(expm(A_align_step_i), silent=TRUE)
          
          if(inherits(expm_attempt, "try-error")) {
            warning("expm(A_align_step) failed for step ", i, " dim ", dim_to_collapse, ". Using trend step only.")
            M_align_step_i <- diag(k)
          } else {
            M_align_step_i <- project_so_k(expm_attempt)
          }
          
          # Compose alignment portion with trend step
          M_i <- M_align_step_i %*% R_trend_step
          M_i <- project_so_k(M_i) # Ensure step is SO(k)
          
          phase_steps_i[[i]] <- M_i
          
        } # End loop through steps i
        
        # Calculate the product actually achieved by these steps for reporting/update
        if(length(phase_steps_i) > 0) {
          prod_calc <- try(Reduce("%*%", phase_steps_i), silent=TRUE)
          if(!inherits(prod_calc, "try-error")) {
            R_product_actual_phase <- project_so_k(prod_calc)
          } else {
            warning("Failed to calculate product of composed steps for phase ", dim_to_collapse)
            R_product_actual_phase <- diag(k) # Fallback product
          }
        } else {
          R_product_actual_phase <- diag(k)
        }
        
      } # End else (logm successful)
    } # End else (dist_to_target > tol)
    
    # Store the steps for this phase
    all_collapse_steps_comp <- c(all_collapse_steps_comp, phase_steps_i)
    
    # Update R_current for the next phase using the actual product achieved
    R_current_comp <- R_product_actual_phase %*% R_start_phase # Apply phase product
    R_current_comp <- project_so_k(R_current_comp)
    
    # Update the trend step for the next phase (use the last generated step)
    if (length(phase_steps_i) > 0) {
      R_trend_step <- project_so_k(phase_steps_i[[length(phase_steps_i)]])
    }
    # else: Keep the existing R_trend_step if phase failed/skipped/already aligned
    
  } # End loop over dimensions
} # End if dimensions needed collapsing

cat("\n--- Composition Collapse Sequence Generation Complete ---\n")

# --- Apply and Animate ---
if(length(full_hcubes)==0) stop("Need initial hypercube vertices.")
collapsed_hcube_vertices_comp <- full_hcubes[[length(full_hcubes)]]
collapse_history_comp <- list(collapsed_hcube_vertices_comp) # Start frame

cat("Applying Composition collapse steps...\n")
for (i in 1:length(all_collapse_steps_comp)) {
  M_step <- all_collapse_steps_comp[[i]]
  # Robustness check
  if (!is.matrix(M_step) || nrow(M_step) != k || ncol(M_step) != k || any(!is.finite(M_step))) {
    warning("Invalid Composition step matrix encountered at index ", i, ". Skipping step.")
    collapse_history_comp[[length(collapse_history_comp) + 1]] <- collapse_history_comp[[length(collapse_history_comp)]]
    next
  }
  prev_verts <- collapse_history_comp[[length(collapse_history_comp)]]
  current_verts <- t(M_step %*% t(prev_verts))
  collapse_history_comp[[length(collapse_history_comp) + 1]] <- current_verts
}
cat("Finished applying Composition collapse steps. Total frames:", length(collapse_history_comp), "\n")


# --- Verification (Composition Approach) ---
cat("\n--- Verifying Final Alignment after Composition ---\n")
final_verts_comp <- collapse_history_comp[[length(collapse_history_comp)]]
final_R_comp <- R_current_comp # Final cumulative rotation

# Check coordinates in the target world columns
for (d_check in start_dim_collapse:end_dim_collapse) {
  final_col_coords <- final_verts_comp[, d_check]
  cat(sprintf("\nCoordinates in world dimension Column %d (Target: +/- 0.5):\n", d_check))
  print(summary(final_col_coords))
  coord_range = range(final_col_coords)
  # Allow slightly more tolerance due to composition approximation
  if(abs(coord_range[1] - (-0.5)) < 0.05 && abs(coord_range[2] - 0.5) < 0.05) {
    cat(sprintf("Verification OK: Coordinates in column %d reasonably collapsed.\n", d_check))
  } else {
    cat(sprintf("Verification WARNING: Coordinates in column %d range [%.3f, %.3f] - alignment might be imperfect.\n", d_check, coord_range[1], coord_range[2]))
  }
}

# Verify the actual final rotation matrix alignment
axis_to_check_local <- numeric(k); axis_to_check_local[end_dim_collapse] <- 1
target_world_axis <- numeric(k); target_world_axis[end_dim_collapse] <- 1
final_axis_orientation_world <- final_R_comp %*% axis_to_check_local
final_dot_prod <- sum(final_axis_orientation_world * target_world_axis)
final_angle_diff_rad <- acos(max(-1, min(1, final_dot_prod)))

cat(sprintf("\nFinal orientation of hypercube axis e_%d relative to world axis w_%d:\n", end_dim_collapse, end_dim_collapse))
cat(sprintf("  Angle difference: %.4e radians (Target: 0)\n", final_angle_diff_rad))
if (final_angle_diff_rad < 1e-3) { # Allow slightly larger tolerance
  cat("Verification OK: Final rotation matrix achieves reasonable alignment.\n")
} else {
  cat("Verification WARNING: Final rotation matrix alignment deviates noticeably (angle: %.4e rad).\n", final_angle_diff_rad)
}

# Check the first step generated in the first collapse phase:
if(length(all_collapse_steps_comp) > 0) {
  first_collapse_step = all_collapse_steps_comp[[1]]
  original_trend_step = project_so_k(rot_mats[[length(rot_mats)]]) # Get initial trend again
  step_diff_norm = norm(first_collapse_step - original_trend_step, "F")
  cat(sprintf("\nDifference (Frobenius norm) between initial trend step and first collapse step: %.4e\n", step_diff_norm))
  # This should be small but likely non-zero, reflecting the start of alignment.
  if(step_diff_norm < 0.1) { # Arbitrary threshold for "close"
    cat("First collapse step appears reasonably close to the initial trend.\n")
  } else {
    cat("First collapse step seems potentially quite different from initial trend.\n")
  }
}