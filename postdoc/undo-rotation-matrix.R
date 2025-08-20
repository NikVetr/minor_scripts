generate_blended_rotation_steps <- function(start_coords, end_coords, R_start, R_target, nframes, easing_fn = NULL) {
  if (is.null(easing_fn)) {
    easing_fn <- function(x) sin(x * pi / 2)
  }
  
  # Generators
  A_start  <- logm(R_start)
  A_target <- logm(R_target)
  
  coords_seq   <- vector("list", nframes)
  cumrot_seq   <- vector("list", nframes)
  rotstep_seq  <- vector("list", nframes)
  
  R_cum <- diag(nrow(R_start))  # identity to start
  
  for (j in 1:nframes) {
    t <- (j - 1) / (nframes - 1)
    phi <- easing_fn(t)
    
    A_blend <- (1 - phi) * A_start + phi * A_target
    R_step <- expm(A_blend / nframes)  # per-frame rotation step
    
    R_cum <- R_step %*% R_cum
    
    coords_seq[[j]] <- t(R_cum %*% t(end_coords))
    cumrot_seq[[j]] <- R_cum
    rotstep_seq[[j]] <- R_step
  }
  
  return(list(
    coords     = coords_seq,
    cumrotmats = cumrot_seq,
    rotsteps   = rotstep_seq
  ))
}


start_idx <- dim_break_pts[4]
end_idx   <- dim_break_pts[5] - 1
subnrots  <- end_idx - start_idx + 1

start_coords <- full_hcubes[[start_idx]]
end_coords   <- full_hcubes[[dim_break_pts[5]]]
forward_rots <- rot_mats[start_idx:end_idx]

R_total      <- Reduce("%*%", forward_rots)
R_target <- solve(R_total)
R_start      <- rot_mats[[end_idx]]  # the last applied rotation in this segment


R_start <- rot_mats[[end_idx]]

# Get R_target via Procrustes
R_target <- {
  X <- full_hcubes[[dim_break_pts[5]]]
  Y <- full_hcubes[[start_idx]]
  svd_res <- svd(t(Y) %*% X)
  R <- svd_res$u %*% t(svd_res$v)
  if (det(R) < 0) {
    svd_res$u[, ncol(svd_res$u)] <- -svd_res$u[, ncol(svd_res$u)]
    R <- svd_res$u %*% t(svd_res$v)
  }
  R
}

res <- generate_blended_rotation_steps(
  start_coords = full_hcubes[[start_idx]],
  end_coords   = full_hcubes[[dim_break_pts[5]]],
  R_start      = R_start,
  R_target     = R_target,
  nframes      = subnrots
)

# Validate
final_coords <- tail(res$coords, 1)[[1]]
cat("Error from start_coords:", mean(abs(final_coords - start_coords)), "\n")

res$rotsteps[[1]]
res$rotsteps[[2]]
R_start


print(round(R_start - res$rotmats[[1]], 5))
print(round(res$rotmats[[1]] - res$rotmats[[2]], 5))
print(round(res$rotmats[[2]] - res$rotmats[[3]], 5))
