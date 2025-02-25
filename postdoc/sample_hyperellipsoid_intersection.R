sample_intersection <- function(L) {
  # Reduce dimensionality
  p <- length(L)
  k <- p - 2
  
  # Select two indices to drop
  drop_indices <- c(which.min(L), which.max(L))
  # drop_indices <- c(which(L > 1)[which.min(L[L > 1])], 
  #                   which(L < 1)[which.max(L[L < 1])])
  remaining_indices <- setdiff(1:p, drop_indices)
  
  # Extract the reduced coefficients
  L_k <- L[remaining_indices]
  L_dropped <- L[drop_indices]
  
  # Sample a random direction from the unit hypersphere in k dimensions
  theta <- rnorm(k)
  # theta <- theta / sqrt(sum(theta^2)) / sqrt(L_k * (L_k > 1) * 1 + (L_k < 1))
  # theta <- theta / sqrt(sum(theta^2)) / sqrt((L_k + 1) / 2)
  # theta <- theta / sqrt(sum(theta^2)) / L_k^2
  theta <- theta / sqrt(sum(theta^2))
  
  # Compute scaling coefficients with constraint
  sts <- sum((theta)^2)
  stes <- sum(L_k * (theta^2))
  sa <- 1 / sqrt(sts)
  sb1 <- 1 / sqrt(stes)
  #can compute it numerically for now
  # sb2_poss <- seq(0, min(sa, sb1), length.out = 200)
  # sb2s <- sapply(sb2_poss, function(i){
  #   remainder_scale <- (1 - i^2*sts) / (1 - i^2*stes)
  #   satisfies_conditions <- max(L_dropped) * remainder_scale > 1 &
  #     min(L_dropped) * remainder_scale < 1
  #   return(satisfies_conditions)
  # })
  # sb2 <- sb2_poss[sum(sb2s)]
  #or solve for it algebraically
  sb2 <- sqrt((L_dropped - 1) / (L_dropped * sts - stes))
  s <- min(c(sa, sb1, sb2))
  # remainder_scale <- (1 - sb2_est^2*sts) / (1 - sb2_est^2*stes)
  # max(L_dropped) * remainder_scale #>1
  # min(L_dropped) * remainder_scale #< 1
  # # min(L_dropped) * A
  # # max(L_dropped) * A
  # s
  # min(sa, sb1, sb2_est)
  # s <- solve_for_constraint(theta, L_k, L_dropped)
  
  # Sample a Beta-distributed variable
  # y <- rbeta(1, shape1 = sum(L_k), shape2 = sum(L_dropped))^(1/k)
  y <- rbeta(1, shape1 = sum(L_k) / sum(L_dropped), shape2 = 1)^(1/k)
  
  # Scale the sample
  x_k <- theta * y * s
  
  # Solve for the two dropped coordinates (x_d1, x_d2)
  # Define the system of equations
  A <- matrix(c(1, 1, L[drop_indices[1]], L[drop_indices[2]]), nrow = 2, byrow = TRUE)
  b <- c(1 - sum(x_k^2), 1 - sum(L_k * x_k^2))
  
  # Solve for the two unknowns (ensuring real solutions exist)
  x_d <- solve(A, b)
  
  # Construct the final sampled point
  x <- numeric(p)
  x[remaining_indices] <- x_k
  x[drop_indices] <- sqrt(x_d) * sample(c(-1,1), 2, T)
  
  return(x)
}

# Example usage
set.seed(1)
p <- 100
L <- rnorm(p)^2  # Example hyperellipsoid coefficients
L <- sort(L / sum(L) * p, T)
nrep <- 1E3
xs <- do.call(rbind, parallel::mclapply(1:nrep, function(repi) 
  sample_intersection(L), mc.cores = 12))
marginal_sds <- apply(xs, 2, sd)
marginal_means <- apply(xs, 2, function(xi) mean(abs(xi)))

par(mfrow = c(2,1), mar = c(4,4,2,2))
plot(marginal_sds, ylim = c(0, max(marginal_sds)),
     xlab = "index of element in vector", ylab = "standard deviation",
     main = c("Standard Deviation of Sampled Elements"))
plot(marginal_means, ylim = c(0, max(marginal_means)),
     xlab = "index of element in vector", ylab = "mean absolute value",
     main = c("Mean Absolute Value of Sampled Elements"))
apply(xs, 1, function(xi) sum(xi^2))
apply(xs, 1, function(xi) sum(xi^2 * L))
