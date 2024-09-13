r <- 0.5
p <- 10

#find range numerically
missing.ind.values <- -100:100/100
R <- diag(p) + r - diag(p) * r 
missing.ind <- c(1,p)
dets <- log10(vapply(missing.ind.values, function(rij){
  R[missing.ind[1], missing.ind[2]] <- R[missing.ind[2], missing.ind[1]] <- rij
  det(R)
}, FUN.VALUE = numeric(1)))
plot(missing.ind.values, dets, type = "l", xlab = "r", ylab = "-log10(det(R))")
possible.values <- range(missing.ind.values[!is.na(dets)])

#find range analytically

# Load necessary library
library(Matrix)

# Function to estimate the missing entry using Schur complement
estimate_missing_entry <- function(R) {
  # Identify the location of the missing entry
  missing_indices <- which(is.na(R), arr.ind = TRUE)
  i <- missing_indices[1]
  j <- missing_indices[2]
  
  # Extract submatrices
  indices <- setdiff(1:nrow(R), c(i, j))
  R11 <- R[indices, indices]
  r1i <- R[indices, i]
  r1j <- R[indices, j]
  
  # Compute the Schur complement
  R11_inv <- solve(R11)
  schur_complement <- 1 - t(r1i) %*% R11_inv %*% r1j
  
  # Estimate the missing entry
  missing_entry <- -schur_complement / sqrt((1 - t(r1i) %*% R11_inv %*% r1i) * (1 - t(r1j) %*% R11_inv %*% r1j))
  
  return(missing_entry)
}

# Example correlation matrix with one missing entry
R <- matrix(c(
  1, 0.8, NA,
  0.8, 1, 0.5,
  NA, 0.5, 1
), nrow = 3, ncol = 3)

# Estimate the missing entry
missing_entry <- estimate_missing_entry(R)
missing_entry <- ifelse(missing_entry > 1, 1, ifelse(missing_entry < -1, -1, missing_entry))

print(paste("Estimated missing entry:", missing_entry))

# Fill the missing entry in the correlation matrix
R[is.na(R)] <- missing_entry

print("Updated correlation matrix:")
print(R)
