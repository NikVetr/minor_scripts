library(expm)

# 1) Helper: put the off-diagonal of G from a vector gamma
#    (length gamma = n(n-1)/2)
vecLower_to_matSym <- function(gamma, n) {
  G <- matrix(0, n, n)
  # Fill lower-triangle (excluding diagonal) with gamma
  # and mirror to keep G symmetric
  idxLower <- lower.tri(G, diag=FALSE)
  G[idxLower] <- gamma
  G <- G + t(G)
  return(G)
}

# 2) Core iterative routine to ensure diag(exp(G)) = 1
#    as proposed in Archakov & Hansen (2021a, 2022).
#    The function modifies diag(G) iteratively using the
#    contraction mapping Gii <- Gii - log(Cii).
ensure_diag_one <- function(G, tol = 1e-12, maxIter = 1000) {
  for (iter in seq_len(maxIter)) {
    C   <- expm(G)         # matrix exponential
    d   <- diag(C)         # current diagonal
    err <- max(abs(d - 1)) # how far from 1's?
    
    if (err < tol) break
    
    # Update Gii to fix the diagonal.  If Cii>1, we
    # subtract log(Cii) from Gii so next iteration
    # shrinks it, etc.
    diag(G) <- diag(G) - log(d)
  }
  return(G)
}

# 3) Full function: from gamma -> correlation matrix
C_of_gamma <- function(gamma, n, tol=1e-12, maxIter=1000) {
  # Step 1: Build G from gamma
  G <- vecLower_to_matSym(gamma, n)
  
  # Step 2: Solve for the diagonal so that diag(exp(G))=1
  G <- ensure_diag_one(G, tol=tol, maxIter=maxIter)
  
  # Step 3: Exponentiate to get the final correlation
  C <- expm(G)
  return(C)
}

# ------------------------------
# EXAMPLE USAGE
# ------------------------------
n <- 50

# We want d = n(n-1)/2 random draws in R^d
d <- n*(n-1)/2

# For demonstration, draw gamma from Normal(0,1)
# (That is purely arbitrary; you can pick any distribution.)
gamma_draw <- rnorm(d, mean=1, sd=0.1)
d_alt <- sample(1:d, size = floor(d/2))
gamma_draw[d_alt] <- rnorm(length(d_alt), mean=-3, sd=0.05)

# Build the correlation matrix
Cmat <- C_of_gamma(gamma_draw, n)
hist(Cmat[upper.tri(Cmat)], breaks = -10:10/10)