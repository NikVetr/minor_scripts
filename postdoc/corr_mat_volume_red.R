# Parameters for simulation
set.seed(123)
n <- 20
p <- 10
r_ij <- 0.5

# Create covariance matrix: 1s on diagonal, 0.5 on off-diagonals
Sigma <- matrix(r_ij, nrow = p, ncol = p)
diag(Sigma) <- 1.0

# Simulate multivariate normal data
x <- matrix(rnorm(n * p), p, n)
L <- t(chol(Sigma))

tx <- L %*% x
tx2 <- L%*% L %*% x
det(x %*% t(x)) * det(Sigma)
det(tx %*% t(tx))
det(x %*% t(x)) * det(Sigma) * det(Sigma)
det(tx2 %*% t(tx2))

library(expm)  # Load the package
L %^% 2 - L %*% L
nrep <- 10
Lrep <- L %^% nrep
apply(L100, 2, mean)
cov2cor(Lrep %*% t(Lrep))
det(Lrep %*% t(Lrep))
