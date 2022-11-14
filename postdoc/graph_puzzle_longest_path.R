library(igraph)

`%^^%` <- function(A, k) { #from https://eric.netlify.app/2017/08/08/taking-powers-of-a-matrix-in-r/
  eig <- eigen(A)
  stopifnot(length(unique(eig$values)) == nrow(A))
  P <- eig$vectors
  D <- diag(eig$values ^ k)
  Ak <- P %*% D %*% solve(P)
  round(Re(Ak))
}

N <- 10
K <- 2

mat <- sample_k_regular(no.of.nodes=N, k=K)
adj <- as.matrix(as_adj(mat))



adj %^^% 10
