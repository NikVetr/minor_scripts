source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")

compute_column_independence <- function(R_mat) {
  qr_decomp <- qr(R_mat)
  R_upper <- qr.R(qr_decomp)  # Extract R from QR decomposition
  independence_scores <- abs(diag(R_upper))  # Larger = more independent
  return(independence_scores / sum(independence_scores))  # Normalize
}

rlkj_orig <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2){
    for (m in 2:(K - 1)) {
      y <- rbeta(1, 0.02, 0.05)^(m) #either going off in a new direction (y near 0), or alligned in with the old (y near)
      # z <- R[1:m, m]# + rnorm(m, 0, 1)
      # prec <- compute_column_independence(R[1:m, 1:m])
      # z <- R[1:m, 1:m] %*% t(gtools::rdirichlet(1, alpha = prec^2)) + rnorm(m, 0, 0.2)
      
      # z <- rowMeans(R[1:m, 1:m]) + rnorm(m, 0, 0.2)
      # z <- R[1:m, 1:m] %*% t(t(rep(1,m))) + rnorm(m, 0, 1)
      z <- R[1:m, 1:m] %*% t(t(colMeans(abs(crossprod(R[1:m, 1:m]))^(1/2) - diag(m)))) + rnorm(m, 0, 1)
      
      z <- z/sqrt(crossprod(z)[1]) #rescale to unit length
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  return(R)
}
par(mfcol = c(2,3))


rcorr <- function (K, eta = 1) {
  R <- matrix(0, K, K)
  R[1,1] <- 1
  z <- rnorm(K)
  z <- z / sqrt(crossprod(z)[1])
  R[,K] <- z
  for (m in (K - 1):2) {
    y <- rbeta(1, 0.2, 0.5)
    z <- rnorm(m-1) * 0.5 + R[1:(m-1),K]
    z <- z / sqrt(crossprod(z)[1])
    R[1:(m-1), m] <- sqrt(y) * z
    R[m + 1, m + 1] <- sqrt(1 - y) * sign(R[m, m+1])
  }
  return(R)
}

par(mfcol = c(2,3))
for(i in 1:3){
  R <- rcorr(50)
  cmat <- crossprod(R)
  hist(cmat[upper.tri(cmat)], main = "", xlab = "pairwise correlation")
  my_heatmap(cmat, plot_labels = F, dim_names = "", reorder_mat = T)
}
mat_order <- order(cmdscale(1.1 * max(abs(R)) - R, k = 1))
# plot(mat_order, order(zK))



# rlkj <- function (K, eta = 1) {
#   alpha <- eta + (K - 2)/2
#   r12 <- 2 * rbeta(1, alpha, alpha) - 1
#   R <- matrix(0, K, K)
#   R[1, 1] <- 1
#   R[1, 2] <- r12
#   R[2, 2] <- sqrt(1 - r12^2)
#   zK <- rbinom(m, 2, prob = 0.5)
#   z <- rnorm(1)
#   if (K > 2) 
#     for (m in 2:(K - 1)) {
#       alpha <- alpha - 0.5
#       y <- rbeta(1, m/2, alpha)
#       # z <- rnorm(m, 0, 1)
#       z <- zK[1:m] + rnorm(m)
#       # z <- rbinom(m, 3, prob = 0.5)
#       z <- z/sqrt(crossprod(z)[1])
#       R[1:m, m + 1] <- sqrt(y) * z
#       R[m + 1, m + 1] <- sqrt(1 - y)
#     }
#   return(list(zK = zK, R = R))
# }
# 
# R_out <- rlkj(50)
# zK <- R_out$zK
# R <- R_out$R
# cmat <- crossprod(R)
# flip_axes <- diag(rbinom(nrow(cmat), 1, 0.5) * 2 - 1)
# cmat <- flip_axes %*% cmat %*% t(flip_axes)
