#### functions ####
random_orthonormal_matrix <- function(k, v = NULL, put_on_end = TRUE) {
  # Suppose 'v' is an existing k x q orthonormal basis (if provided).
  if (!is.null(v)) {
    if (nrow(v) != k) stop("The input matrix v must have k rows.")
    q <- ncol(v)
  } else {
    q <- 0
  }
  
  V <- matrix(0, k, k)
  if (q > 0) V[, 1:q] <- v
  
  # Build remaining columns via Gram–Schmidt
  for (i in (q+1):k) {
    repeat {
      v_temp <- rnorm(k)  # random vector, *not* normalized yet
      
      # Orthogonalize against existing columns
      if (i > 1) {
        for (j in 1:(i-1)) {
          v_temp <- v_temp - sum(v_temp * V[, j]) * V[, j]
        }
      }
      
      # Check that v_temp has nontrivial norm
      nrm <- sqrt(sum(v_temp^2))
      if (nrm > 1e-8) {
        # Normalize and accept
        v_temp <- v_temp / nrm
        V[, i] <- v_temp
        break
      }
    }
  }
  
  # Optionally put the existing basis vectors at the end
  if (q > 0 & put_on_end) {
    col_order <- c((q+1):k, 1:q)
  } else {
    col_order <- 1:k
  }
  V[, col_order]
}

unitvec <- function(p){
  vec <- rnorm(p)
  return(vec / norm(vec, "2"))
}

rlkj <- function (K, eta = 1) {
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
  sb2_sq <- (L_dropped - 1) / (L_dropped * sts - stes)
  sb2 <- sqrt(sb2_sq[sb2_sq > 0])
  s <- min(c(sa, sb1, sb2))
  # remainder_scale <- (1 - sb2^2*sts) / (1 - sb2^2*stes)
  # max(L_dropped) * remainder_scale #>1
  # min(L_dropped) * remainder_scale #< 1
  # s <- solve_for_constraint(theta, L_k, L_dropped)
  
  # Sample a Beta-distributed variable
  trunc_rnorm <- function(n, mean = 0, sd = 1, lb = -Inf, ub = Inf) {
    alpha <- (lb - mean) / sd
    beta <- (ub - mean) / sd
    p_lb <- pnorm(alpha)
    p_ub <- pnorm(beta)
    u <- runif(n, p_lb, p_ub)
    samples <- qnorm(u)
    samples <- mean + sd * samples
    return(samples)
  }
  
  y <- rbeta(1, shape1 = sum(L_k) / sum(L_dropped), shape2 = 1)^(1/k)
  yqs <- runif(1, -1, 1)
  yq <- abs(yqs)
  # y <- runif(1,0,1)^(1/k)
  # y <- abs(trunc_rnorm(1, 0, 1, lb = -1, ub = 1)) #hmm so this works poorly...
  y <- rbeta(1, shape1 = sum(L_k^2), shape2 = sum(L_dropped^2))^(1/k)
  
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
  x[drop_indices[1]] <- sqrt(x_d[1])
  x[drop_indices[2]] <- sqrt(x_d[2]) * sign(yqs)
  
  return(x)
}



#so now in this new basis defined by V, the ellipse equation is 
# D[1] * a^2 + D[2] * b^2 ... = 1
# and the unit hypersphere is still a^2 + b^2 ... = 1
#so we can multiply the latter by D[1] and subtract it from the former to get
#coordinates for the intersection less one dimension

# #TODO need additional constraint here!
# #not actually sampled from the intersection?
# #as sum(vec^2) can be > 1
# #just hack it with rejection sampling for now
# vec <- rep(10, prem)
# coefs <- (D[-1] - D[1]) / (1-D[1]) # = 1
# while(sum(vec^2) > 1){
#   vec <- unitvec(prem-1) / sqrt(coefs)   
# }
# evec <- c(sqrt(1 - sum(vec^2)), vec)


intersection_constraint <- function(L, V_curr, B_curr = NULL){
  
  ptot <- length(L)
  pcurr <- ncol(V_curr)
  prem <- ptot - pcurr
  
  #using a random matrix
  if(pcurr == 0){
    B <- diag(p)
  } else {
    if(is.null(B_curr)){
      #using QR decomp
      tmp_qr <- qr(V_curr)
      U_full <- qr.Q(tmp_qr, complete=TRUE)
      B <- U_full[ , (pcurr+1):ptot, drop=FALSE]
      # U <- cbind(B,
      #            U_full[ , 1:pcurr,      drop=FALSE])  
    } else {
      B <- B_curr
    }
  }
  
  # E <- diag(L)
  # Q <- t(U) %*% E %*% U 
  # Q <- Q[1:prem, 1:prem]
  
  #equivalently
  Q <- t(B) %*% diag(L) %*% B
  
  # Perform eigendecomposition of Q
  eig <- eigen(Q)
  QV <- eig$vectors  # Eigenvectors
  QD <- eig$values  # Eigenvalues as a diagonal matrix
  
  #ensure top eigenvector element is positive (to avoid teleporting)
  if(pcurr != 0){
    QV <- t(t(QV) * sign(QV[1,]))  
  }
  
  #sample this directly
  evec <- sample_intersection(L = QD)
  
  #so evec is our vector in the ellipsoid basis. We can transform it back to
  #the anti-plane basis with
  Qvec <- QV %*% evec
  
  #so this point has this coordinate in the old new basis U
  # new_pt <- c(Qvec, rep(0, pcurr)) #last elements set to 0 bc of how we defined the basis
  #and we can transform it back into the original coordinate system with
  # new_vec <- c(U %*% new_pt)
  
  #equivalently
  new_vec <- c(B %*% Qvec)
  # print(sum(new_vec^2))
  # print(sum(new_vec^2*L))
  
  #rank-1 update to obtain new complement (ie, B for the new V_curr)
  alpha <- crossprod(B, new_vec)
  V_new <- cbind(V_curr, new_vec)
  #remove any v-component from old complement columns
  B_new <- B - new_vec %*% t(alpha)
  # re-orthonormalize and shrink dimension by 1
  tmp_qr <- qr(B_new)
  # B_new <- qr.Q(tmp_qr, complete=TRUE)[ , 1:(prem-1), drop=FALSE]
  #slightly faster
  B_temp <- qr.Q(tmp_qr)  # dimension n x (n-m)
  B_new  <- B_temp[, 1:(prem-1), drop=FALSE]
    
  return(list(v = new_vec,
              B_new = B_new))
  
}
#switch up signs? if reusing old matrix. Getting lots more + correlations otherwise,
#probs because of how we generate the 'new_vec's with solving for the strictly + final elements
# sign_matrix <- ifelse2(prem > 2, diag(sample(c(-1,1), prem-1, T)), sample(c(-1,1), prem-1, T))
# B_new = B_new %*% sign_matrix


ifelse2 <- function(test, yes, no){
  if(test){
    return(yes)
  } else {
    return(no)
  }
}

sample_corrmat <- function(L, fresh_QR = F){
  p <- length(L)
  V <- diag(p) * 0
  #first row vector
  # u1 <- unitvec(p-1)
  # c1 <- (L[1] - L[-1]) / (L[1] - 1)
  # v1 <- u1 / sqrt(c1)
  # V[,1] <- c(sqrt(1 - sum(v1^2)), v1)
  intersect_out <- intersection_constraint(L, V_curr = matrix(0, p, 0))
  V[,1] <- intersect_out$v
  B_curr <- intersect_out$B_new
  # norm(V[,1], "2")
  # sum(V[,1]^2*L)
  
  #now for the rest
  for(i in 2:(p-1)){
    intersect_out <- intersection_constraint(L, V_curr = V[,1:(i-1), drop = F], 
                                             B_curr = ifelse2(fresh_QR, NULL, B_curr))
    V[,i] <- intersect_out$v
    B_curr <- intersect_out$B_new
    # print(norm(V[,i], "2"))
    # print(sum(V[,i]^2*L))
  }
  V[,p] <- B_curr
  eV <- t(V)
  
  #hyperellipsoid constraint
  # rowSums(eV^2)
  # colSums(eV^2)
  # rowSums(eV^2 %*% diag(L))
  # sum(eV[1,]^2 * diag(L))
  # max(abs((eV %*% t(eV))[upper.tri(diag(p))]))
  # max(abs((t(eV) %*% eV)[upper.tri(diag(p))]))
  
  #reconstruct correlation matrix
  C_samp <- eV %*% diag(L) %*% t(eV)
  # diag(C_samp)
  # range(C_samp[upper.tri(C_samp)])
  # eigen(C_samp)$values - L
  mean(abs(eigen(C_samp)$values - L) / L * 100)
  # plot(eigen(C_samp)$values, L)
  return(C_samp)
}


#### confirm properties ####
if(F){

p <- 10
R <- rlkj(p)
eigensystem <- eigen(R)
L <- diag(eigensystem$values)
V <- eigensystem$vectors

#unit hypersphere constraint on columns
colSums(V^2)

#hyperellipsoid constraint
rowSums(V^2)
rowSums(V^2 %*% L)
sum(V[1,]^2 * diag(L))

#### generation ####

#now to generate an eigenvector matrix from eigenvalues
p <- 4
L <- c(2.1,0.9,0.75,0.25)
L <- sort(unitvec(p)^2, T) * p

#first row vector
set.seed(1)
u1 <- unitvec(3)
c1 <- (L[1] - L[-1]) / (L[1] - 1)

v1 <- c(0, u1 / sqrt(c1))
v1[1] <- sqrt(1 - sum(v1^2))

#check
norm(v1, "2")
sum(v1^2*L)

#second row vector
# v2 <- c(0, 0, unitvec(2) / sqrt(c1[2:3]) * runif(1, 0, 1)^(1/2))
# v2[2] <- -sqrt(1 - sum(v2[3:4]^2 * c1[2:3]))
# v2[1] <- -sum(v2[-1] * v1[-1]) / v1[1]
# v2
# 
# sum(v2 * v1)
# norm(v2, "2")


U <- random_orthonormal_matrix(4, t(t(v1)))

#this is a new basis where the third vector defines the row we already solved for
#so in the new basis, the second row will always have a 0 in the fourth position
#the unit hypersphere has the same equation as before
#the hyperellipsoid now has a quadratic equation given by the matrix
E <- diag(L)       # 4x4
Q <- t(U) %*% E %*% U   # 4x4 matrix, represents equation for the hyperellipsoid as
# {a, b, c, d} Q {a, b, c, d}^T, where abcd are coordinates in the new basis
#and d = 0 to enforce the hyperplane constraint
Q <- Q[1:3, 1:3]

# Perform eigendecomposition of Q
eig <- eigen(Q)
V <- eig$vectors  # Eigenvectors
D <- eig$values  # Eigenvalues as a diagonal matrix

#so now in this new basis defined by V, the hyperellipsoid equation is 
# D[1] * a^2 + D[2] * b^2 + D[3] * c^2 = 1
# and the unit hypersphere is still a^2 + b^2 + c^2 = 1
#so we can multiply the latter by D[1] and subtract it from the former to get
bc_coefs <- (D[-1] - D[1]) / (1-D[1]) # = 1
bcvec <- unitvec(2) / sqrt(bc_coefs)
abcvec <- c(sqrt(1 - sum(bcvec^2)), bcvec)
norm(abcvec, "2")
t(abcvec) %*% diag(D) %*% t(t(abcvec))

#so abcvec is our vector in the hyperellipsoid basis. We can transform it back to
#the anti-plane basis with
Qvec <- V %*% abcvec
norm(Qvec, "2")
t(Qvec) %*% Q %*% Qvec

#so this point has this coordinate in the old new basis U
nv1_pt <- c(Qvec, 0) #last element set to 0 bc of how we defined the basis

#and we can transform it back into the original coordinate system with
v2 <- c(U %*% nv1_pt)
sum(v2 * v1)
norm(v2, "2")

v1
v2

#OK! So now we have two row vectors. We need to find a new point in the
#unit hypersphere and the hyperellipsoid that is orthogonal to both
#let's generate a new basis orthogonal to v1 and v2!

U <- random_orthonormal_matrix(4, cbind(v1, v2))
#and repeat the above procedure
Q <- t(U) %*% E %*% U
Q <- Q[1:2, 1:2]

# Perform eigendecomposition of Q
eig <- eigen(Q)
V <- eig$vectors  # Eigenvectors
D <- eig$values  # Eigenvalues as a diagonal matrix

#so now in this new basis defined by V, the ellipse equation is 
# D[1] * a^2 + D[2] * b^2 = 1
# and the unit hypersphere is still a^2 + b^2 = 1
#so we can multiply the latter by D[1] and subtract it from the former to get
c_coefs <- (D[-1] - D[1]) / (1-D[1]) # = 1
cvec <- unitvec(1) / sqrt(c_coefs)
bcvec <- c(sqrt(1 - sum(cvec^2)), cvec)
norm(bcvec, "2")
t(bcvec) %*% diag(D) %*% t(t(bcvec))

#so bcvec is our vector in the ellipsoid basis. We can transform it back to
#the anti-plane basis with
Qvec <- V %*% bcvec
norm(Qvec, "2")
t(Qvec) %*% Q %*% Qvec

#so this point has this coordinate in the old new basis U
nv1v2_pt <- c(Qvec, 0, 0) #last element set to 0 bc of how we defined the basis

#and we can transform it back into the original coordinate system with
v3 <- c(U %*% nv1v2_pt)
# v3 <- intersection_constraint(L, V_curr = cbind(v1, v2))
sum(v3 * v1)
sum(v2 * v3)
norm(v3, "2")

#so now we have vectors corresponding to the first three rows of of our eigenvector matrix
#the last one is full constrained by these
#and can be solved for by solving a linear system of equations
#because its dot product with all of them has to equal 0
#but that still gives us the choice of sign, doesn't it
#so maybe we still need the hyperellipsoid constraint?
#we can check both
V <- t(random_orthonormal_matrix(4, cbind(v1, v2, v3), F))

#unit hypersphere constraint on columns
colSums(V^2)

#hyperellipsoid constraint
rowSums(V^2)
rowSums(V^2 %*% L)
sum(V[1,]^2 * diag(L))

#reconstruct correlation matrix
C_samp <- V %*% diag(L) %*% t(V)
diag(C_samp)
eigen(C_samp)$values
L

}
#### now run it for real ####


source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")

par(mfcol = c(3,3))
p <- 10
for(i in 1:3){
  L <- sort(unitvec(p)^2, T)
  
  #block structure?
  # nb <- 3
  # prop_var <- 0.75
  # L <- c(rep(p * prop_var / nb,nb), rep(p * (1 - prop_var) / (p-nb), p-nb))
  # L <- L * runif(p, 1-1E-3, 1+1E-3) #dislodge equal values
  
  #near identity
  # L <- sort(abs(rnorm(p) + 10), T)
  
  #stairstep structure eigenvalues?
  # L <- exp(sort(cumsum(rbinom(p, 1, 0.1)) + abs(rnorm(p)) * 0.2, T))
  
  #resembling an existing matrix (with compound symmetry)
  # rij <- -1/(p+1)
  # rij <- 0.3
  # Rmat_true <- diag(p) + rij - diag(p) * rij
  # n <- p
  # Rmat <- cor(t(t(chol(Rmat_true)) %*% matrix(rnorm(n*p), p, n)))
  # L <- eigen(Rmat)$values
  
  L <- L / sum(L) * p
  cmat <- sample_corrmat(L, fresh_QR = F)
  print(range(diag(cmat)))
  print(range(L-eigen((cmat))$values))
  plot(L, type = "l")
  hist(cmat[upper.tri(cmat)], main = "", xlab = "pairwise correlation", breaks = -10:10/10)
  my_heatmap(cmat, plot_labels = F, dim_names = "", reorder_mat = T, col_scale = c(-1,1), 
             black_diag = T)
}

#### marginal distributions ####
#look at marginal distributions of correlations
p <- 10
nrep <- 2E3
cmats <- abind::abind(parallel::mclapply(1:nrep, function(i){
  L <- sort(unitvec(p)^2, T)
  L <- L / sum(L) * p
  cmat <- sample_corrmat(L, fresh_QR = F)
}, mc.cores = 8), along = 3)
par(mfrow = c(p,p), mar = c(1,1,1,1) * 0.5)
for(i in 1:p){
  for(j in 1:p){
    if(i <= j){
      plot.new()
    } else {
      hist(cmats[i,j,], breaks = -10:10/10, 
           main = paste0("(", i, ", ", j, ")"),
           xlab = "", ylab = "", yaxt = "n", 
           xaxt = "n", xpd=NA, cex.main = 0.9)
      segments(-1,0,1,0, xpd = NA)
      segments(-1,0,-1,par("usr")[4], xpd = NA)
      segments(par("usr")[2],0,par("usr")[2],par("usr")[4], 
               xpd = NA)
    }
  }
}
