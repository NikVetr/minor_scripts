# set.seed(1)
p <- 80 #analogous to cell types
n <- 100 #analogous to genes or analytes
d <- matrix(rnorm(p*n), ncol = n*p, nrow = 1) #uncorrelated obs

#correlated ps?
r_p <- 0.7
r_n <- 0.5
R_p <- diag(p) + r_p - diag(p) * r_p
R_n <- diag(n) + r_n - diag(n) * r_n
# R <- kronecker(R_p, R_n)
RR <- kronecker(chol(R_p), chol(R_n))
d <- d %*% RR
d <- matrix(c(d), ncol = p, nrow = n)

#check PC score correlations
par(mfrow = c(4,1))
breaks = -50:50/50
hist(cor(t(prcomp(d)$x)), breaks = breaks)
hist(cor(t(prcomp(t(d))$x)), breaks = breaks)

#try consecutive PCA?
dpca_1 <- prcomp(t(prcomp(d)$x))$x
dpca_2 <- prcomp(t(prcomp(t(d))$x))$x

hist(cor(t(dpca_1)), breaks = breaks)
hist(cor(t(dpca_2)), breaks = breaks)

#try eigendecomp approach?
# cov(d %*% eigen(cov(d))$vectors %*% diag(1/sqrt(eigen(cov(d))$values)))
# diag(kronecker(diag(eigen(R_p)$values), diag(eigen(R_n)$values)))
R_p_2use <- R_p
R_n_2use <- R_n
R_p_2use <- as.matrix(Matrix::nearPD(cov(d))$mat)
R_n_2use <- as.matrix(Matrix::nearPD(cov(t(d)))$mat)
eigenvalues <- rep(eigen(R_p_2use)$values, each = n) * rep(eigen(R_n_2use)$values, p)
bad_L <- eigenvalues < 1E-6
eigenvalues[bad_L] <- Inf
eigenvectors <- kronecker(eigen(R_p_2use)$vectors, eigen(R_n_2use)$vectors)
d_eig <- matrix(d, ncol=p*n) %*% eigenvectors %*% diag(1/sqrt(eigenvalues))
# d_eig[, bad_L] <- 0
d_kron <- matrix(c(d_eig), ncol = p, nrow = n)
hist(cor(d_kron, use = "pair"), breaks = breaks)
hist(cor(t(d_kron), use = "pair"), breaks = breaks)

#reconstruct data?
recon_eigenvalues <- eigenvalues
recon_eigenvalues[bad_L] <- 0
d_eig_recon <- d_eig %*% diag(sqrt(recon_eigenvalues)) %*% t(eigenvectors)
d_eig_recon <- matrix(c(d_eig_recon), ncol = p, nrow = n)
par(mfrow = c(1,1))
plot(d_eig_recon, matrix(d, ncol=p*n)); abline(0,1,col=2,lwd=2)
hist(cor(d_eig_recon)[upper.tri(R_p)]); abline(v = r_p, col = "red", lwd = 3)
hist(cor(t(d_eig_recon))[upper.tri(R_n)]) ; abline(v = r_n, col = "red", lwd = 3)

#shuffle data indices and check if correlations are preserved
d_eig_recon <- sample(d_eig) %*% diag(sqrt(recon_eigenvalues)) %*% t(eigenvectors)
d_eig_recon <- matrix(c(d_eig_recon), ncol = p, nrow = n)
par(mfrow = c(1,1))
plot(d_eig_recon, matrix(d, ncol=p*n)); abline(0,1,col=2,lwd=2)
hist(cor(d_eig_recon)[upper.tri(R_p)]); abline(v = r_p, col = "red", lwd = 3)
hist(cor(t(d_eig_recon))[upper.tri(R_n)]) ; abline(v = r_n, col = "red", lwd = 3)


#try only single pca + permutation + reconstruction 
#to check if covariances are preserved
ed <- eigen(cov(d))
eval <- ed$values
evec <- ed$vectors
dp <- d %*% evec %*% diag(1/sqrt(eval))
dps <- do.call(cbind, lapply(1:p, function(i) dp[sample(1:n), i]))
dpr <- dps %*% diag(sqrt(eval)) %*% t(evec)

#compare the two
breaks = -20:20/20
hist(cor(d)[upper.tri(diag(p))], breaks = breaks, col = 2)
hist(cor(dpr)[upper.tri(diag(p))], add = T, 
     breaks = breaks, col = adjustcolor(1,0.25))
plot(log(apply(d,2,var)), log(apply(dpr,2,var))); abline(0,1)

hist(cor(d)[upper.tri(diag(p))], breaks = breaks, col = 2)
hist(cor(dpr)[upper.tri(diag(p))], add = T, 
     breaks = breaks, col = adjustcolor(1,0.25))


hist(cor(t(d))[upper.tri(diag(p))], breaks = breaks, col = 2)
hist(cor(t(dpr))[upper.tri(diag(p))], add = T, 
     breaks = breaks, col = adjustcolor(1,0.25))

plot(log(apply(d,2,var)), log(apply(dpr,2,var))); abline(0,1)

