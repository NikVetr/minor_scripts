n <- 1E2 #analogous to number of cell types
p <- 10 #analogous to number of metabolytes
r <- 0.7 #correlation among cell types
R <- diag(n) + r - diag(n) * r
L <- t(chol(R))
x <- matrix(rnorm(n*p), n, p)
x <- L %*% x #z-scores now multivariate normal with covariance matrix given by R
cx <- do.call(rbind, lapply(1:n, function(i) { #residualized scores
  fit <- lm(x[i,] ~ 1 + x[1,]) 
  resid <- x[i,] - x[1,] * fit$coefficients[2]
}))
apply(x, 1, var) - apply(cx, 1, var) #variance within each individual was reduced by this
head(cx) #first row had itself residualized out and is equal to 0 (up to machine precision)
head(cov(t(cx))) #it also has 0 covariance with other rows
mean(apply(x, 2, var) - apply(cx, 2, var)) #variance within each metabolyte not really reduced

#equivalent to partial correlation wrt 1st dimension 
# (all the prewritten functions don't let you specify conditioning dim afaict)
ppcor::pcor(t(x[1:3,]))$estimate[2:3,2:3]
cor(t(cx[2:3,]))

#can also compute this pseudo-partial correlation matrix directly
ucov <- cov(t(x))
head(cov(t(cx)[,-1])) #same as next line
head(ucov[-1, -1, drop = F] - t(ucov[1, -1, drop = F]) %*% solve(ucov[1, 1, drop = F]) %*% ucov[1, -1, drop = F])

#or the residuals directly
t(t(x) - t(t(x[1,])) %*% ucov[1,,drop=F] / ucov[1,1]) #same as next line
cx
