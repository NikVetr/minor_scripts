n <- 1E3
d <- 1E1
x <- matrix(rnorm(n*d), n, d)
cor(x)
rij <- 0.7
r <- matrix(c(1, rij, rij, 1), 2, 2)
L <- t(chol(r))

pairs <- t(combn(d, 2))[1:2,]
for(i in 1:nrow(pairs)){
  x[,pairs[i,]] <- t(L %*% t(x[,pairs[i,]]))
}

#need to add in the partial correlations, not the marginal
#because earlier transformations already align pairs
cor(x[,c(2,3)])
