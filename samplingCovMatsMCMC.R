n <- 20
cov <- matrix(rnorm(n^2), n, n); cov <- cov %*% t(cov)
cor <- cov2cor(cov)
cor
proposed.elements <- rep(0, dim(cor)[1])
cor <- cbind(rbind(cor, t(proposed.elements)), c(proposed.elements,1)) #initialize
any(eigen(cor)$values < 0)
cor

niter <- 1e7
thin <- 1e3
correlation.coefficients <- matrix(0, ncol = dim(cor)[1], nrow = niter/thin)
for(i in 1:niter){
  samp.indices <- sample(1:(dim(cor)[1]-1), size = sample(dim(cov)[1] - 1))
  samp.width <- runif(1,0,0.1)
  displacements <- runif(length(samp.indices), -1*sampWidth, sampWidth)
  proposed.elements <- cor[dim(cor)[1],]
  proposed.elements[samp.indices] <- proposed.elements[samp.indices] + displacements
  proposed.elements[proposed.elements > 1] <- proposed.elements[proposed.elements > 1] * -1 + 2
  proposed.elements[proposed.elements < -1] <- proposed.elements[proposed.elements < -1] * -1 - 2
  propCor <- cor
  propCor[dim(cor)[1],] <- propCor[,dim(cor)[1]] <- proposed.elements
  if(det(propCor) >= 0){ #auto accept valid matrices, since ratios all = 1
    cor <- propCor
  }
  if (i %% thin == 0){
    print(paste0((i / niter) * 100, "% complete"))
    print(cor[dim(cor)[1],])
    correlation.coefficients[i/thin,] <- cor[dim(cor)[1],]
  }
}
any(eigen(cor)$values < 0)

par(mfrow = c(1,2))
for(k in 1:n){
  plot(correlation.coefficients[,k], type = "l")
}

