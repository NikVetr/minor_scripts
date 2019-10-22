cov <- matrix(rnorm(26^2), 26)
cov <- cov %*% t(cov)
vecs <- eigen(cov)$vectors
vals <- eigen(cov)$values
rownames(cov) <- colnames(cov) <- LETTERS
cov
round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 

#reordering eigenvectors and eigenvalues results in same composition to cov
# reords <- sample(1:dim(vecs)[2])
# vecs  <- vecs[,reords]
# vals <- vals[reords]
# round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 

#construct a permutation matrix
per <- matrix(0, 26, 26)
newOrd <- sample(1:26)
for(i in 1:26){
  per[i,newOrd[i]] <- 1
}
#verify integrity of eigenvalues under both alternative permutations (these are not the same permutation)
eigen(t(per)%*% cov %*% per)$values
eigen((per)%*% cov %*% t(per))$values

#this permutation is equivalent to re-indexing in the usual way (see next line)
newCov <- per %*% cov %*% t(per)
cov[newOrd, newOrd] - newCov #same matrix
newVecs <- eigen(newCov)$vectors
newVecs2 <- (per) %*% vecs
newVals <- eigen(newCov)$values

#this gives many but not all of the same eigenvectors
round((per) %*% vecs - newVecs, 5)

#but wait, are they the right eigenvectors?

#they're orthogonal
round(newVecs2 %*% t(newVecs2) - diag(rep(1,26)), 5)

#and length  one
apply(X = ((per) %*% vecs) * ((per) %*% vecs), MARGIN = 2, FUN = sum)

#and appropriately recompose the permuted covariance matrix
round(newVecs2 %*% diag(newVals) %*% t(newVecs2) - newCov, digits = 8) 

#as do the directly computer eigenvectors
round(newVecs %*% diag(newVals) %*% t(newVecs) - newCov, digits = 8) 

#so why aren't they equal??
round(newVecs2 - newVecs, 5)

#wait are they just reversed in sign?
round(abs(newVecs2) - abs(newVecs), 5)

#totally!
(newVecs2 - newVecs) / newVecs2
