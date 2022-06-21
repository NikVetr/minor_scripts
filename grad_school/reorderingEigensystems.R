cov <- matrix(rnorm(26^2), 26)
cov <- cov %*% t(cov)
vecs <- eigen(cov)$vectors
vals <- eigen(cov)$values
rownames(cov) <- colnames(cov) <- LETTERS
cov
round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 

reords <- sample(1:dim(vecs)[2])
vecs  <- vecs[,reords]
vals <- vals[reords]
round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 

newOrd <- sample(LETTERS)
newCov <- cov[newOrd, newOrd]
newCov <- cov[newOrd,]
newCov <- newCov[,newOrd]
newCov == cov

per <- matrix(0, 26, 26)
for(i in 1:26){
  per[i,(27-i)] <- 1
}
eigen(per%*% cov %*% per)$values

per <- matrix(0, 26, 26)
newOrd <- sample(1:26)
for(i in 1:26){
  per[i,newOrd[i]] <- 1
}
eigen(t(per)%*% cov %*% per)$values
eigen((per)%*% cov %*% t(per))$values



newCov <- per %*% cov %*% t(per)
cov[newOrd, newOrd] - newCov #same matrix
newVecs <- eigen(newCov)$vectors
newVals <- eigen(newCov)$values

round((per) %*% vecs - newVecs, 5)
round(matrix(rep((per) %*% vecs[,3], 26)) - newVecs, 5)

newVecs

for(i in 1:1000){
Ls <- sample(LETTERS, size = 2, replace = T)
if(newCov[Ls[1], Ls[2]] == cov[Ls[1], Ls[2]]){cat("T")} 
}
newVecs <- eigen(newCov)$vectors
newVals <- eigen(newCov)$values
vals
abs(newVecs) - abs(vecs)

for(i in 1:26){
  vec <- vecs[,i]
  print(c(LETTERS[i],sum(abs(abs(vec) - abs(t(newVecs))) < 0.001)))
}
