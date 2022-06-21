checkmat <- function(d, corr){
  r <- diag(d) * (1-corr) + corr
  r[1,2] <- r[2,1] <-0
  all(eigen(r)$values >= 0)
}


checkmat <- function(d, corr){
  r <- diag(d) * (1-corr) + corr
  r[1,2] <- r[2,1] <-0
  det(r)
}

checkmat(3, 0.7)

d = 3
corr = sqrt(1/2) - sqrt(1/2)
r <- diag(d) * (1-corr) + corr
r[1,2] <- r[2,1] <- 0.9
all(eigen(r)$values >= 0)


rs_to_check <- 1:1000/1000
dims_to_check <- 3:100

out <- sapply(setNames(dims_to_check, dims_to_check), function(d) 
  max(rs_to_check[sapply(rs_to_check, function(corr) checkmat(d, corr))]))

plot(names(out), out, type = "l", xlab = "number of variables (aka dimensionality of correlation matrix)",
     ylab = "maximum correlation of all but that one zero correlation", cex.lab = 1.5, lwd =2, col = "darkblue")
box(lwd=2)


lkj1 <- function(){
  d = 3
  eta = 1
  Beta = eta + (d-2)/2
  r12 <- 2*rbeta(1,Beta,Beta)-1
  r <- matrix(c(1,r12,r12,1),2,2)
  Beta <- Beta - 1/2
  y <- rbeta(1, 2, Beta)
  u <- rnorm(2)
  u <- u / sqrt(sum(u^2))
  w <- sqrt(y) * u
  z <- t(chol(r)) %*% w
  r <- rbind(cbind(r, z), cbind(t(z), 1))
  r  
}

x <- replicate(1000000, lkj1())
hist(x[2,3,x[1,2,] > 0.5 & x[1,2,] < 0.6 & x[1,3,] > 0.3 & x[1,3,] < 0.4])


r
detr <- function(r){
  1 - r[2,3]^2 - r[1,2]^2 + r[2,3]*r[1,3]*r[1,2] + r[1,2]*r[1,3]*r[2,3] - r[1,3]^2
}

detr <- function(r12, r13, r23){
  1 - r23^2 - r12^2 + r23*r13*r12 + r12*r13*r23 - r13^2
}

corrs <- -1000:1000/1000

dens <- detr(0.5, -0.25, corrs)
dens[dens < 0] <- 0
eta = 2
dens <- dens^(eta-1)
dens <- dens / sum(diff(corrs) * (dens[-length(dens)] + diff(dens)/2))
plot(corrs, dens, type = "l")

a = -0.15
b = matrix(c(1,0.6,0.7,0.6,1,a,0.7,a,1),3,3)
eigen(b)
