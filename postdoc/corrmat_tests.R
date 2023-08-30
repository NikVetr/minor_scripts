library(corpcor)


vals <- -1000:1000/1000
out <- sapply(vals, function(b){
  x <- diag(3) * 0.2 + 0.8
  a <- 0.5
  c <- 0.9
  vec <- c(a,b,c)
  x <- cbind(rbind(x, c(vec)), c(vec, 1))
  all(eigen(x)$values >= 0)
})

range(vals[out])


k <- 10
x <- diag(k) * 0.2 + 0.8
n <- 250
d <- matrix(rnorm(n * k), k, n)
d <- apply(d, 1, scale)
y <- t(t(chol(x)) %*% t(d %*% solve(chol(cov(d)))))
x
cov(y)

y

d2rot <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
len <- function(x) sqrt(sum(x^2))

dkrot <- function(k, thetas) {
  angles <- data.frame(cbind(t(combn(1:k, 2)), theta = thetas))
  rotmats <- lapply(angles$theta, d2rot)
  big_rotmats <- lapply(1:choose(k, 2), function(i){
    bigrm <- diag(k)
    bigrm[angles$V1[i], angles$V1[i]] <- rotmats[[i]][1,1]
    bigrm[angles$V1[i], angles$V2[i]] <- rotmats[[i]][1,2]
    bigrm[angles$V2[i], angles$V1[i]] <- rotmats[[i]][2,1]
    bigrm[angles$V2[i], angles$V2[i]] <- rotmats[[i]][2,2]
    bigrm
  })
  rotmat <- Reduce("%*%", big_rotmats)
  rotmat
}

rotmat <- dkrot(k, c(0.3, rep(0, choose(k,2)-1)))

y <- y %*% rotmat
cor(y)

x
px <- cor2pcor(x)
px[1,2] <- px[2,1] <- -0.9
pcor2cor(px)
