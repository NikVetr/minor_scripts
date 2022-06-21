dim <- 50
cov <- matrix(runif(n = 100, min = 0.1, max = 0.2), dim, dim); cov <- cov + t(cov); diag(cov) <- 1
cov
eigen(cov)$vectors %*% diag(eigen(cov)$values) %*% t(eigen(cov)$vectors)
solve(eigen(cov)$vectors)
plot(eigen(cov)$values)
cov

vec <- rnorm(n = 10, mean = 0, sd = 1)
vec <- vec / sum(vec^2)^0.5

vec2 <- rnorm(n = 9, mean = 0, sd = 1)
vec2[10] <- -1 * sum(vec2 * vec[-10]) / vec[10]
vec2 <- vec2 / sum(vec2^2)^0.5
sum(vec2 * vec)

vec3 <- rnorm(n = 8, mean = 0, sd = 1)
vec3[9:10] <- -1 * solve(matrix(data = c(vec[9:10], vec2[9:10]), nrow = 2, ncol = 2, byrow = T)) %*% c(sum(vec3 * vec[-c(9:10)]), sum(vec3 * vec2[-c(9:10)])) 
sum(vec * vec3)
sum(vec2 * vec3)

dim <- 10
eigvecs <- matrix(0, dim, dim)
for(i in 1:dim){
  vec <- rnorm(n = dim + 1 - i, mean = 0, sd = 1)
  if(i == 1){
  vec <- abs(vec)
  vec <- vec * rep(x = c(1,-1), length.out =(dim+1-i))
  vec <- -abs(vec)
  }
  if(i > 1){
    if (i == 2){
      vec[(dim+2-i):dim] <- -1 * solve(eigvecs[1:(i-1), (dim+2-i):dim]) %*% sum((vec * eigvecs[1:(i-1), -((dim+2-i):dim)])) 
    } else {
      vec[(dim+2-i):dim] <- -1 * solve(eigvecs[1:(i-1), (dim+2-i):dim]) %*% sapply(1:(i-1), function(x) sum(vec * eigvecs[x, -((dim+2-i):dim)]))
    }
  }
  vec <- vec / sum(vec^2)^0.5 #convert to unit vector
  eigvecs[i,] <- vec
}
eigvecs <- t(eigvecs)
round(eigvecs %*% t(eigvecs), 5)
eigvals <- dexp(x = sort(rexp(n = dim, rate = 10), decreasing = F), rate = 10)
eigvals[1] <- 50*eigvals[1]
cov <- eigvecs %*% diag(eigvals) %*% t(eigvecs)
cor <- cov2cor(cov)
cor
is.positive.definite(cov)
rmvnorm(100, sigma = cov)
