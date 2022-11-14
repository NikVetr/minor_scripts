p <- 3
n <- 2^(p+1)

test_assignment <- function(p){
  
  bc <- rnorm(p, 1, 1)
  
  x <- as.matrix(expand.grid(rep(list(0:1), p+1)))
  xc <- x[,-1]
  x1f <- x[,1] #balanced across subcats
  x2f <- sample(x1f) #random across subcats
  x3f <- xc_max <- xc[,which.max(abs(bc))]
  x3f[xc_max == 0] <- sample(rep(c(0,1), 2^(p-1)))
  x3f[xc_max == 1] <- sample(rep(c(0,1), 2^(p-1)))
  
  
  yc <- xc %*% t(t(bc))
  bf <- 1
  e <- rnorm(n, sd = sqrt(var(yc) * 9))
  ye <- yc + e
  y1e <- x1f * bf + ye
  y2e <- x2f * bf + ye
  y3e <- x3f * bf + ye
  
  return(c(balanced = summary(lm(y1e ~ cbind(x1f, xc)))$coefficients[2,4],
           random = summary(lm(y2e ~ cbind(x2f, xc)))$coefficients[2,4],
           targeted = summary(lm(y3e ~ cbind(x3f, xc)))$coefficients[2,4]))
}

nrep <- 1E3
out <- t(replicate(nrep, test_assignment(p)))
pairs(log(out))
hist(apply(log(out[,c(1,3)]), 1, diff)) #positive values indicate lower pvals in balanced group
mean(apply(log(out[,c(1,3)]), 1, diff) > 0)

