n <- 10
ni <- 2000

x <- matrix(c(rnorm(n), rep(0, n*ni - n)), n, ni)
y <- matrix(c(rnorm(n), rep(0, n*ni - n)), n, ni) * 0
z <- matrix(c(rnorm(n), rep(0, n*ni - n)), n, ni) * 0

b <- runif(3, 0.2, 0.4)
for(i in 2:ni){
  y[,i] <- y[,i-1] + b[1] * (x[,i] - x[,i-1])
  z[,i] <- z[,i-1] + b[2] * (y[,i] - y[,i-1])
  x[,i] <- x[,i-1] + b[3] * (z[,i] - z[,i-1])
}

# pairs(cbind(x[,ni], y[,ni], z[,ni]))
par(mfrow = c(3,1))
for(j in list(x,y,z)){
  plot(j[2,], type = "l", ylim = range(j))
  for(i in 1:n){lines(j[i,])}  
}

