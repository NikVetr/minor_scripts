#perform this step once
n <- 1E5

#initialize to something
x0 <- seq(0,1,length.out=n)
x0 <- runif(n)

r <- rank(x0)
q <- pbeta(x0, r, n + 1 - r)

#iterate
niter <- 0
if(niter > 0){
  
  for(i in 1:niter){
    x <- q
    r <- rank(x)
    q <- pbeta(x, r, n + 1 - r)
  }
    
}

par(mfrow = c(2,1), mar = c(4,4,2,1))
breaks <- 0:100/100
hist(x0, breaks = breaks)
hist(q, breaks = breaks)
mean(x)
mean(q)


#try replicating?

nrep <- 1E4
n <- 1E2
q <- replicate(nrep, {
  x0 <- runif(n)
  r <- rank(x0)
  q <- pbeta(x0, r, n + 1 - r)
  q
})
par(mfrow = c(2,1), mar = c(4,4,2,1))
breaks <- 0:100/100
hist(x0, breaks = breaks)
hist(q, breaks = breaks)
