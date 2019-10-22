library(mvtnorm)

meandiff <- vector()
for(i in 1:1000){
  print(i)
  simD <- rmvnorm(3,sigma = cov)
  meandiff[i] <- mean(cov(simD) - cov)
}
hist(meandiff)


c <- matrix(c(1,.5,.7,.5,1,.6,.7,.6,1), nrow = 3)


niter <- 1e5
df <- 300
wishes <- list()
for(i in 1:niter){
  print(i)
  wishes[[i]] <- rwish(df, diag(3))
}

wishez <- list()
for(i in 1:niter){
  print(i)
  simdata <- rmvnorm(n = df, sigma = diag(3))
  wishez[[i]] <- cov(simdata) * df
}


var1 <- sapply(1:length(wishes), function (x) wishes[[x]][1,1])
var2 <- sapply(1:length(wishez), function (x) wishez[[x]][1,1])

hist(var1, breaks = 100)
hist(var2, breaks = 100)

corres <- list()
for(i in 1:niter){
  corres[[i]] <- cov2cor(wishes[[i]])
}

var1 <- sapply(1:length(wishes), function (x) wishes[[x]][1,1])
cor1 <- sapply(1:length(corres), function (x) corres[[x]][1,2])

plot(cor1, log(var1))

#Wishart samples
rwishart <- function(r,R)
{
  X <- rmvnorm(r,sig=R)
  t(X)%*%X
}


#create an MCMC sampler for unknown off diag components

c <- matrix(c(1,.5,.7,.5,1,.6,.7,.6,1), nrow = 3)

data <- rmvnorm(n = 100, sigma = c)

sub <- c[1:2, 1:2]
iter <- 1000
matrices <- list()
likelihoods <- list()
priors <- list()
matrices[[1]] <- diag(3)
subDF <- 10
subP <- dwish(W = matrices[[1]][1:2, 1:2], v = subDF, S = sub * subDF)
wholeDF <- 3
cP <- dwish(W = matrices[[1]], v = wholeDF, S = wholeDF * diag(3))
priors[[1]] <- subL * cL
for(i in 2:iter){
  current <- matrices[[i-1]]
  proposed <- cov2cor(rwish(v = 10, S = current))
  currentL <- 
}
